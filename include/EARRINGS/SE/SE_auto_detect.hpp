#pragma once
#include <future>
#include <streambuf>
#include <biovoltron/algo/align/tailor/tailor.hpp>
#include <biovoltron/algo/align/tailor/bidirectional_index.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/file_io/fastq.hpp>
#include <experimental/vector>
#include <EARRINGS/graph.hpp>
#include <EARRINGS/common.hpp>
#include <EARRINGS/assemble_adapters.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <range/v3/all.hpp>
#include <omp.h>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <memory>
#include <mutex>
#include <optional>
#include <atomic>
#include <utility>
#include <vector>
#include <span>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>

using namespace EARRINGS;
namespace EARRINGS {

// Reads whose exact-match phase Tailor::search_batch interleaves at a time:
// while one read waits for index memory, the next one is extended. Results are
// the same as searching the reads one by one; only the speed depends on it.
inline constexpr auto search_window = size_t{16};

auto make_input_view(std::istream &ifs) -> ranges::any_view<biovoltron::FastqRecord<>> {
    if (is_fastq) {
        return ranges::istream_range<biovoltron::FastqRecord<>>(ifs);
    } else {
        return ranges::istream_range<biovoltron::FastaRecord<>>(ifs)
            | ranges::views::transform([](const auto &rec) {
                return biovoltron::FastqRecord<>{rec.name, rec.seq, {}};
            });
    }
}

// Align reads from the start of the input until 3000 tails from distinct
// reads are found (or the input ends), and return the most common head length
// and the tails. Duplicate reads (e.g. adapter dimers) add their tail once, so
// a few abundant sequences cannot dominate the assembly. Each thread takes the
// next chunk of reads from a shared reader, so no thread waits for the slowest
// read of a batch. Results are then cut back, in input order, to the shortest
// prefix of reads that holds 3000 such tails, so they do not depend on the
// thread count or timing.
template<class IFStream, class Tailor>
std::pair<size_t, std::vector<std::string>> tailor_pipeline(
    IFStream &&ifs, Tailor &&tailor
) {
    constexpr auto npos = std::numeric_limits<uint32_t>::max();
    constexpr auto chunk_size = size_t{128};
    constexpr auto target_tails = size_t{3000};

    struct Hit {
        size_t read_idx;
        size_t head_pos;
        std::optional<std::string> tail;
        std::string seq;  // aligned read sequence, kept only for reads with a tail
    };
    std::unordered_set<std::string> seen_reads;
    std::mutex seen_mutex;
    std::vector<std::vector<Hit>> thread_hits(thread_num);

    auto records_view = make_input_view(ifs);
    auto it = records_view.begin();
    auto next_read_idx = size_t{};
    std::mutex reader_mutex;
    std::atomic<size_t> num_tails{0};

    #pragma omp parallel num_threads(thread_num)
    {
        auto& hits = thread_hits[omp_get_thread_num()];
        std::vector<biovoltron::FastqRecord<>> chunk;
        chunk.reserve(chunk_size);

        while (num_tails.load(std::memory_order_relaxed) < target_tails) {
            chunk.clear();
            auto first_idx = size_t{};
            {
                std::lock_guard lock(reader_mutex);
                first_idx = next_read_idx;
                for (; chunk.size() < chunk_size && it != records_view.end();
                     ++it, ++next_read_idx) {
                    chunk.emplace_back(*it);
                }
            }
            if (chunk.empty()) break;

            // The whole chunk is searched either way (the stop condition is
            // checked between chunks), so searching it as a batch searches the
            // same reads and gives the same alignments as one search() each.
            const auto results = tailor.search_batch(chunk, search_window);
            for (auto i = size_t{}; i < chunk.size(); ++i) {
                const auto& [fwd_aln, rev_aln] = results[i];
                const auto& alignment = fwd_aln.hits.empty() ? rev_aln : fwd_aln;
                if (alignment.hits.empty()) continue;

                auto& hit = hits.emplace_back(Hit{first_idx + i, alignment.head_pos, std::nullopt});
                if (alignment.tail_pos != npos) {
                    hit.tail = alignment.seq.substr(alignment.tail_pos);
                    hit.seq = alignment.seq;
                    auto is_new = false;
                    {
                        std::lock_guard lock(seen_mutex);
                        is_new = seen_reads.insert(hit.seq).second;
                    }
                    if (is_new) num_tails.fetch_add(1, std::memory_order_relaxed);
                }
            }
        }
    }

    // Every chunk handed out was finished, so the hits cover a contiguous
    // prefix of the input; keep the part up to the read that completes
    // target_tails tails, taking each distinct read's tail once.
    std::vector<Hit> all_hits;
    for (auto& hits : thread_hits) std::ranges::move(hits, std::back_inserter(all_hits));
    std::ranges::sort(all_hits, {}, &Hit::read_idx);

    std::unordered_map<size_t, size_t> head_lens;
    std::vector<std::string> tails;
    tails.reserve(std::min(target_tails, num_tails.load()));
    std::unordered_set<std::string> taken_reads;
    for (auto& hit : all_hits) {
        if (tails.size() >= target_tails) break;
        if (hit.head_pos != npos) ++head_lens[hit.head_pos];
        if (hit.tail && taken_reads.insert(std::move(hit.seq)).second)
            tails.push_back(std::move(*hit.tail));
    }

    const auto head_len = head_lens.empty()
        ? size_t{}
        : ranges::max_element(head_lens, {}, &std::pair<const size_t, size_t>::second)->first;

    return {head_len, tails};
}

// Open the SE reads (plain or gzip) as an istream.
inline std::unique_ptr<std::istream> open_se_reads(const std::string& reads_path) {
    if (is_gz_input) {
        auto ifs = std::make_unique<boost::iostreams::filtering_istream>();
        ifs->push(boost::iostreams::gzip_decompressor());
        ifs->push(boost::iostreams::file_source(reads_path, std::ios_base::binary));
        if (!ifs->good()) throw std::runtime_error("Can't open input gz file: " + reads_path);
        return ifs;
    }
    auto ifs = std::make_unique<std::ifstream>(reads_path);
    if (!ifs->good()) throw std::runtime_error("Can't open input file: " + reads_path);
    return ifs;
}

// Run one independent Tailor search per read and write a SAM file (<output>.sam).
// The reads are the input (untrimmed) reads: bases skipped by -A that the head
// extension cannot cover are a leading soft-clip, and the trailing soft-clip is the 3' tail plus any 3' tags and
// adapter (sides swapped on '-' strand records). One record per reported locus
// (NH tag = locus count), unmapped reads included. Records are in input order.
// Opt-in via --sam; cap is the mapped-read limit for this pass (0 = unlimited,
// i.e. the whole file), independent of adapter detection's own preview window.
//
// Known limitation (not a SAM-format violation, left as-is): MAPQ is a crude
// 255 - tail_len and is not lowered for multi-mappers.
template <class Tailor, class Index>
void write_alignment_sam(const std::string& reads_path, const Tailor& tailor, const Index& bidir_index, std::size_t cap) {
    const auto sam_path = output_prefix() + ".sam";
    std::ofstream ofs(sam_path);
    if (!ofs) throw std::runtime_error("Cannot open " + sam_path + " for writing");

    ofs << "@HD\tVN:1.6\tSO:unsorted\n";
    // The index stores the sequence IDs (the FASTA header up to the first whitespace), which is
    // what @SQ SN and RNAME need.
    for (const auto& cb : bidir_index.chr_bounds)
        ofs << "@SQ\tSN:" << cb.chrom << "\tLN:" << bidir_index.get_chr_size(cb.chrom) << '\n';
    ofs << "@PG\tID:EARRINGS\tPN:EARRINGS\tVN:" << std::string(GET_EARRINGS_VERSION) << '\n';

    // The library's Record operator<< appends a trailing tab after every field;
    // strip it (and any stray newline) so the line is spec-clean.
    auto emit = [](std::ostringstream& oss, const auto& sam) {
        std::ostringstream tmp;
        tmp << sam;
        auto s = tmp.str();
        while (!s.empty() && (s.back() == '\t' || s.back() == '\n')) s.pop_back();
        oss << s << '\n';
    };
    // QNAME must not contain whitespace, but a read name can carry a comment after a space.
    auto first_token = [](std::string& s) {
        if (auto p = s.find_first_of(" \t"); p != std::string::npos) s.resize(p);
    };

    auto rfs = open_se_reads(reads_path);
    auto records_view = make_input_view(*rfs);
    auto it = records_view.begin();

    constexpr auto chunk_size = size_t{16};
    const auto batch_size = chunk_size * thread_num;
    size_t n_reads = 0, n_mapped = 0, n_unmapped = 0;

    while (it != records_view.end() && (cap == 0 || n_mapped < cap)) {
        std::vector<biovoltron::FastqRecord<>> records;
        records.reserve(batch_size);
        for (size_t i = 0; i < batch_size && it != records_view.end(); ++i, ++it)
            records.emplace_back(*it);

        std::vector<std::string> lines(records.size());
        std::vector<char> mapped(records.size(), 0);

        // Each thread takes chunk_size reads at a time (as the per-read
        // schedule(dynamic, chunk_size) did) and searches them as one batch;
        // lines stay in input order.
        const auto n_chunks = (records.size() + chunk_size - 1) / chunk_size;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(thread_num)
        for (size_t c = 0; c < n_chunks; ++c) {
            const auto first = c * chunk_size;
            const auto count = std::min(chunk_size, records.size() - first);
            for (size_t i = first; i < first + count; ++i) {
                first_token(records[i].name);
                if (records[i].qual.empty()) records[i].qual.assign(records[i].seq.size(), 'I');
            }
            const auto results = tailor.search_batch(
                std::span<const biovoltron::FastqRecord<>>{records.data() + first, count},
                search_window);

            for (size_t k = 0; k < count; ++k) {
                const auto i = first + k;
                const auto& rec = records[i];
                const auto& [fwd_aln, rev_aln] = results[k];
                auto fwd_sams = biovoltron::aln_to_sam_list(fwd_aln);
                auto rev_sams = biovoltron::aln_to_sam_list(rev_aln);

                std::ostringstream oss;
                if (fwd_sams.empty() && rev_sams.empty()) {
                    oss << rec.name << "\t4\t*\t0\t0\t*\t*\t0\t0\t"
                        << rec.seq << '\t' << rec.qual << '\n';
                } else {
                    mapped[i] = 1;
                    for (const auto& sam : fwd_sams) emit(oss, sam);
                    for (const auto& sam : rev_sams) emit(oss, sam);
                }
                lines[i] = oss.str();
            }
        }

        for (size_t i = 0; i < lines.size(); ++i) {
            ofs << lines[i];
            ++n_reads;
            if (mapped[i]) ++n_mapped; else ++n_unmapped;
        }
    }

    std::cout << "SAM written to " << sam_path << " (" << n_reads << " reads, "
              << n_mapped << " mapped, " << n_unmapped << " unmapped";
    if (cap != 0) std::cout << "; capped at " << cap << " mapped";
    std::cout << ")\n";
}

template <class IndexType>
std::pair<size_t, std::pair<std::string, bool>>
seat_adapter_auto_detect_impl(std::string &reads_path, IndexType &bidir_index) {
    // Hard ceiling on multi-mapping candidates, applied even for -M 0
    // (deliberately unlimited), to avoid exhausting memory on a degenerate read.
    // Tighter than Biovoltron's own 1,000,000-position last-resort ceiling in
    // Tailor::raws2alignment, so that one should never trigger for EARRINGS.
    constexpr size_t HARD_MAX_MULTI_CEILING = 100000;

    biovoltron::Tailor tailor{bidir_index};
    tailor.seed_len = seed_len;
    tailor.allow_seed_mismatch = !no_mismatch;
    tailor.max_multi = std::min(
        min_multi == 0 ? std::numeric_limits<size_t>::max() : min_multi,
        HARD_MAX_MULTI_CEILING
    );
    tailor.skipped_5prime_len = skipped_5prime_len;
    tailor.fix_ext_len = fix_ext_len;
    tailor.min_head_pos = tags5_total_len;

    size_t head_len;
    std::vector<std::string> tails;
    if (is_gz_input) {
        boost::iostreams::filtering_istream ifs;

        ifs.push(boost::iostreams::gzip_decompressor());
        auto src = std::make_shared<boost::iostreams::file_source>(reads_path, std::ios_base::binary);
        if (!src->is_open()) throw std::runtime_error("Can't open input gz file normally\n");

        ifs.push(*src);
        if (!ifs.good()) throw std::runtime_error("Can't open input gz stream normally\n");

        std::tie(head_len, tails) = tailor_pipeline(ifs, tailor);
    } else {
        std::ifstream ifs(reads_path);
        if (!ifs.is_open() || !ifs.good()) throw std::runtime_error("Can't open input file normally\n");

        std::tie(head_len, tails) = tailor_pipeline(ifs, tailor);
    }

    // Optional full per-read alignment output, reusing the loaded index / Tailor.
    if (sam_cap) write_alignment_sam(reads_path, tailor, bidir_index, *sam_cap);

    std::string adapter3;
    std::pair<std::string, bool> adapter3_info;
    if (user_adapter3) {
        adapter3_info = {*user_adapter3, !user_adapter3->empty() && detect_low_complexity(*user_adapter3, init_kmer_size)};
    } else {
        adapter3_info = is_sensitive
            ? assemble_adapters<true>(tails, init_kmer_size, 5)
            : assemble_adapters<false>(tails, init_kmer_size, 3);
    }

    // head_len is either 0 (no 5' extension found in any read) or, thanks to
    // tailor.min_head_pos == tags5_total_len, at least tags5_total_len. Guard
    // the subtraction so the "not found" case reports 0 instead of underflowing.
    const auto adapter5_len = head_len > tags5_total_len ? head_len - tags5_total_len : size_t{};
    std::cout << "5' adapter length: " << adapter5_len << '\n';
    if (head_len == 0 && tags5_total_len > 0) {
        std::cout << "[WARN] declared 5' tag structure but no 5' region was detected; "
                  << "5' tags will not be extracted\n";
    }

    adapter3 = std::get<0>(adapter3_info);

    // Default: assume the whole declared 3' tag region survives on the read
    // after skewer (peeled per-read by trim_tags3). Overridden below when the
    // tails show that some tag bases were absorbed into adapter3.
    tags3_kept_on_read = tags3_total_len;
    tags3_absorbed_prefix.clear();

    if (user_adapter3) {
        // Taken as given, like the declared tags: nothing is estimated from the tails.
        if (adapter3.empty()) std::cout << "no 3' adapter (-a \"\"), 3' adapter trimming is off\n";
        else std::cout << "3' adapter given: " << adapter3 << '\n';
    } else if (adapter3 == "" && trim_default_adapter3) {
        std::cout << "unable to detect adapter, use default adapter\n";
        adapter3 = DEFAULT_ADAPTER1;
    } else if (adapter3 == "") {
        std::cout << "[WARN] unable to detect adapter, 3' adapter trimming is off (use -a to give the adapter)\n";
        // The tags lie before an adapter that was not found, so their position is unknown.
        if (!tag_structure3.empty()) {
            std::cout << "[WARN] 3' tags will not be extracted\n";
            tags3_kept_on_read = 0;
        }
    } else {
        // is low complexity
        adapter3 = std::get<1>(adapter3_info)
            ? adapter3.substr(0, 16)
            : adapter3.substr(0, 32);

        std::cout << "3' adapter found: " << adapter3 << '\n';

        // The estimate below is only meaningful when adapter3 was actually
        // assembled from this data; a given or default adapter (above) has
        // no relationship to these reads' tails.
        if (!tag_structure3.empty()) {
            const auto [est_len, reliable] = estimate_tags3_len(tails, adapter3);
            estimated_tags3_len = est_len;
            std::cout << "estimated tags3 length: " << estimated_tags3_len << '\n';

            // estimated_tags3_len = how many tag bases sit OUTSIDE adapter3
            // (still on the read after skewer). The remaining innermost bases
            // were folded into adapter3 by the de Bruijn assembly because they
            // are ~constant across reads (a variable tag like a UMI never
            // assembles in).
            constexpr size_t MIN_REMAINING_ADAPTER_LEN = 10;
            if (!reliable) {
                // Signal too weak/flat to locate the insert/tag boundary
                // (e.g. reads whose own 3' tails spread the offsets). The user
                // declared the structure, so take it as given: all declared
                // tag bases sit on the read between the insert and adapter3.
                // Dropping the tags instead would leave them on every insert,
                // where they read as a tail.
                std::cout << "[WARN] 3' tag position could not be estimated reliably "
                          << "(estimate " << estimated_tags3_len << " nt); using the declared "
                          << tags3_total_len << " nt\n";
                tags3_kept_on_read = tags3_total_len;
                tags3_absorbed_prefix.clear();
            } else if (estimated_tags3_len >= tags3_total_len) {
                // Nothing absorbed: the whole tag region is on the read.
                tags3_kept_on_read = tags3_total_len;
                tags3_absorbed_prefix.clear();
            } else {
                const auto swallowed_len = tags3_total_len - estimated_tags3_len;

                // Strip the swallowed bases off adapter3 - so skewer removes only
                // the genuine adapter and leaves the whole tag region on the read
                // for per-read extraction - ONLY when there is a variable portion
                // that must come from the read AND enough real adapter remains for
                // skewer to anchor on. Otherwise keep adapter3 fused: skewer then
                // anchors on the longer sequence (a more precise insert cut) and
                // the absorbed bases, being ~constant by definition, are rebuilt
                // from adapter3's prefix.
                if (estimated_tags3_len > 0 && adapter3.size() >= swallowed_len + MIN_REMAINING_ADAPTER_LEN) {
                    std::cout << "3' tags partially absorbed into the assembled adapter; "
                              << "stripping " << swallowed_len << " leading base(s) from it\n";
                    adapter3 = adapter3.substr(swallowed_len);
                    tags3_kept_on_read = tags3_total_len;
                    tags3_absorbed_prefix.clear();
                } else {
                    tags3_kept_on_read = estimated_tags3_len;
                    tags3_absorbed_prefix = adapter3.substr(0, std::min(swallowed_len, adapter3.size()));
                    std::cout << "3' tags absorbed into the assembled adapter; "
                              << "reconstructing " << tags3_absorbed_prefix.size()
                              << " constant base(s) ('" << tags3_absorbed_prefix
                              << "') from it\n";
                }
            }
        }
    }

    std::get<0>(adapter3_info) = adapter3;

    return {head_len, adapter3_info};
}

// Frees the index of the last detection. Returning gigabytes to the OS takes
// about a second for a genome index; on this thread it overlaps the trimming
// instead of delaying it. Joined at exit.
inline std::jthread index_release;

template <class IndexType>
std::pair<size_t, std::pair<std::string, bool>>
seat_adapter_auto_detect_with(std::string &reads_path, std::ifstream &fm_ifs, const std::string &fm_path) {
    auto bidir_index = std::make_unique<IndexType>(INDEX_LOOKUP_LEN);
    // With the path, the index's large arrays are read by several threads.
    // The lookup table is skipped (2 GB in an index built with length 14):
    // only get_range(seed) uses it, and detection searches through a
    // BidirectionalSearcher. Without it the file's lookup length need not
    // match INDEX_LOOKUP_LEN.
    bidir_index->load(fm_ifs, fm_path, /*with_lookup=*/false);
    auto result = seat_adapter_auto_detect_impl(reads_path, *bidir_index);
    index_release = std::jthread([index = std::move(bidir_index)]() mutable { index.reset(); });
    return result;
}

// .table files start with a 1-byte flag (written by `build` in EARRINGS.cpp)
// saying whether the index is Dense (SA_INTV bucket <= 100Kbp) or Sampled,
// then a size_flag byte saying whether it's the 32-bit or 64-bit
// instantiation -- together they determine which compile-time
// BidirectionalIndex instantiation load() must be called on.
std::pair<size_t, std::pair<std::string, bool>> seat_adapter_auto_detect(std::string &reads_path) {
    try {
        const auto fm_path = index_prefix + ".table";
        std::ifstream fm_ifs{fm_path};

        char dense_flag = 0;
        char size_flag = 0;
        fm_ifs.read(&dense_flag, 1);
        fm_ifs.read(&size_flag, 1);
        if (size_flag) {
            if (dense_flag) return seat_adapter_auto_detect_with<DenseIndex64>(reads_path, fm_ifs, fm_path);
            return seat_adapter_auto_detect_with<SampledIndex64>(reads_path, fm_ifs, fm_path);
        }
        if (dense_flag) return seat_adapter_auto_detect_with<DenseIndex32>(reads_path, fm_ifs, fm_path);
        return seat_adapter_auto_detect_with<SampledIndex32>(reads_path, fm_ifs, fm_path);
    } catch (const std::exception& e) {
        // e.g. an index built by an older version, or an unreadable input
        std::cerr << "Error: " << e.what() << "\n";
        std::exit(EXIT_FAILURE);
    }
}

// Appends what `oss << rec << '\n'` would write, without a stream: the record in FASTA or
// FASTQ format and a newline. One of these per read is hot when the reads go through a pipe.
template <class Record>
void append_record(std::string& out, const Record& rec) {
    out += Record::START_SYMBOL;
    out += rec.name;
    out += '\n';
    out += rec.seq;
    if constexpr (std::is_same_v<Record, biovoltron::FastqRecord<>>) {
        out += '\n';
        out += Record::DELIM;
        out += '\n';
        out += rec.qual;
    }
    out += '\n';
}

// Reads `head`, then whatever is left in `src`.
class ChainedStreambuf : public std::streambuf {
  public:
    ChainedStreambuf(std::string head, std::streambuf* src) : head_(std::move(head)), src_(src), buf_(1 << 16) {
        setg(head_.data(), head_.data(), head_.data() + head_.size());
    }

  protected:
    int_type underflow() override {
        if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
        if (traits_type::eq_int_type(src_->sgetc(), traits_type::eof())) return traits_type::eof();
        const auto avail = std::max<std::streamsize>(1, src_->in_avail());
        const auto n = src_->sgetn(buf_.data(), std::min<std::streamsize>(avail, buf_.size()));
        if (n <= 0) return traits_type::eof();
        setg(buf_.data(), buf_.data(), buf_.data() + n);
        return traits_type::to_int_type(*gptr());
    }

  private:
    std::string head_;
    std::streambuf* src_;
    std::vector<char> buf_;
};

// Fast path of trim_heads_to_stream for records in the usual layout: 4-line FASTQ
// ("@id\nseq\n+...\nqual\n", qual as long as seq) or 2-line FASTA. Reads large chunks, splits
// them into records, cuts the 5' ends on several threads and writes them in order, while the next
// chunk is prepared. Stops at the
// first record of another layout (or the unterminated last one) and returns the unread bytes
// from there on, for the general parser to continue with what is left in `src`.
inline std::string trim_heads_fast(std::streambuf* src, size_t head_len, size_t start_pos, bool keep_tags,
                                   std::vector<std::vector<std::string>>& tags5, std::vector<std::string>& names5,
                                   FILE* wfp, size_t chunk_bytes = size_t{4} << 20) {
    const bool fastq = is_fastq;
    const char start_symbol = fastq ? '@' : '>';
    const size_t lines_per_record = fastq ? 4 : 2;
    const auto num_tags5 = tag_structure5.size();
    const auto umi_idx = umi_index(tag_structure5);
    const bool add_umi = record_tags && num_tags5 > 0 && umi_idx != num_tags5;
    // Beyond this length every tag's start is inside the read (an empty or shortened tag at the
    // very end included), as std::string::substr requires; the general parser reports the rest.
    size_t min_len = 0;
    if (keep_tags && num_tags5 > 0) {
        min_len = start_pos;
        for (size_t i = 0; i + 1 < num_tags5; ++i) min_len += tag_structure5[i].second;
    }
    size_t umi_off = start_pos;
    for (size_t i = 0; i < umi_idx && i < num_tags5; ++i) umi_off += tag_structure5[i].second;

    struct Rec { const char *name, *seq, *qual; uint32_t name_len, seq_len; };
    // skewer runs its own threads next to this; a few are enough to keep up with it.
    const size_t num_threads = std::clamp<size_t>(thread_num, 1, 4);
    std::string buf;
    std::vector<Rec> recs;
    std::vector<std::string> outs(num_threads);
    size_t carry = 0;
    bool eof = false;

    // One chunk is sent while the next is prepared; get() waits for it and reports its failure.
    std::vector<std::string> sending(num_threads);
    std::future<void> sent;
    auto send = [&] {
        for (const auto& o : sending) {
            // Fails (EPIPE) once skewer has stopped reading, e.g. when it
            // could not open its output.
            if (!o.empty() && fwrite(o.data(), 1, o.size(), wfp) != o.size())
                throw std::runtime_error("Cannot pass reads to skewer");
        }
    };

    for (;;) {
        // Fill the chunk, keeping what was read if the source fails (a corrupt gzip file).
        buf.resize(carry + chunk_bytes);
        size_t total = carry;
        try {
            while (total < buf.size()) {
                if (src->sgetc() == std::char_traits<char>::eof()) { eof = true; break; }
                const auto avail = std::max<std::streamsize>(1, src->in_avail());
                const auto n = src->sgetn(buf.data() + total, std::min<std::streamsize>(avail, buf.size() - total));
                if (n <= 0) { eof = true; break; }
                total += static_cast<size_t>(n);
            }
        } catch (...) {
            eof = true;   // the general parser meets the same error
        }
        buf.resize(total);

        // Split the complete records; what follows one must start the next (or be the end).
        recs.clear();
        size_t used = 0;
        bool regular = true;
        while (used < total) {
            const char* base = buf.data();
            size_t line_start[4], line_len[4];
            auto p = used;
            bool complete = true;
            for (size_t l = 0; l < lines_per_record; ++l) {
                const auto* nl = p < total ? static_cast<const char*>(memchr(base + p, '\n', total - p)) : nullptr;
                if (!nl) { complete = false; break; }
                line_start[l] = p;
                line_len[l] = static_cast<size_t>(nl - (base + p));
                p = static_cast<size_t>(nl - base) + 1;
            }
            if (!complete || (p == total && !eof)) break;
            if (base[used] != start_symbol || (p < total && base[p] != start_symbol) || line_len[1] < min_len ||
                (fastq && (base[line_start[2]] != '+' || line_len[3] != line_len[1]))) {
                regular = false;
                break;
            }
            recs.push_back({base + line_start[0] + 1, base + line_start[1], fastq ? base + line_start[3] : nullptr,
                            static_cast<uint32_t>(line_len[0] - 1), static_cast<uint32_t>(line_len[1])});
            used = p;
        }

        const auto n = recs.size();
        const auto first = names5.size();
        if (keep_tags) {
            names5.resize(first + n);
            tags5.resize(first + n);
        }
        for (auto& o : outs) o.clear();
        #pragma omp parallel num_threads(num_threads)
        {
            const auto tid = static_cast<size_t>(omp_get_thread_num());
            const auto nt = static_cast<size_t>(omp_get_num_threads());
            auto& o = outs[tid];
            for (auto idx = n * tid / nt; idx < n * (tid + 1) / nt; ++idx) {
                const auto& r = recs[idx];
                const std::string_view name(r.name, r.name_len), seq(r.seq, r.seq_len);
                const auto id_end = std::min(name.find(' '), name.size());
                if (keep_tags) {
                    names5[first + idx].assign(name.substr(0, id_end));
                    auto& rec_tags = tags5[first + idx];
                    rec_tags.resize(num_tags5);
                    for (size_t i = 0, pos = start_pos; i < num_tags5; ++i) {
                        rec_tags[i].assign(seq.substr(pos, tag_structure5[i].second));
                        pos += tag_structure5[i].second;
                    }
                }
                o += start_symbol;
                o.append(name.substr(0, id_end));
                if (add_umi) {
                    o += ':';
                    o.append(seq.substr(umi_off, tag_structure5[umi_idx].second));
                }
                o.append(name.substr(id_end));
                o += '\n';
                const auto cut = std::min<size_t>(head_len, seq.size());
                o.append(seq.substr(cut));
                if (fastq) {
                    o += "\n+\n";
                    o.append(r.qual + cut, seq.size() - cut);
                }
                o += '\n';
            }
        }
        if (sent.valid()) sent.get();
        std::swap(outs, sending);
        sent = std::async(std::launch::async, send);

        buf.erase(0, used);
        if (!regular || eof) {
            sent.get();
            return buf;
        }
        carry = buf.size();
    }
}

void trim_heads_to_stream(std::string& reads_path, size_t head_len, std::vector<std::vector<std::string>>& tags5, std::vector<std::string>& names5, FILE* wfp) {
    std::unique_ptr<std::istream> ifs;
    if (is_gz_input) {
        auto gz_ifs = std::make_unique<boost::iostreams::filtering_istream>();
        auto src = std::make_shared<boost::iostreams::file_source>(reads_path, std::ios_base::binary);
        if (!src->is_open()) throw std::runtime_error("Can't open input gz file");
        gz_ifs->push(boost::iostreams::gzip_decompressor());
        gz_ifs->push(*src);
        ifs = std::move(gz_ifs);
    } else {
        auto* fs = new std::ifstream(reads_path);
        if (!fs->is_open() || !fs->good()) throw std::runtime_error("Can't open input file");
        ifs.reset(fs);
    }

    const auto num_tags5 = tag_structure5.size();
    // names5 and tags5 only feed the tags table (see export_tags_tsv in main), which is
    // written only with --record_tags and a declared tag structure; otherwise keeping a name
    // and a vector for every read is wasted time and memory (the 5' tags are still cut off).
    const bool keep_tags = record_tags && (!tag_structure5.empty() || !tag_structure3.empty());
    const auto start_pos = head_len - tags5_total_len;
    const auto umi_idx = umi_index(tag_structure5);
    if (record_tags && num_tags5 > 0 && umi_idx != num_tags5) umi5_in_read_id = true;

    // Whatever the fast path did not take, with any layout the Biovoltron parser accepts.
    ChainedStreambuf rest_buf(trim_heads_fast(ifs->rdbuf(), head_len, start_pos, keep_tags, tags5, names5, wfp), ifs->rdbuf());
    std::istream rest(&rest_buf);

    auto extract_tags5 = [&]<typename Record>() {
        Record rec;
        std::string out;
        while (rest >> rec) {
            std::vector<std::string> rec_tags;
            if (keep_tags && num_tags5 > 0) {
                rec_tags.resize(num_tags5);
                auto pos = start_pos;
                for (auto i = size_t{}; i < num_tags5; ++i) {
                    rec_tags[i] = rec.seq.substr(pos, tag_structure5[i].second);
                    pos += tag_structure5[i].second;
                }
            }
            // The original read ID keys the tag table; export_tags_tsv adds
            // the UMI back to match the names trim_tags3 reads from skewer.
            if (keep_tags) names5.push_back(read_qname(rec.name));
            if (record_tags && num_tags5 > 0 && umi_idx != num_tags5) add_umi_to_read_id(rec.name, rec_tags[umi_idx], true);
            if (keep_tags) tags5.push_back(std::move(rec_tags));

            rec.seq.erase(0, head_len);
            if constexpr (std::is_same_v<Record, biovoltron::FastqRecord<>>) {
                rec.qual.erase(0, head_len);
            }

            out.clear();
            append_record(out, rec);
            // Fails (EPIPE) once skewer has stopped reading, e.g. when it
            // could not open its output.
            if (fwrite(out.data(), 1, out.size(), wfp) != out.size())
                throw std::runtime_error("Cannot pass reads to skewer");
        }
    };
    if (is_fastq) {
        extract_tags5.template operator()<biovoltron::FastqRecord<>>();
    } else {
        extract_tags5.template operator()<biovoltron::FastaRecord<>>();
    }
}

}
