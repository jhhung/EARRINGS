#pragma once
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
#include <atomic>
#include <utility>
#include <vector>
#include <unordered_map>
#include <utility>

using namespace EARRINGS;
namespace EARRINGS {

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

template<class IFStream, class Tailor>
std::pair<size_t, std::vector<std::string>> tailor_pipeline(
    IFStream &&ifs, Tailor &&tailor, size_t num_reads
) {
    constexpr auto npos = std::numeric_limits<uint32_t>::max(); 
    constexpr auto target_tails = size_t{3000};

    std::unordered_map<size_t, size_t> head_lens;
    std::vector<std::string> tails;
    tails.reserve(target_tails);

    constexpr auto chunk_size = size_t{16};
    const auto batch_size = chunk_size * thread_num;

    auto records_view = make_input_view(ifs);
    auto it = records_view.begin();

    while (ifs.good() && tails.size() < target_tails) {
        std::vector<biovoltron::FastqRecord<>> records;
        records.reserve(batch_size);
        for (auto i = size_t{}; i < batch_size && it != records_view.end(); ++i, ++it) {
            records.emplace_back(*it);
        }

        #pragma omp parallel for schedule(dynamic, chunk_size) num_threads(thread_num)
        for (auto&& record : records) {
            auto aln_pair = tailor.search(record);
            const auto& [fwd_aln, rev_aln] = aln_pair;
            const auto& alignment = fwd_aln.hits.empty() ? rev_aln : fwd_aln;

            if (!alignment.hits.empty()) {
                if (alignment.head_pos != npos) {
                    #pragma omp critical (head)
                    ++head_lens[alignment.head_pos];
                }

                if (alignment.tail_pos != npos) {
                    auto tail = alignment.seq.substr(alignment.tail_pos);
                    #pragma omp critical (tail)
                    tails.push_back(std::move(tail));
                }
            }
        }
    }

    const auto head_len = head_lens.empty()
        ? size_t{}
        : ranges::max_element(head_lens, {}, &std::pair<const size_t, size_t>::second)->first;

    return {head_len, tails};
}

std::pair<size_t, std::pair<std::string, bool>> seat_adapter_auto_detect(std::string &reads_path) {
    biovoltron::BidirectionalIndex<SA_INTV> bidir_index;
    std::ifstream fm_ifs{index_prefix + ".table"}, rev_fm_ifs{index_prefix + ".rev_table"};
    bidir_index.load(fm_ifs, rev_fm_ifs);
    biovoltron::Tailor tailor{bidir_index};
    tailor.seed_len = seed_len;
    tailor.allow_seed_mismatch = !no_mismatch;
    tailor.max_multi = min_multi == 0 ? std::numeric_limits<size_t>::max() : min_multi;
    tailor.skipped_5prime_len = skipped_5prime_len;
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

        std::tie(head_len, tails) = tailor_pipeline(ifs, tailor, DETECT_N_READS);
    } else {
        std::ifstream ifs(reads_path);
        if (!ifs.is_open() || !ifs.good()) throw std::runtime_error("Can't open input file normally\n");

        std::tie(head_len, tails) = tailor_pipeline(ifs, tailor, DETECT_N_READS);
    }

    std::string adapter3;
    std::pair<std::string, bool> adapter3_info;
    adapter3_info = is_sensitive
        ? assemble_adapters<true>(tails, init_kmer_size, 5)
        : assemble_adapters<false>(tails, init_kmer_size, 3);

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

    if (adapter3 == "") {
        std::cout << "unable to detect adapter, use default adapter\n";
        adapter3 = DEFAULT_ADAPTER1;
    } else {
        // is low complexity
        adapter3 = std::get<1>(adapter3_info)
            ? adapter3.substr(0, 16)
            : adapter3.substr(0, 32);

        std::cout << "3' adapter found: " << adapter3 << '\n';

        // The estimate below is only meaningful when adapter3 was actually
        // assembled from this data; DEFAULT_ADAPTER1 (the fallback above) has
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
                // Signal too weak/flat to locate the insert/tag boundary.
                std::cout << "[WARN] 3' tag position could not be estimated reliably "
                          << "(try adjusting --seed_len / --prune_factor); "
                          << "3' tags will not be extracted\n";
                tags3_kept_on_read = 0;
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
    const auto start_pos = head_len - tags5_total_len;
    const size_t umi_idx = ranges::find_if(tag_structure5, [](const auto& p) { return p.first == "UMI"; }) - tag_structure5.begin();

    auto extract_tags5 = [&]<typename Record>() {
        Record rec;
        while (*ifs >> rec) {
            std::vector<std::string> rec_tags;
            if (num_tags5 > 0) {
                rec_tags.resize(num_tags5);
                auto pos = start_pos;
                for (auto i = size_t{}; i < num_tags5; ++i) {
                    rec_tags[i] = rec.seq.substr(pos, tag_structure5[i].second);
                    pos += tag_structure5[i].second;
                }
                if (umi_idx != num_tags5) rec.name += ":" + rec_tags[umi_idx];
            }
            // Record the name as written downstream (post UMI append) so it
            // matches what trim_tags3 reads back from skewer's output.
            names5.push_back(read_qname(rec.name));
            tags5.push_back(std::move(rec_tags));

            rec.seq.erase(0, head_len);
            if constexpr (std::is_same_v<Record, biovoltron::FastqRecord<>>) {
                rec.qual.erase(0, head_len);
            }

            std::ostringstream oss;
            oss << rec << '\n';
            const auto s = oss.str();
            fwrite(s.c_str(), 1, s.size(), wfp);
        }
    };
    if (is_fastq) {
        extract_tags5.template operator()<biovoltron::FastqRecord<>>();
    } else {
        extract_tags5.template operator()<biovoltron::FastaRecord<>>();
    }
}

}
