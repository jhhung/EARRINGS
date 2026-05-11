#pragma once
#include <biovoltron/algo/align/tailor/tailor.hpp>
#include <biovoltron/algo/align/tailor/bidirectional_index.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/file_io/fastq.hpp>
#include <experimental/vector>
#include <EARRINGS/graph.hpp>
#include <EARRINGS/common.hpp>
#include <EARRINGS/assemble_adapters.hpp>
#include <EARRINGS/SE/format_reader.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <range/v3/all.hpp>
#include <omp.h>
#include <string>
#include <cmath>
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
        constexpr EARRINGS::format_reader_fn<biovoltron::FastqRecord<>> fastq_reader{};
        return ifs | fastq_reader();
    } else {
        constexpr EARRINGS::format_reader_fn<biovoltron::FastaRecord<>> fasta_reader{};
        return ifs | fasta_reader() | ranges::view::transform([](const auto &rec) {
            return biovoltron::FastqRecord<>{rec.name, rec.seq, std::string(rec.seq.size(), 'I')};
        });
    }
}

template<class IFStream, class Tailor>
std::pair<size_t, std::vector<std::string>> tailor_pipeline(
    IFStream &&ifs, Tailor &&tailor, size_t num_reads
) {
    std::unordered_map<size_t, size_t> head_lens;
    std::vector<std::string> tails;
    tails.reserve(3000);

    constexpr size_t chunk_size = 16;
    const size_t batch_size = chunk_size * thread_num;

    auto records_view = make_input_view(ifs);
    auto it = records_view.begin();

    while (ifs.good() && tails.size() < 3000) {
        std::vector<biovoltron::FastqRecord<>> records;
        records.reserve(batch_size);
        for (int i = 0; i < batch_size && it != records_view.end(); ++i, ++it) {
            records.emplace_back(*it);
        }

        #pragma omp parallel for schedule(dynamic, chunk_size) num_threads(thread_num)
        for (auto&& record : records) {
            auto aln_pair = tailor.search(record);
            const auto& [fwd_aln, rev_aln] = aln_pair;
            const auto& alignment = fwd_aln.hits.empty() ? rev_aln : fwd_aln;

            if (!alignment.hits.empty()) {
                if (alignment.head_pos != static_cast<uint32_t>(-1)) {
                    #pragma omp critical (head)
                    ++head_lens[alignment.head_pos];
                }

                if (alignment.tail_pos != static_cast<uint32_t>(-1)) {
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
    tailor.max_multi = 1;
    tailor.skipped_5prime_len = skipped_5prime_len;

    size_t head_len;
    std::vector<std::string> tails;
    if (is_gz_input) {
        boost::iostreams::filtering_istream ifs;

        ifs.push(boost::iostreams::gzip_decompressor());
        auto src = std::make_shared<boost::iostreams::file_source>(reads_path, std::ios_base::binary);
        if (!src->is_open())
            throw std::runtime_error("Can't open input gz file normally\n");

        ifs.push(*src);
        if (!ifs.good())
            throw std::runtime_error("Can't open input gz stream normally\n");

        std::tie(head_len, tails) = tailor_pipeline(ifs, tailor, DETECT_N_READS);
    } else {
        std::ifstream ifs(reads_path);
        if (!(ifs.is_open() && ifs.good()))
            throw std::runtime_error("Can't open input file normally\n");

        std::tie(head_len, tails) = tailor_pipeline(ifs, tailor, DETECT_N_READS);
    }

    std::string adapter3;
    std::pair<std::string, bool> adapter3_info;
    if (is_sensitive) {
        adapter3_info = assemble_adapters<true>(tails, init_kmer_size, 5);
    } else {
        adapter3_info = assemble_adapters<false>(tails, init_kmer_size, 3);
    }

    const auto adapter5_len = head_len - tags5_total_len;
    std::cout << "5' adapter length: " << adapter5_len << '\n';

    adapter3 = std::get<0>(adapter3_info);

    if (adapter3 == "") {
        std::cout << "unable to detect adapter, use default adapter\n";
        adapter3 = DEFAULT_ADAPTER1;
    } else {
        // is low complexity
        if (std::get<1>(adapter3_info)) {
            adapter3 = adapter3.substr(0, 16);
        } else {
            adapter3 = adapter3.substr(0, 32);
        }

        std::cout << "3' adapter found: " << adapter3 << '\n';
    }

    std::get<0>(adapter3_info) = adapter3;

    if (!tag_structure3.empty()) estimated_tags3_len = estimate_tags3_len(tails, adapter3);

    return {head_len, adapter3_info};
}

std::string trim_heads_to_tmpfile(std::string& reads_path, size_t head_len, std::vector<std::vector<std::string>>& tags5) {
    const auto ts = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    const std::string tmpfile = "/tmp/EARRINGS_trimmed_" + std::to_string(getpid()) + "_" + std::to_string(ts) + ".tmp";
    
    std::unique_ptr<std::istream> ifs;
    if (is_gz_input) {
        auto gz_ifs = std::make_unique<boost::iostreams::filtering_istream>();
        auto src = std::make_shared<boost::iostreams::file_source>(reads_path, std::ios_base::binary);
        if (!src->is_open()) {
            throw std::runtime_error("Can't open input gz file");
        }
        gz_ifs->push(boost::iostreams::gzip_decompressor());
        gz_ifs->push(*src);
        ifs = std::move(gz_ifs);
    } else {
        auto* fs = new std::ifstream(reads_path);
        if (!fs->is_open() || !fs->good()) {
            throw std::runtime_error("Can't open input file");
        }
        ifs.reset(fs);
    }
    
    std::ofstream ofs(tmpfile);
    const int num_tags5 = tag_structure5.size();
    size_t start_pos = head_len - tags5_total_len;

    const auto umi_idx = ranges::find_if(tag_structure5, [](const auto& p) { return p.first == "UMI"; }) - tag_structure5.begin();

    auto extract_tags5 = [&]<typename Record>() {
        constexpr EARRINGS::format_reader_fn<Record> reader{};
        std::vector<Record> records;
        for (auto&& rec : (*ifs) | reader()) {
            records.push_back(std::move(rec));
        }
        const size_t num_records = records.size();
        tags5.resize(num_records);

        #pragma omp parallel
        {
            std::ostringstream oss;

            #pragma omp for nowait
            for (size_t idx = 0; idx < num_records; ++idx) {
                auto& rec = records[idx];
                auto& rec_tags = tags5[idx];

                if (num_tags5 > 0) {
                    rec_tags.resize(num_tags5);
                    size_t pos = start_pos;
                    for (size_t i = 0; i < num_tags5; ++i) {
                        const auto& len = tag_structure5[i].second;
                        rec_tags[i] = rec.seq.substr(pos, len);
                        pos += len;
                    }
                    if (umi_idx != num_tags5) rec.name += ":" + rec_tags[umi_idx];
                }

                rec.seq.erase(0, head_len);
                if constexpr (std::is_same_v<Record, biovoltron::FastqRecord<>>) rec.qual.erase(0, head_len);
                oss << rec << '\n';
            }

            #pragma omp critical
            ofs << oss.str();
        }
    };
    if (is_fastq) {
        extract_tags5.template operator()<biovoltron::FastqRecord<>>();
    } else {
        extract_tags5.template operator()<biovoltron::FastaRecord<>>();
    }

    return tmpfile;
}

}
