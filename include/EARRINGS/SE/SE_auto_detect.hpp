#pragma once
#include <biovoltron/algo/align/tailor/tailor.hpp>
#include <biovoltron/algo/align/tailor/index.hpp>
#include <experimental/vector>
#include <EARRINGS/graph.hpp>
#include <EARRINGS/common.hpp>
#include <EARRINGS/assemble_adapters.hpp>
#include <EARRINGS/SE/format_reader.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <range/v3/all.hpp>
#include <string>
#include <cmath>
#include <fstream>
#include <atomic>
#include <utility>

using namespace EARRINGS;
namespace EARRINGS
{
    auto make_input_view(std::istream &ifs) -> ranges::any_view<biovoltron::FastqRecord<>> {
        if (is_fastq) {
            constexpr EARRINGS::format_reader_fn<biovoltron::FastqRecord<>> fastq_reader{};
            return ifs | fastq_reader();
        } else {
            constexpr EARRINGS::format_reader_fn<biovoltron::FastaRecord<>> fasta_reader{};
            return ifs
                   | fasta_reader()
                   | ranges::view::transform([](const auto &rec) {
                return biovoltron::FastqRecord<>{rec.name, rec.seq, std::string(rec.seq.size(), 'I')};
            });
        }
    }

    template<class IFStream, class Tailor>
    std::vector<std::string> tailor_pipeline(IFStream &&ifs, size_t thread_num, Tailor &&tailor, size_t num_reads) {
        std::vector<std::string> tails;
        tails.reserve(num_reads);
        size_t counter(0);
        while (ifs.good() && tails.size() < 3000 && counter < 3) {
            auto alignment_pairs = make_input_view(ifs)
                | ranges::view::transform([&tailor](auto &&query) {
                    return tailor.search(
                        is_fastq
                            ? static_cast<biovoltron::FastqRecord<>>(query)
                            : biovoltron::FastqRecord<>{query.name, query.seq, std::string(query.seq.size(), 'I')}
                    );
                })
                | ranges::to_vector;

            for (auto &&pair : alignment_pairs) {
                const auto &alignment = pair.first.hits.empty() ? pair.second : pair.first;
                if (!alignment.hits.empty()) {
                    tails.emplace_back(alignment.seq.substr(alignment.seq.length() - alignment.tail_pos - 1));
                }
            }
            counter++;
        }

        return tails;
    }


    std::pair<std::string, bool> seat_adapter_auto_detect(std::string &reads_path, size_t thread_num = 1) {
        std::vector<std::string> tails;
        biovoltron::Index fm_index{12}, rc_fm_index{12};
        std::ifstream fm_ifs{index_prefix + ".table"}, rc_fm_ifs{index_prefix + ".rc_table"};
        fm_index.load(fm_ifs);
        rc_fm_index.load(rc_fm_ifs);
        biovoltron::Tailor tailor{fm_index, rc_fm_index};
        tailor.seed_len = seed_len;
        tailor.allow_seed_mismatch = !no_mismatch;

        if (is_gz_input) {
            boost::iostreams::filtering_istream ifs;

            ifs.push(boost::iostreams::gzip_decompressor());
            auto&& src(boost::iostreams::file_source(reads_path, std::ios_base::binary));
            if (!src.is_open())
                throw std::runtime_error("Can't open input gz file normally\n");

            ifs.push(src);
            if (!ifs.good())
                throw std::runtime_error("Can't open input gz stream normally\n");

            tails = tailor_pipeline(ifs, thread_num, tailor, DETECT_N_READS);
        } else {
            std::ifstream ifs(reads_path);
            if (!(ifs.is_open() && ifs.good()))
                throw std::runtime_error("Can't open input file normally\n");

            tails = tailor_pipeline(ifs, thread_num, tailor, DETECT_N_READS);
        }
        // std::cerr << "total number of tails sampled: " << tails.size() << "\n";

        std::string adapter;
        std::pair<std::string, bool> adapter_info;
        if (is_sensitive) {
            adapter_info = assemble_adapters<true>(tails, init_kmer_size, 5);
        } else {
            adapter_info = assemble_adapters<false>(tails, init_kmer_size, 3);
        }
        adapter = std::get<0>(adapter_info);

        if (adapter == "") {
            std::cout << "unable to detect adapter, use default adapter\n";
            adapter = DEFAULT_ADAPTER1;
        } else {
            // is low complexity
            if (std::get<1>(adapter_info))
            {
                adapter = adapter.substr(0, 16);
            }
            else
            {
                adapter = adapter.substr(0, 32);
            }

            std::cout << "adapter found: " << adapter << '\n';
        }

         std::get<0>(adapter_info) = adapter;

         return adapter_info;
    }
}
