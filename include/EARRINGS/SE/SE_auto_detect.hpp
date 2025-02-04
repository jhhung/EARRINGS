#pragma once
#include <string>
#include <biovoltron/algo/align/tailor/tailor.hpp>
#include <biovoltron/algo/align/tailor/index.hpp>
#include <Tailor/range.hpp>
#include <experimental/vector>
#include <EARRINGS/graph.hpp>
#include <EARRINGS/common.hpp>
#include <EARRINGS/assemble_adapters.hpp>
#include <Nucleona/range/v3_impl.hpp>
#include <Nucleona/parallel/thread_pool.hpp>
#include <Nucleona/parallel/asio_pool.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <cmath>
#include <fstream>
#include <atomic>
#include <utility>

using namespace EARRINGS;
namespace EARRINGS
{
template<class FORMAT>
void get_reads(std::istream& is, std::vector<FORMAT>& buf, size_t num_reads)
{
    FORMAT f;
	while(FORMAT::get_obj(is, f))
	{
		buf.emplace_back(f);
		if (buf.size() == num_reads)
			break;
	} 

}

using trueType = std::bool_constant<true>;
using falseType = std::bool_constant<false>;

template<class IFStream, class Tailor>
std::vector<std::string> tailor_pipeline(IFStream&& ifs
                                       , size_t thread_num
                                       , Tailor&& tailor
                                       , size_t num_reads)
{
    auto tp = nucleona::parallel::make_asio_pool(thread_num);
    constexpr pipeline::range::format_reader_fn<biovoltron::FastqRecord<>> format_reader{};

    std::vector<std::string> tails;
    tails.reserve(num_reads);
    size_t counter(0);
    while (ifs.good() && tails.size() < 3000 && counter < 3) {
        ifs
        | format_reader()
        | ranges::view::transform([&tailor](auto &&query) {
            if (is_fastq) {
                return tailor.search(query);
            } else {
                return tailor.search(biovoltron::FastqRecord<>{query.name, query.seq, std::string(query.seq.size(), 'I')});
            }
        })
        | ranges::view::transform([&tails](auto &&alignment) {
            if (!alignment.hits.empty()) {  // alignment.tail_pos >= 0
                tails.emplace_back(alignment.seq.substr(alignment.seq.length() - alignment.tail_pos - 1));
            }
            return alignment;
        })
        | nucleona::range::endp;
        counter++;
    }

    return tails;
}


std::pair<std::string, bool> seat_adapter_auto_detect( 
                                      std::string& reads_path
                                    , size_t thread_num = 1
                                    )
{
    std::vector<std::string> tails;
    biovoltron::Index fm_index{12}, rc_fm_index{12};
    std::ifstream fm_ifs{index_prefix + ".table"}, rc_fm_ifs{index_prefix + ".rc_table"};
    fm_index.load(fm_ifs);
    rc_fm_index.load(rc_fm_ifs);
    biovoltron::Tailor tailor{fm_index, rc_fm_index};
    tailor.seed_len = seed_len;
    tailor.allow_seed_mismatch = !no_mismatch;

    if (is_gz_input)
    {
        boost::iostreams::filtering_istream ifs;

        ifs.push(boost::iostreams::gzip_decompressor());
        auto&& src(boost::iostreams::file_source(reads_path, std::ios_base::binary));
        if (!src.is_open())
            throw std::runtime_error("Can't open input gz file normally\n");

        ifs.push(src);
        if (!ifs.good())
            throw std::runtime_error("Can't open input gz stream normally\n");

        tails = tailor_pipeline(ifs, thread_num, tailor, DETECT_N_READS);
    }
    else
    {
        std::ifstream ifs(reads_path);
        if (!(ifs.is_open() && ifs.good()))
            throw std::runtime_error("Can't open input file normally\n");

        tails = tailor_pipeline(ifs, thread_num, tailor, DETECT_N_READS);
    }
    // std::cerr << "total number of tails sampled: " << tails.size() << "\n";

    std::string adapter;
    std::pair<std::string, bool> adapter_info;
    if (is_sensitive)
    {
        adapter_info = assemble_adapters<true>(tails, init_kmer_size, 5);
    }
    else
    {
        adapter_info = assemble_adapters<false>(tails, init_kmer_size, 3);
    }
    adapter = std::get<0>(adapter_info);

    if (adapter == "")
    {
        std::cout << "unable to detect adapter, use default adapter\n";
        adapter = DEFAULT_ADAPTER1;
    }
    else
    {
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