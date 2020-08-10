#pragma once
#include <string>
#include <Tailor/tailor.hpp>
#include <Tailor/range.hpp>
#include <experimental/vector>
#include <EARRINGS/graph.hpp>
#include <EARRINGS/common.hpp>
#include <EARRINGS/assemble_adapters.hpp>
#include <Nucleona/range/v3_impl.hpp>
#include <Nucleona/parallel/thread_pool.hpp>
#include <Nucleona/parallel/asio_pool.hpp>
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

template<class FORMAT, class TailorMain>
std::vector<std::string> tailor_pipeline(std::string& reads_path
                                       , size_t thread_num
                                       , TailorMain&& tailor_mapping
                                       , size_t num_reads) 
{
    const auto& aligner = tailor_mapping.get_table();
    const auto& paras = tailor_mapping.get_paras();
    auto tp = nucleona::parallel::make_asio_pool(thread_num);
    using DataType = typename std::remove_const_t<std::remove_reference_t<decltype(aligner)>>::FASTQ;
    constexpr pipeline::range::format_reader_fn<DataType> format_reader{};  

    std::ifstream input(reads_path);
    std::vector<std::string> tails;
    tails.reserve(num_reads);
    size_t counter(0);
    while(input.good() && tails.size() < 3000 && counter < 3)
    {
        input
        | format_reader()
        | pipeline::range::align(aligner)
        | ranges::view::transform(
            [&tails](auto&& AlignedReads_v) {
            for(auto&& i : AlignedReads_v)
            {
                if (i.tail_pos_ >= 0)
                {
                    tails.emplace_back(i.fq_.seq.substr(i.fq_.seq.length() - i.tail_pos_ - 1));
                }
            }
            return AlignedReads_v;
            }
        )	
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
    if (is_fastq)
    {
        tailor::TailorMain<falseType::value> tailor_mapping(thread_num, seed_len, min_multi, index_prefix, enable_mismatch);
        tails = tailor_pipeline<tailor::Fastq>(reads_path, thread_num, tailor_mapping, DETECT_N_READS);
    }
    else
    {
        tailor::TailorMain<trueType::value> tailor_mapping(thread_num, seed_len, min_multi, index_prefix, enable_mismatch);
        tails = tailor_pipeline<tailor::Fasta>(reads_path, thread_num, tailor_mapping, DETECT_N_READS);
    }

    // std::cerr << "total number of tails sampled: " << tails.size() << "\n";

    std::string adapter;
    std::pair<std::string, bool> adapter_info;
    if (is_sensitive)
    {
        adapter_info = assemble_adapters<true>(tails, 10, 5);
    }
    else
    {
        adapter_info = assemble_adapters<false>(tails, 10, 3);
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