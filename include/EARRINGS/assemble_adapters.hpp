#pragma once
#include <EARRINGS/graph.hpp>
#include <EARRINGS/common.hpp>
#include <string>
#include <vector>
#include <experimental/vector>
#include <array>
#include <string>
#include <fstream>
#include <unordered_map>

namespace EARRINGS
{
static constexpr std::array<char, 4> bases = {'A', 'T', 'C', 'G'};
static constexpr size_t MAX_KMER = 35;
bool detect_low_complexity(const std::string& adapter,
                           const size_t kmer_size, 
                           const float threshold = 0.7)
{
    size_t same_char(0);

    for (const auto base : bases)
    {
        same_char = 0;
        for (const char c : adapter.substr(0, kmer_size))
        {
            if (base == c)
            {
                same_char++;
            }
        }

        if ((float)same_char / kmer_size >= threshold)
        {
            return true;
        } 
    }

    return false;
}

size_t estimate_umi_size(const std::vector<std::string>& tails,
                         const std::string& adapter_seq)
{
    // len, count
    std::string sub_seq = adapter_seq.substr(0, 5);
    std::unordered_map<size_t, size_t> umi_size_cnt;
    for (const auto& t : tails)
    {
        if (t.size() < sub_seq.size())
        {
            continue;
        }

        for (size_t i = 0; i < t.size(); ++i)
        {
            auto sub_tail = t.substr(i, 5);
            if (sub_tail == sub_seq)
            {
                umi_size_cnt[i]++;
                break;
            }
        }   
    }
    
    size_t max_len = 0, max_value = 0;
    for (const auto& ele : umi_size_cnt)
    {
        // std::cerr << ele.first << "\t" << ele.second << "\n";
        if (ele.second > max_value)
        {
            max_value = ele.second;
            max_len = ele.first;
        }
    }

    return max_len;
}

// max_try is set to 5 in sensitive mode
template<bool IS_SENSITIVE>
std::pair<std::string, bool> assemble_adapters(
        std::vector<std::string>& tails, 
        const size_t kmer_size = 10,
        const size_t max_try = 3)
{
	std::vector<std::string> adapters;
    auto original_tails = tails;
    size_t tail_size = tails.size();
    kyutora::GraphWrapper g(kmer_size);
    
    for (size_t i(0); i < max_try; ++i)
    {
        // removing tails smaller than kmer size
        std::experimental::erase_if(tails, [&g](const auto& s){return s.size() < g.get_kmer_size();});
        
        // min_percentage of the adapter: prune_factor / max_try 
        float percentage = (prune_factor / (i + 1));

        // lowering PRUNE_FACTOR if adapters not found
        kyutora::PRUNE_FACTOR = size_t(std::ceil(tail_size * percentage));
        
        // prune factor too low, abort
        if constexpr (!IS_SENSITIVE)
        {
            if (kyutora::PRUNE_FACTOR < 10)
            {
                return std::make_pair("", false);
            }
        }
        
        g.build(tails);
        // g.print();
        g.find_sources_and_sinks();
        g.find_paths();
        adapters = g.get_adapters();


        if (!adapters.empty()) break; 
    }


    if (adapters.empty())
    {
        if (kmer_size >= MAX_KMER)
        {
            return std::make_pair("", false);
        }

        return assemble_adapters<IS_SENSITIVE>(original_tails, kmer_size + 5, max_try);
    }
    
    // increase k-mer when low complexity adapters are assembled
    if (detect_low_complexity(adapters[0], kmer_size))
    {
        if (kmer_size >= MAX_KMER)
        {
            return std::make_pair(adapters[0], true);
        }
        else
        {
            auto tmp_adapter = std::get<0>(assemble_adapters<IS_SENSITIVE>(original_tails, kmer_size + 5, max_try));
            if (tmp_adapter.empty())
            {
                return std::make_pair(adapters[0], true);
            }
            return std::make_pair(tmp_adapter, true);
        }
    } 
    
    return std::make_pair(adapters[0], false); 
}


}