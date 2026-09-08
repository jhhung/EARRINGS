#pragma once
#include <EARRINGS/graph.hpp>
#include <EARRINGS/common.hpp>
#include <string>
#include <vector>
#include <experimental/vector>
#include <array>
#include <string>
#include <unordered_map>
#include <range/v3/all.hpp>

namespace EARRINGS {

static constexpr std::array<char, 4> bases = {'A', 'T', 'C', 'G'};
static constexpr size_t MAX_KMER = 35;
bool detect_low_complexity(
    const std::string& adapter,
    const size_t kmer_size, 
    const float threshold = 0.7
) {
    size_t same_char(0);

    for (const auto& base : bases) {
        same_char = 0;
        for (const char c : adapter.substr(0, kmer_size)) {
            if (base == c) ++same_char;
        }

        if ((float) same_char / kmer_size >= threshold) return true; 
    }

    return false;
}

// Returns {estimated tag3 length, reliable}. The length is the offset at which
// adapter3's prefix first appears across the tails (i.e. how many tag bases sit
// outside the assembled adapter). `reliable` is false when the signal is too
// weak or too flat to base the insert/tag boundary on.
std::pair<size_t, bool> estimate_tags3_len(
    const std::vector<std::string>& tails,
    const std::string& adapter_seq
) {
    constexpr size_t CHECK_LEN = 10;
    if (adapter_seq.size() < CHECK_LEN) return {0, false};

    const std::string_view adapter_prefix(adapter_seq.data(), CHECK_LEN);
    std::unordered_map<size_t, size_t> tags3_len_counts;
    for (const auto& tail : tails) {
        if (const auto pos = tail.find(adapter_prefix); pos != std::string::npos) {
            ++tags3_len_counts[pos];
        }
    }

    if (tags3_len_counts.empty()) return {0, false};

    size_t matched = 0;
    for (const auto& [pos, cnt] : tags3_len_counts) matched += cnt;

    const auto best = ranges::max_element(tags3_len_counts, {}, &std::pair<const size_t, size_t>::second);

    // The estimate is meaningful only when one offset clearly dominates: a sharp
    // peak means a fixed insert/tag boundary, a flat spread means we are guessing.
    const bool reliable = best->second * 2 >= matched;   // top offset holds >= half the hits

    return {best->first, reliable};
}

// max_try is set to 5 in sensitive mode
template<bool IS_SENSITIVE>
std::pair<std::string, bool> assemble_adapters(
    std::vector<std::string>& tails, 
    const size_t kmer_size = init_kmer_size,
    const size_t max_try = 3
) {
	std::vector<std::string> adapters;
    auto original_tails = tails;
    size_t tail_size = tails.size();
    kyutora::GraphWrapper g(kmer_size);
    
    for (size_t i(0); i < max_try; ++i) {
        // removing tails smaller than kmer size
        std::experimental::erase_if(tails, [&g](const auto& s){return s.size() < g.get_kmer_size();});
        
        // min_percentage of the adapter: prune_factor / max_try 
        float percentage = (prune_factor / (i + 1));

        // lowering PRUNE_FACTOR if adapters not found
        kyutora::PRUNE_FACTOR = size_t(std::ceil(tail_size * percentage));
        
        // prune factor too low, abort
        if constexpr (!IS_SENSITIVE) {
            if (kyutora::PRUNE_FACTOR < 10) {
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


    if (adapters.empty()) {
        if (kmer_size >= MAX_KMER) {
            return std::make_pair("", false);
        }
        return assemble_adapters<IS_SENSITIVE>(original_tails, kmer_size + kmer_step, max_try);
    }
    
    // increase k-mer when low complexity adapters are assembled
    if (detect_low_complexity(adapters[0], kmer_size)) {
        if (kmer_size >= MAX_KMER) {
            return std::make_pair(adapters[0], true);
        } else {
            auto tmp_adapter = std::get<0>(assemble_adapters<IS_SENSITIVE>(original_tails, kmer_size + kmer_step, max_try));
            if (tmp_adapter.empty()) {
                return std::make_pair(adapters[0], true);
            }
            return std::make_pair(tmp_adapter, true);
        }
    } 
    
    return std::make_pair(adapters[0], false); 
}

}