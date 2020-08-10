#pragma once
#include <vector>
#include <cmath>
namespace nucleona::algo {

constexpr struct Segment {
    struct Result {
        std::size_t pos;
        std::size_t size;
    };
    std::vector<Result> operator()(
        std::size_t total, 
        std::size_t segment_num
    ) const {
        auto init_segment_size  = total / segment_num;
        auto remain             = total % segment_num;
        auto anchor_dist        = (float)segment_num / remain;

        std::vector<std::size_t> segment_size(segment_num, init_segment_size);
        for(decltype(remain) i = 0; i < remain; i ++ ) {
            segment_size[(std::size_t)std::round(i * anchor_dist)] += 1;
        }

        std::vector<Result> res;
        std::size_t pos = 0;
        for(auto& s : segment_size) {
            res.push_back(Result{pos, s});
            pos += s;
        }
        return res;
    }
} segment;
    
} 
