#pragma once
#include <Nucleona/range.hpp>
namespace nucleona::algo::string_sort {
struct MultikeyQuickSort {
    MultikeyQuickSort(std::size_t max_depth)
    : max_depth_(max_depth)
    {}
    template<class Executor, class StrRng>
    void call_impl(
        Executor&           extor,
        StrRng&&            __str_rng, 
        const std::size_t&  depth,
        const std::size_t&  lower,
        const std::size_t&  upper
    ) const {
        auto&& str_rng = ranges::view::slice(__str_rng, lower, upper);
        auto str_num = upper - lower;
        if(str_num <= 1 || depth >= max_depth_) return;
        auto pivot = pick_pivot(str_rng, depth);
        auto [i, j] = partition(str_rng, depth, pivot, str_num);
        auto fut_l = extor.submit(
            [&](){
                call_impl(extor, __str_rng, depth,     lower,     lower + i);
            }
        );
        auto fut_m = extor.submit(
            [&](){
                call_impl(extor, __str_rng, depth + 1, lower + i, lower + j);
            }
        );
        call_impl(extor, __str_rng, depth,     lower + j, lower + str_num);
        fut_l.sync();
        fut_m.sync();
    }
    template<class Executor, class StrRng>
    void operator()(Executor& extor, StrRng&& str_rng, std::size_t depth) const {
        return call_impl(extor, str_rng, depth, 0, ranges::distance(str_rng));
    }
    template<class StrRng>
    auto pick_pivot(StrRng&& str_rng, std::size_t depth) const {
        auto& str = *str_rng.begin();
        return str[depth];
    }
    template<class StrRng, class Pivot>
    auto partition(StrRng&& str_rng, std::size_t depth, Pivot pivot, std::size_t str_num) const {
        std::size_t i = 0;
        std::size_t j = 0;
        std::size_t n = str_num - 1;
        while(j <= n) {
            auto&& str_i = str_rng[i];
            auto&& str_j = str_rng[j];
            auto&& v_i = str_i[depth];
            auto&& v_j = str_j[depth];
            if( v_j < pivot) {
                std::swap(str_i, str_j);
                i = i + 1;
                j = j + 1;
            } else if( v_j > pivot) {
                std::swap(str_j, str_rng[n]);
                n = n - 1;
            } else {
                j = j + 1;
            }
        }
        return std::make_tuple(i, j);

    }
private:
    std::size_t max_depth_;
};

}
