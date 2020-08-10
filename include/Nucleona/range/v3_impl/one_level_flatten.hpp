#pragma once
#include "common.h"
#include <Nucleona/tuple.hpp>
namespace nucleona::range {
    constexpr struct OneLevelFlattened {
        decltype(auto) operator()() const {
            return ranges::view::for_each([](auto&& v){
                return ranges::view::for_each(FWD(v), [](auto&& iv){
                    return ranges::yield(FWD(iv));
                });
            });
        }
    } one_level_flattened;
}