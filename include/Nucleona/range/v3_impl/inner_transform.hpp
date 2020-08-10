#pragma once
#include "common.h"
#include <Nucleona/tuple.hpp>
#include <vector>
namespace nucleona::range {

struct InnerTransformFn {
    template<class FUNC>
    auto operator()(FUNC&& func) const {
        return ranges::view::transform([
            _param = nucleona::make_tuple(FWD(func))
        ](auto&& inner_rng){
            auto& _func = std::get<0>(_param);
            using Elem = std::decay_t<decltype(_func(*inner_rng.begin()))>;
            std::vector<Elem> result;
            result.reserve(ranges::distance(inner_rng));
            for(auto&& v : inner_rng) {
                result.emplace_back(_func(v));
            }
            return result;
        });
    }
};
RANGES_INLINE_VARIABLE(InnerTransformFn, inner_transformed);

}