#pragma once
#include "common.h"
#include <Nucleona/range/utils.hpp>
namespace nucleona::range {

template<class RNG>
auto make_any_range(RNG&& rng) {
    using Value = ValueT<std::decay_t<RNG>>;
    return ranges::any_view<Value>(FWD(rng));
}

}