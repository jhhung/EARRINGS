#pragma once
#include "common.h"
#include <Nucleona/tuple.hpp>
namespace nucleona::range {


struct IndexFn
{
    template<
        class Rng, 
        class Idx = int,
        CONCEPT_REQUIRES_(ranges::InputRange<Rng>())
    >
    auto operator()(
        Rng&& rng,
        Idx start = 0
    ) const {
        using Idx_ = decltype(ranges::distance(rng));
        return ranges::view::zip(
            ranges::view::ints((Idx_)start, (Idx_)(start + ranges::distance(rng))), 
            FWD(rng)
        );
    }
    template<class Idx = int>
    auto operator()(Idx start = 0) const {
        return ranges::make_pipeable([
            ___this = this, start
        ](auto&& rng){
            return (*___this)(FWD(rng), start);
        });
    }
};
RANGES_INLINE_VARIABLE(IndexFn, indexed)

}