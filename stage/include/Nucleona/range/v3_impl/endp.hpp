#pragma once
#include "common.h"
namespace nucleona::range {

    struct EndP : public ranges::pipeable<EndP>
    {
        template<
            class Rng, 
            CONCEPT_REQUIRES_(ranges::Range<Rng>())
        >
        void operator()(Rng&& rng) const {
            for(auto&& v : rng) {}
        }
    };

    RANGES_INLINE_VARIABLE(EndP, endp)
}