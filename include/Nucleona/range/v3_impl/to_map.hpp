#pragma once
#include "common.h"
#include <Nucleona/range/utils.hpp>
#include <Nucleona/concept.hpp>
#include <map>
namespace nucleona::range {

struct ToMapFn {

    template<class Rng>
    using PairRangeConcept = ranges::view::keys_fn::Concept<Rng>;

    template<class Rng, class Less, CONCEPT_REQUIRES_(PairRangeConcept<Rng>())>
    auto operator()(Rng&& rng, Less&& less) const {
        using Pair = ValueT<std::remove_reference_t<Rng>>;
        using First = typename Pair::first_type;
        using Second = typename Pair::second_type;
        std::map<First, Second, Less> res(FWD(less));
        for(auto&& p : rng) {
            res.emplace(p.first, p.second);
        }
        return res;
    }
    template<class Less>
    auto operator()(Less&& less) const {
        return ranges::make_pipeable([
            _param = nucleona::make_tuple(*this, FWD(less))
        ](auto&& rng){
            auto& [_func, _less] =_param;
            return _func(FWD(rng), std::forward<Less>(_less));
        });
    }
    auto operator()() const {
        return ranges::make_pipeable([
            _param = nucleona::make_tuple(*this)
        ](auto&& rng){
            using Rng = decltype(rng);
            using Pair = ValueT<std::remove_reference_t<Rng>>;
            using Key = typename Pair::first_type;
            auto& [_func] =_param;
            return _func(FWD(rng), std::less<Key>());
        });
    }
};

RANGES_INLINE_VARIABLE(ToMapFn, to_map)
}