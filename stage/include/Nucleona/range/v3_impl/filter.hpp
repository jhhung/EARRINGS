#pragma once
#include "common.h"
#include <Nucleona/tuple.hpp>
namespace nucleona::range {
    template<class Pred>
    decltype(auto) filtered(Pred&& pred) {
        return ranges::view::for_each([
            _param = nucleona::make_tuple(FWD(pred))
        ](auto&& v){
            auto& [ _pred ] = _param;
            bool tmp = _pred(v);
            return ranges::yield_if(tmp, FWD(v));
        });
    }
}