#pragma once
#include "common.h"
#include <Nucleona/tuple.hpp>
namespace nucleona::range {

template<class OutS, class Delim = const char*>
decltype(auto) printed(OutS&& os, Delim&& del = "") {
    return ranges::view::transform([
        _param = nucleona::make_tuple(
            FWD(os), FWD(del)
        )
    ](auto&& v){
        auto& [_os, _del] = _param;
        _os << v << _del;
        return FWD(v);
    });
}

}