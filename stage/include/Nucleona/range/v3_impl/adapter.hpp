#pragma once
#include "common.h"
namespace nucleona::range {
    template<class T>
    decltype(auto) make_adapter( T&& obj ) {
        return FWD(obj) | ranges::view::all ;
    }
}