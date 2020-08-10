#pragma once
#include "common.h"
namespace nucleona::range {
    template<class F>
    decltype(auto) transformed(F&& f) {
        return ranges::view::transform(FWD(f));
    }
}