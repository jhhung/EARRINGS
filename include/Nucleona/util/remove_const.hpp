#pragma once
#include <utility>
namespace nucleona::util {

constexpr struct RemoveConst {
    template<class T>
    decltype(auto) operator()(T&& o) const {
        return FWD(o);
    }
    template<class T>
    decltype(auto) operator()(const T& o) const {
        return const_cast<T&>(o);
    }
} remove_const;

}