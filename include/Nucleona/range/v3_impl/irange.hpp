#pragma once
#include "common.h"
namespace nucleona::range {

template<class I>
decltype(auto) irange_0( I&& i) {
    return ranges::view::ints((std::decay_t<I>)0, std::forward<I>(i) );
}
// template<class... T>
// decltype(auto) irange( T&&... o)
// {
//     return boost::irange( std::forward<T>(o)... );
// }

}