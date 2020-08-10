#pragma once
#include <tuple>
namespace nucleona{ namespace tuple{

template<std::size_t i, class TUP>
decltype(auto) forward( TUP&& tup )
{
    return std::forward<
        std::tuple_element_t<
              i
            , std::decay_t<TUP>
        >>( std::get<i>(tup) );
}

}}
