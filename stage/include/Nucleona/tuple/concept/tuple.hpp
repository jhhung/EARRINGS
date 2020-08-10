#pragma once
#include <Nucleona/concept/core.hpp>
#include <tuple>
namespace nucleona{ namespace tuple{ namespace concept{

struct Tuple
{
    template<class... T>
    int helper( const std::tuple<T...>& t );
    template<class T>
    auto requires( T o ) NUCLEONA_EXPRS(
        helper(o)
    );
};


}}}
