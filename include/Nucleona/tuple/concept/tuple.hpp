#pragma once
#include <Nucleona/concept/core.hpp>
#include <tuple>
namespace nucleona{ namespace tuple{ namespace  concept_{

struct Tuple
{
    template<class... T>
    int helper( const std::tuple<T...>& t );
    template<class T>
    auto requires_( T o ) NUCLEONA_EXPRS(
        helper(o)
    );
};


}}}
