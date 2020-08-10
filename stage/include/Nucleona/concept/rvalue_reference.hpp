#pragma once
#include <type_traits>
#include <Nucleona/concept/core.hpp>
namespace nucleona{ namespace concept {

struct RvalueReference
{
    template<class T>
    auto requires( T o ) NUCLEONA_EXPRS(
        nucleona::concept::is_true< 
            std::is_rvalue_reference< T >::value
        >
    );
};

}}
