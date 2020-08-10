#pragma once
#include <type_traits>
#include <Nucleona/range/core.hpp>
namespace nucleona{ namespace concept {

struct LvalueReference
{
    template<class T>
    auto requires( T o ) NUCLEONA_EXPRS(
        nucleona::concept::is_true< 
            std::is_lvalue_reference< T >::value
        >
    );
};

}}
