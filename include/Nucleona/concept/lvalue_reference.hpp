#pragma once
#include <type_traits>
#include <Nucleona/range/core.hpp>
namespace nucleona{ namespace  concept_ {

struct LvalueReference
{
    template<class T>
    auto requires_( T o ) NUCLEONA_EXPRS(
        nucleona::concept_::is_true<
            std::is_lvalue_reference< T >::value
        >
    );
};

}}
