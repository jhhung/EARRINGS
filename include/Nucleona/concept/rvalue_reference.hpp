#pragma once
#include <type_traits>
#include <Nucleona/concept/core.hpp>
namespace nucleona{ namespace  concept_ {

struct RvalueReference
{
    template<class T>
    auto requires_( T o ) NUCLEONA_EXPRS(
        nucleona::concept_::is_true<
            std::is_rvalue_reference< T >::value
        >
    );
};

}}
