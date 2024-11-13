#pragma once
#include <Nucleona/concept/core.hpp>
#include <type_traits>
namespace nucleona{ namespace  concept_{

struct Reference
{
    template<class T>
    auto requires_( T o ) NUCLEONA_EXPRS(
        nucleona::concept_::is_true<
            std::is_reference<T>::value
        >
    );
};


}}
