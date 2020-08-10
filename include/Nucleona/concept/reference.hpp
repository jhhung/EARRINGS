#pragma once
#include <Nucleona/concept/core.hpp>
#include <type_traits>
namespace nucleona{ namespace concept{

struct Reference
{
    template<class T>
    auto requires( T o ) NUCLEONA_EXPRS(
        nucleona::concept::is_true<
            std::is_reference<T>::value
        >
    );
};


}}
