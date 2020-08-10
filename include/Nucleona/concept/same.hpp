/**
 * @file    Nucleona/concept/same.hpp
 * @author  Chia-Hua Chang
 * @brief   Define the concept to check the input type is same
 */
#pragma once
#include <type_traits>
#include <Nucleona/concept/core.hpp>
namespace nucleona {
namespace concept {

/**
 * @brief Used to check the input type is same.
 * @details Concept wrap of std::is_same
 */
struct Same
{
    template < class T1, class T2 >
    auto requires ( T1 t1, T2 t2 )
        NUCLEONA_EXPRS ( (is_true< std::is_same< T1, T2 >::value >));
};
}
}
