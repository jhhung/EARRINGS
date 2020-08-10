/**
 * @file    Nucleona/concept/convertable.hpp
 * @author  Chia-Hua Chang
 * @brief   Define the type cast check concept.
 */
#pragma once
#include <Nucleona/concept/core.hpp>
namespace nucleona {
namespace concept {
/**
 * @brief Check the 2 input type is static castable or not
 */
struct StaticCastable
{
    template < class FROM, class TO >
    auto requires ( FROM&& f, TO&& t )
        NUCLEONA_EXPRS ( static_cast< TO > ( f ) );
};
}
}
