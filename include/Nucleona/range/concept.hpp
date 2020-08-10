#pragma once 
#include <Nucleona/range/core.hpp>
#include <Nucleona/concept.hpp>
#include <Nucleona/concept/iterator.hpp>
namespace nucleona{ namespace range{ namespace concept {
namespace nc_ = nucleona::concept;
struct Range
{
    template<class T>
    auto requires( T o ) NUCLEONA_EXPRS(
          nc_::requires<nc_::Iterator>( o.begin() )
        , nc_::requires<nc_::Iterator>( o.end()   )
    );
};
}
}}
