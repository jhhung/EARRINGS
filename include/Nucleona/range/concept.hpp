#pragma once 
#include <Nucleona/range/core.hpp>
#include <Nucleona/concept.hpp>
#include <Nucleona/concept/iterator.hpp>
namespace nucleona{ namespace range{ namespace  concept_ {
namespace nc_ = nucleona::concept_;
struct Range
{
    template<class T>
    auto requires_( T o ) NUCLEONA_EXPRS(
          nc_::requires_<nc_::Iterator>( o.begin() )
        , nc_::requires_<nc_::Iterator>( o.end()   )
    );
};
}
}}
