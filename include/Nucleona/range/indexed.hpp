#pragma once
#include <Nucleona/tuple.hpp>
#include <Nucleona/range/transform.hpp>
#include <Nucleona/language.hpp>

namespace nucleona{ namespace range{

template< class I, class V >
struct IndexPair
{
    I i;
    V v;
};

template< class I, class V >
auto make_index_pair( I&& i, V&& v )
{
    return IndexPair<I, V>{ FWD(i), FWD(v) };
}

constexpr struct Index
{
    template<class RNG, class START = int >
    decltype(auto) operator()( RNG&& rng, START&& start = 0 ) const
    {

        return nucleona::range::transform( 
              FWD( rng )
            , [ i = FWD(start) ]( auto&& obj ) mutable
              {
                  return make_index_pair( nucleona::copy(i++), FWD( obj ) );
              }
        );
    } 
} index;

constexpr struct IndexPipeOp
{
    template<class START = int >
    decltype(auto) operator()( START&& start = 0 ) const
    {
        return nucleona::make_tuple( index, FWD(start) );
    }
} indexed;

}}
