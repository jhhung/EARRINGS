#pragma once
#include <Nucleona/range/core.hpp>
#include <Nucleona/range/concept.hpp>
#include <Nucleona/concept/core.hpp>
#include <Nucleona/mpl/integer_sequence.hpp>
#include <Nucleona/concept/rvalue_reference.hpp>
#include <Nucleona/tuple/apply.hpp>
#include <Nucleona/tuple/concept.hpp>

namespace nucleona{ namespace range{ namespace pipe_operator {
namespace nc_ = nucleona::concept;
namespace nr_ = nucleona::range;

template< class RNG, class ROP >
struct RangeTransform
{
    RNG rng;
    ROP rop;
    template<class... ARGS>
    decltype(auto) operator()( ARGS&&... args )
    {
        return FWD(rop)(FWD(rng), FWD(args)...);
    }
};
template< class RNG, class ROP >
auto make_range_transform( RNG&& rng, ROP&& rop )
{
    return RangeTransform< RNG, ROP >{ FWD(rng), FWD(rop) };
}


template<
      class RNG
    , class RNG_OP
    , NUCLEONA_CONCEPT_REQUIRE_( nr_::concept::Range, RNG )
    , NUCLEONA_CONCEPT_REQUIRE_( nucleona::tuple::concept::Tuple, RNG_OP )
>
decltype(auto) operator| ( RNG&& rng, RNG_OP&& op )
{
    return nucleona::tuple::detail::apply_impl(
          make_range_transform( FWD(rng), nucleona::tuple::forward<0>( op ) )
        , FWD(op)
        , nucleona::mpl::make_index_for<
              1
            , std::tuple_size<std::decay_t<RNG_OP>>::value
            , 1
        >() 
    );    
}


}}}

