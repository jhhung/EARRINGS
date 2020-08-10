#pragma once
#include <Nucleona/range/irange.hpp>
#include <initializer_list>
#include <Nucleona/range/transform.hpp>
namespace nucleona{ namespace range{

constexpr struct At
{
    template<class RNG, class IDX_RNG>
    decltype(auto) operator()( 
          RNG&& rng
        , IDX_RNG&& ids
    ) const
    {
        auto idsnums = ids.size();
        return ( irange_0( idsnums ) 
        | transformed( 
            [ 
                  ids_ = FWD(ids) 
                , rng_ = FWD(rng)
            ]( const typename IDX_RNG::value_type& i ) mutable 
                -> decltype(auto)
            {
                return rng_[ ids_[i] ];
            }
        ) );
    }

} at_;
constexpr struct AtPipOp
{
    template< class INT >
    decltype(auto) operator()( 
          const std::initializer_list<INT>& ids
    ) const 
    {
        std::vector<INT> idxs( ids );
        return nucleona::make_tuple( at_, std::move( idxs ) );
    }
    template< class INT, class IDX_RNG >
    decltype(auto) operator()( IDX_RNG&& idx_rng ) const 
    {
        return nucleona::make_tuple( at_, FWD(idx_rng) );
    }
} at;


}}
