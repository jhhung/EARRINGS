#pragma once
#include <Nucleona/tuple.hpp>
namespace nucleona{ namespace range{

constexpr struct EndPipe
{
    template<class RNG>
    auto operator()( RNG&& rng ) const
    {
        for ( auto&& c : rng ){}
    }
} endpipe;

constexpr static auto endp = nucleona::make_tuple( endpipe );

}}
