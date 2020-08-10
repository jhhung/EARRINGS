#pragma once
#include <Nucleona/tuple.hpp>
#include <Nucleona/range/transform.hpp>
#include <Nucleona/language.hpp>

namespace nucleona{ namespace range{

constexpr struct Print
{
    template<class RNG, class OS, class DEL = const char* >
    decltype(auto) operator()( RNG&& rng, OS&& os, DEL&& del = "" ) const
    {
        struct Pack // workaround for lambda
        { 
            OS os; 
            DEL del;
        } pack { FWD(os), FWD(del) };

        return nucleona::range::transform( 
              FWD( rng )
            , [
                  p = std::move(pack)
              ]( auto&& obj ) mutable
              {
                  p.os << obj << p.del;
                  return obj;
              }
        );
    } 
} print;

constexpr struct PrintPipeOp
{
    template<class OS, class DEL = const char* >
    decltype(auto) operator()( OS&& os, DEL&& del = "" ) const
    {
        return nucleona::make_tuple( print, FWD(os), FWD(del) );
    }
} printed;

}}
