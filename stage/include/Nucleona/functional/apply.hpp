#pragma once
#include <Nucleona/language.hpp>
namespace nucleona{ namespace functional{

template<class FUNC, class... ARGS>
auto apply( FUNC&& func, ARGS&&... args )
{
    return FWD(func)( FWD(args)... );
}

}}
