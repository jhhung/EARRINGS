#pragma once
#include <Nucleona/tuple/apply.hpp>
#include <Nucleona/tuple/concept.hpp>
#include <Nucleona/tuple/forward.hpp>
#include <Nucleona/language.hpp>
namespace nucleona{

template<class... T>
constexpr auto make_tuple( T&&... o )
{
    return std::tuple<T...>{ FWD(o)... };
}


}
