#pragma once
#include <Nucleona/functional/apply.hpp>
#include <Nucleona/tuple/forward.hpp>
#include <Nucleona/mpl/integer_sequence.hpp>
namespace nucleona{ namespace tuple{
namespace detail{


template<class FUNC, class TUP, std::size_t... n>
decltype(auto) apply_impl( FUNC&& func, TUP&& tup, std::integer_sequence<std::size_t, n...> idx )
{
    return FWD(func)( FWD( nucleona::tuple::forward<n>(tup) )...  );
}


}
template<class FUNC, class TUP>
decltype(auto) apply( FUNC&& func, TUP&& tup )
{
    return std::apply(FWD(func), FWD(tup));
    // return apply_impl( 
    //       FWD(func)
    //     , FWD(tup)
    //     , nucleona::mpl::make_index_for<0, std::tuple_size< std::decay_t<TUP> >::value >()
    // );
}

}}
