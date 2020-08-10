#pragma once
#include <boost/range/iterator_range.hpp>
#include <Nucleona/concept/reference.hpp>
#include <Nucleona/language.hpp>
namespace nucleona{ namespace range{

struct IteratorEnd{};


template<class CORE>
struct RangeProto
{
    CORE core;
};

template<class CORE, class ITR>
struct Range : 
      public RangeProto<CORE>
{
    using iterator = ITR;
    using const_iterator = ITR;
    using value_type = typename std::iterator_traits< iterator >::value_type;
    using reference = typename std::iterator_traits< iterator >::reference;
    using iterator_category = typename std::iterator_traits< iterator>::iterator_category;
    using difference_type = typename std::iterator_traits< iterator >::difference_type;

    Range( CORE&& core )
    : RangeProto<CORE>{ FWD( core ) }
    {}
    auto begin() const
    {
        return const_cast< Range<CORE, ITR>& >(*this).begin();
    }
    auto end() const 
    {
        return const_cast< Range<CORE, ITR>& >(*this).end();
    }
    auto begin()
    {
        return ITR( RangeProto<CORE>::core );
    }
    auto end()
    {
        return ITR( IteratorEnd{}, RangeProto<CORE>::core );
    }
    template<class I>
    decltype(auto) operator[]( I&& i )
    {
        return *(begin() + FWD(i));
    }
    template<class I>
    decltype(auto) operator[]( I&& i ) const
    {
        return const_cast< Range<CORE, ITR>& >(*this).operator[](FWD(i));
    }
};

template<class T, class... ARGS> // check T match iterator concept and has create
auto make( ARGS&&... args )
{
    auto&& core = T::make_range_core( FWD(args)... );
    return Range< 
          decltype( T::make_range_core( FWD(args)... ) ) 
        , T
    >( FWD(core) );
}




}
}
