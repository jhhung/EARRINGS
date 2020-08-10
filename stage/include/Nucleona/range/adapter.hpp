#pragma once
#include <Nucleona/range/core.hpp>
#include <Nucleona/range/concept.hpp>
namespace nucleona{ namespace range{

template<class B_RANG, class ITER = typename std::decay_t<B_RANG>::iterator >
struct Adapter 
: public RangeProto< B_RANG >
{
    using iterator = ITER;
    using reference = typename std::iterator_traits< iterator >::reference;

    Adapter( B_RANG&& base_range )
    : RangeProto< B_RANG >{ FWD( base_range ) } 
    {}
    auto begin() const
    {
        return const_cast< Adapter<B_RANG, ITER>& >(*this).begin();
    }
    auto end() const 
    {
        return const_cast< Adapter<B_RANG, ITER>& >(*this).end();
    }
    auto begin()
    {
        return this->core.begin();
    }
    auto end()
    {
        return this->core.end();
    }
};

template < class T, NUCLEONA_CONCEPT_REQUIRE_( concept::Range, T ) >
auto make_adapter( T&& br )
{
    return Adapter< T >( FWD( br ) );
}



}}
