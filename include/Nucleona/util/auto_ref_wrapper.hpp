#pragma once
#include <utility>
#include <Nucleona/concept.hpp>
#include <Nucleona/util.hpp>
namespace nucleona{ namespace util{

template<class T, class E = nucleona::concept::Enable> 
struct AutoRefWrapper 
{
    static_assert( 
          std::is_same<E, nucleona::concept::Enable>::value
        , "AutoRefWrapper static assert fail" 
    );
};

template<class T> 
struct AutoRefWrapper< T, NUCLEONA_CONCEPT_REQUIRE_NOT__( nucleona::concept::LvalueReference, T ) >
{
    using Storage = std::remove_const_t<T>;
    T data;
    bool enable;
    AutoRefWrapper( T&& d )
    : data( std::move( d ) )
    , enable( true )
    {}
    AutoRefWrapper( std::nullptr_t )
    : data()
    , enable( false )
    {}
    AutoRefWrapper()
    : data()
    , enable( false )
    {}
    auto set( std::nullptr_t )
    {
        nucleona::clear(data);
        enable = false;
    }

    auto set( T&& o )
    {
        return data = std::move(o);
    }
    const T& get() const
    {
        return data;
    }
    T& get()
    {
        return data;
    }
    auto empty() const
    {
        return !enable;
    }
};

template< class T > 
struct AutoRefWrapper< T, NUCLEONA_CONCEPT_REQUIRE__( nucleona::concept::LvalueReference, T ) > 
{
    using Data = std::remove_reference_t<T>;
    using Storage = Data*;
    Data* data;
    AutoRefWrapper( T& d )
    : data( &d )
    {}
    AutoRefWrapper( std::nullptr_t )
    : data( nullptr )
    {}
    AutoRefWrapper()
    : data( nullptr )
    {}
    auto set( std::nullptr_t )
    {
        data = nullptr;
    }
    auto set( Data& o )
    {
        return data = &o;
    }
    const auto& get() const 
    {
        return *data;
    }
    auto& get() {
        return *data;
    }
    auto empty() const
    {
        return data == nullptr;
    }
};
template<class T1, class T2>
bool equal( const AutoRefWrapper<T1>& t1, const AutoRefWrapper<T2>& t2 )
{
    if ( t1.empty() == t2.empty() )
    {
        if ( t1.empty() ) return true;
        else return t1.get() == t2.get() ;
    }
    else return false;
}

template<class T>
using AutoRef = AutoRefWrapper< T >;

}}
