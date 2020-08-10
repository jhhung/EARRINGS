/**
 * @file    Nucleona/concept/iterator.hpp
 * @author  Chia-Hua Chang
 * @brief   Define the iterator concept
 */
#pragma once
#include <vector>
#include <Nucleona/concept/convertable.hpp>
#include <Nucleona/concept/core.hpp>
#include <Nucleona/concept/same.hpp>
namespace nucleona {
namespace concept {
/**
 * @brief Check input type is a iterator or not
 */
struct Iterator
{
    template<class T>
    auto requires( T o ) NUCLEONA_EXPRS(
          is_true<  std::is_copy_constructible<T>::value    >
        , is_true<  std::is_copy_assignable<T>::value       >
        , is_true<  std::is_destructible<T>::value          >
        , std::swap( o, o )
        , std::declval< typename std::iterator_traits<std::decay_t<T>>::value_type >()
        , std::declval< typename std::iterator_traits<std::decay_t<T>>::difference_type >()
        , std::declval< typename std::iterator_traits<std::decay_t<T>>::reference >()
        , std::declval< typename std::iterator_traits<std::decay_t<T>>::pointer >()
        , std::declval< typename std::iterator_traits<std::decay_t<T>>::iterator_category >()
        , *o
        , requires<Same>( ++o, std::declval<T&>() )
    );
};
/**
 * @brief Check input type is a input iterator or not
 */
struct InputIterator
{
    template < class T >
    auto requires ( T i ) NUCLEONA_EXPRS (
        nucleona::concept::requires< nucleona::concept::StaticCastable > ( i != i, bool{} ),
        nucleona::concept::requires< nucleona::concept::StaticCastable > ( i == i, bool{} ), *i,
        nucleona::concept::requires< Same > (
            i.operator-> (), std::declval< typename std::iterator_traits<
                                 std::decay_t< T > >::pointer > () ),
        nucleona::concept::requires< nucleona::concept::Same > ( ++i, std::declval< T& > () ), (void)i++, (void)++i,
        nucleona::concept::requires< nucleona::concept::Same > ( *i++, std::declval< typename std::iterator_traits<
                                     std::decay_t< T > >::value_type > () ) );
};

/**
 * @brief Check input type is random access or not
 */
struct RandomAccessIterator
{
    template < class T >
    auto requires ( T&& o ) NUCLEONA_EXPRS (
        *o, o++, o--, ++o, --o, o + int(), o - int(), o += int(), o -= int() );
};
}
}
