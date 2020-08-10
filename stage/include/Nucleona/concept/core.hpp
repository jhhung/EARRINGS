/**
 * @file    Nucleona/concept/core.hpp
 * @author  Chia-Hua Chang
 * @brief   The core part of concept
 */
#pragma once
#include <type_traits>
#include <Nucleona/concept/expr_validator.hpp>
#include <Nucleona/concept/has_requires.hpp>
#include <boost/mpl/eval_if.hpp>
#include <Nucleona/util.hpp>
/**
 * @def NUCLEONA_EXPRS()
 * This macro provides convenient way to create expression validate check.
 */
#ifdef _MSC_VER
    #define NUCLEONA_EXPRS(...) -> decltype( nucleona::concept::valid_expr( EX(NUCLEONA_VALID_EXPRS(__VA_ARGS__) )) )
#else
    #define NUCLEONA_EXPRS(...) -> decltype( nucleona::concept::valid_expr( NUCLEONA_VALID_EXPRS(__VA_ARGS__) ) )
#endif
/**
 * @def NUCLEONA_CONCEPT_REQUIRE__( CONCEPT_ID, ... )
 * This macro wrapped the EnableIf, which is used to judge a compile path should
 * remove or not by input concept and args
 */
#define NUCLEONA_CONCEPT_REQUIRE__( CONCEPT_ID, ... )            \
    nucleona::concept::EnableIf<                                 \
        nucleona::concept::satisfied< CONCEPT_ID, __VA_ARGS__ >, \
        nucleona::concept::Enable, CONCEPT_ID, __VA_ARGS__ >

/**
 * @def NUCLEONA_CONCEPT_REQUIRE_NOT__( CONCEPT_ID, ... )
 * Inverse check of NUCLEONA_CONCEPT_REQUIRE__
 */
#define NUCLEONA_CONCEPT_REQUIRE_NOT__( CONCEPT_ID, ... )         \
    nucleona::concept::EnableIf<                                  \
        !nucleona::concept::satisfied< CONCEPT_ID, __VA_ARGS__ >, \
        nucleona::concept::Enable, CONCEPT_ID, __VA_ARGS__ >

/**
 * @def NUCLEONA_CONCEPT_REQUIRE_( CONCEPT_ID, ... )
 * This macro has same purpose with NUCLEONA_CONCEPT_REQUIRE__ , but this macro
 * usually put in a template arguments declartion.
 */
#ifdef _MSC_VER
    #define NUCLEONA_CONCEPT_REQUIRE_( CONCEPT_ID, ... ) \
        EX(NUCLEONA_CONCEPT_REQUIRE__(CONCEPT_ID, __VA_ARGS__ )) = nucleona::concept::yes
#else
    #define NUCLEONA_CONCEPT_REQUIRE_( CONCEPT_ID, ... ) \
        NUCLEONA_CONCEPT_REQUIRE__(CONCEPT_ID, __VA_ARGS__ ) = nucleona::concept::yes
#endif

/**
 * @def NUCLEONA_CONCEPT_REQUIRE_( CONCEPT_ID, ... )
 * Inverse check of NUCLEONA_CONCEPT_REQUIRE_
 */
#ifdef _MSC_VER
    #define NUCLEONA_CONCEPT_REQUIRE_NOT_( CONCEPT_ID, ... ) \
        EX(NUCLEONA_CONCEPT_REQUIRE_NOT__(CONCEPT_ID,__VA_ARGS__)) = nucleona::concept::yes
#else
    #define NUCLEONA_CONCEPT_REQUIRE_NOT_( CONCEPT_ID, ... ) \
        NUCLEONA_CONCEPT_REQUIRE_NOT__(CONCEPT_ID,__VA_ARGS__) = nucleona::concept::yes
#endif

/**
 * @def NUCLEONA_CONCEPT_REQUIRE( CONCEPT_ID, ... )
 * This macro has same purpose with NUCLEONA_CONCEPT_REQUIRE_ , but is a
 * template declartion.
 * This macro usually used to check a member method, static member variable or
 * nested type declartion.
 */
#ifdef _MSC_VER
    #define NUCLEONA_CONCEPT_REQUIRE( CONCEPT_ID, ... ) \
        template<EX(NUCLEONA_CONCEPT_REQUIRE_(CONCEPT_ID, __VA_ARGS__ ))>
#else
    #define NUCLEONA_CONCEPT_REQUIRE( CONCEPT_ID, ... ) \
        template<NUCLEONA_CONCEPT_REQUIRE_(CONCEPT_ID, __VA_ARGS__ )>
#endif

/**
 * @def NUCLEONA_CONCEPT_REQUIRE_NOT( CONCEPT_ID, ... )
 * Inverse check of NUCLEONA_CONCEPT_REQUIRE
 */
#ifdef _MSC_VER
    #define NUCLEONA_CONCEPT_REQUIRE_NOT( CONCEPT_ID, ... ) \
        template<EX(NUCLEONA_CONCEPT_REQUIRE_NOT_(CONCEPT_ID,__VA_ARGS__))>
#else
    #define NUCLEONA_CONCEPT_REQUIRE_NOT( CONCEPT_ID, ... ) \
        template<NUCLEONA_CONCEPT_REQUIRE_NOT_(CONCEPT_ID,__VA_ARGS__)>
#endif
/** 
 * @def NUCLEONA_CONCEPT_ASSERT( CONCEPT_ID, ... )
 * This macro create static assert of a concept check.
 */
#define NUCLEONA_CONCEPT_ASSERT( CONCEPT_ID, ... )                           \
    static_assert ( nucleona::concept::satisfied< CONCEPT_ID, __VA_ARGS__ >, \
        __FILE__ ":" __LINE__ ":" #CONCEPT_ID " not satisfied" )
/**
 * @def NUCLEONA_CONCEPT_ASSERT_NOT( CONCEPT_ID, ... )
 * Inverse check of NUCLEONA_CONCEPT_ASSERT
 */
#define NUCLEONA_CONCEPT_ASSERT_NOT( CONCEPT_ID, ... )                        \
    static_assert ( !nucleona::concept::satisfied< CONCEPT_ID, __VA_ARGS__ >, \
        __FILE__ ":" __LINE__ ":" #CONCEPT_ID " not satisfied" )

namespace nucleona {
namespace concept {
/**
 * This is a dummy enum type for better look of concept check expression.
 */
enum Enable
{
    yes
};
using ValidExprRet = void;

/**
 * @brief Expression check of concept
 * @details All expression in arguments should compiled and could return a non
 * void result or this function will compile failed.
 */
template < class... T >
constexpr ValidExprRet valid_expr ( T... o );
/**
 * @brief Convert expression check result( compile pass or failed ) into boolean
 * variable
 * @details For user friendly, the concept check result should be a boolean
 * state, and can reuse the state easier.
 * It calls the type HasRequires to check the compile path exist or not.
 */
template < class CONCEPT, class... T >
constexpr bool satisfied =
    nucleona::concept::HasRequires< CONCEPT, ValidExprRet ( T... ) >::value;

/**
 * @brief Convert boolean state back to the compile pass or failed state.
 * @details An interface transform for std::enable_if.
 */
template < bool b, std::enable_if_t< b, int > = 0 >
constexpr bool  is_true = b;

/**
 * @brief Alias function call of the requires in concept type.
 * @details Same purpose of CONCEPT::requires, object based concept check.
 */
template < class CONCEPT, class... T >
constexpr auto requires ( T... o )
    NUCLEONA_EXPRS ( std::declval< CONCEPT > ().requires ( o... ) );
/**
 * @brief Customized enable if
 * @details Same purpose of std::enable_if, but provide batter error output for
 * concept check
 */
template < bool b, class T, class CONCEPT, class... CONCEPT_ARGS >
struct EnableIfProto
{
};
template < class T, class CONCEPT, class... CONCEPT_ARGS >
struct EnableIfProto< true, T, CONCEPT, CONCEPT_ARGS... >
{
    using Type = T;
};

template < class T, class CONCEPT, class... CONCEPT_ARGS >
struct EnableIfProto< false, T, CONCEPT, CONCEPT_ARGS... >
{
};

/**
 * @brief Customized enable if
 * @details wrapped of EnableIfProto
 */
template < bool b, class T, class CONCEPT, class... CONCEPT_ARGS >
using EnableIf = typename EnableIfProto< b, T, CONCEPT, CONCEPT_ARGS... >::Type;
}
/**
 * @brief Template interface for concept check
 * @details Convert the concept expression compile result into boolean state,
 * and pass to enable if to make the function or class compiled or not.
 */
template < class CONCEPT, class... T >
using concept_require =
    concept::EnableIf< nucleona::concept::satisfied< CONCEPT, T... >,
        concept::Enable, CONCEPT, T... >;
/**
 * @brief Template interface for concept check
 * @details Inverse check of concept_require
 */
template < class CONCEPT, class... T >
using concept_require_not =
    concept::EnableIf< !nucleona::concept::satisfied< CONCEPT, T... >,
        concept::Enable, CONCEPT, T... >;
}
