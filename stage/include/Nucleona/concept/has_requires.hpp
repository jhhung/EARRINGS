/**
 * @file    Nucleona/concept/has_requires.hpp
 * @author  Chia-Hua Chang
 * @brief   Check input type has method "requires", or not.
 */
#pragma once
#include <type_traits>
namespace nucleona {
namespace concept {
/**
 * @brief The default type of HasRequires, is a error assert section which
 * should never be compiled.
 */
template < typename, typename T >
struct HasRequires
{
    static_assert ( std::integral_constant< T, false >::value,
        "Second template parameter should be a function type." );
};
/**
 * @brief The main check section of the HasRequires
 * @details This struct check the class C has the method "requires" and match
 * the input function type "Ret(Args...)", or not
 * This is one of the main part of the concept check.
 * @tparam C Concept need to check
 * @tparam Ret The return type which the requires function should match
 * @tparam Args The arguments types which the requires function should have
 */
template < typename C, typename Ret, typename... Args >
struct HasRequires< C, Ret ( Args... ) >
{
  private:
    template < typename T >
    static constexpr auto check ( T* ) ->
        typename std::is_same< 
            decltype ( 
                std::declval< T > ().template requires< Args... > (
                    std::declval< Args > ()... 
                ) 
            )
            , Ret 
        >::type;

    template < typename >
    static constexpr std::false_type check ( ... );

    typedef decltype ( check< C > ( 0 ) ) type;

  public:
    static constexpr bool value = type::value;
};
}
}
