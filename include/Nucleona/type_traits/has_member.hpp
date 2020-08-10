#pragma once
#include <type_traits>
#include <Nucleona/type_traits/type.hpp>
namespace nucleona { namespace type_traits {

template <typename... Args> 
struct ambiguate : public ToClass<Args>::Type... 
{};
template<typename A, typename = void>
struct got_type : std::false_type {};

template<typename A>
struct got_type<A> : std::true_type {
        typedef A type;
};
template<typename T, T>
struct sig_check : std::true_type {};
template<typename Alias, typename AmbiguitySeed>
struct has_member {
    template<typename C> static char ((&f(decltype(&C::value))))[1];
    template<typename C> static char ((&f(...)))[2];
    static_assert(
            (sizeof(f<AmbiguitySeed>(0)) == 1)
            , "Member name specified in AmbiguitySeed is different from member name specified in Alias, or wrong Alias/AmbiguitySeed has been specified."
            );

    static bool const value = sizeof(f<Alias>(0)) == 2;
};

}}

#define CREATE_MEMBER_CHECK(member)                                         \
template<typename T, typename = std::true_type>                             \
struct Alias_##member;                                                      \
template<typename T>                                                        \
struct Alias_##member <                                                     \
    T, std::integral_constant<bool, nucleona::type_traits::got_type<decltype(&T::member)>::value>  \
> { static const decltype(&T::member) value; };                             \
struct AmbiguitySeed_##member { char member; };                             \
template<typename T>                                                        \
struct has_member_##member {                                                \
    static const bool value                                                 \
        = nucleona::type_traits::has_member<                                                       \
            Alias_##member<nucleona::type_traits::ambiguate<T, AmbiguitySeed_##member>>            \
            , Alias_##member<AmbiguitySeed_##member>                        \
        >::value                                                            \
    ;                                                                       \
};                                                                          \
template<class T>                                                           \
constexpr bool has_member_##member##_v = has_member_##member<T>::value;       

