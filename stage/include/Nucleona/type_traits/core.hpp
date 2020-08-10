#pragma once
#include <vector>
#include <list>
#include <array>
#include <forward_list>
#include <deque>
#include <Nucleona/type_traits/has_member.hpp>
#include <Nucleona/type_traits/is_callable.hpp>
namespace nucleona { namespace type_traits {
template<class T>
struct IsSequenceContainerBase
{
    const static bool value = false;
};
template<class T, class... R>
struct IsSequenceContainerBase<
    std::vector<T, R...>
>
{
    const static bool value = true;
};
template<class T, class... R>
struct IsSequenceContainerBase<
    std::list<T, R...>
>
{
    const static bool value = true;
};
template<class T, std::size_t n>
struct IsSequenceContainerBase<
    std::array<T, n>
>
{
    const static bool value = true;
};
template<class T, class... R>
struct IsSequenceContainerBase<
    std::forward_list<T, R...>
>
{
    const static bool value = true;
};
template<class T, class... R>
struct IsSequenceContainerBase<
    std::deque<T, R...>
>
{
    const static bool value = true;
};
template<class T>
struct IsSequenceContainerBase<
    std::initializer_list<T>
>
{
    const static bool value = true;
};
template<class T>
struct IsSequenceContainer
{
    constexpr static bool value = IsSequenceContainerBase<std::decay_t<T>>::value;
};

template<class T>
constexpr static bool is_sequence_container_v = IsSequenceContainer<T>::value;


template <typename T, typename U>
struct decay_equiv : 
        std::is_same<std::decay_t<T>, std::decay_t<U>>::type 
{};
template <typename T, typename U>
constexpr auto decay_equiv_v = decay_equiv<T, U>::value;

template<class T, class ENABLE = void>
struct GetValueTypeBase
{
    typedef typename T::value_type value_type;
};
template<class T>
struct GetValueTypeBase<T, typename std::enable_if<!is_sequence_container_v<T>>::type>
{
    typedef void value_type;
};
template<class T>
using get_value_type_t = typename GetValueTypeBase<std::decay_t<T>>::value_type;
template<bool pred, class T>
using disable_if_t = std::enable_if_t<!pred, T>;
template<class TT, class T>
using TypeCheckP = std::enable_if_t<
      std::is_same<std::decay_t<TT>, T>::value
>;
template<class TT, class T>
using NTypeCheckP = std::enable_if_t<
      !std::is_same<std::decay_t<TT>, T>::value
>;
template<class T>
struct ReverseURefDeduce
{
    using Type = T;
};
template<class T>
struct ReverseURefDeduce<T&>
{
    using Type = T&;
};
template<class T>
struct ReverseURefDeduce<T&&>
{
    using Type = T;
};
template<class T>
using ReverseURefDeduceT = typename ReverseURefDeduce<T>::Type;
template<class T>
using RURDT = ReverseURefDeduceT<T>;
template<class T>
using DERET_T = RURDT<T>;
#define DERET_V(v) RURDT<decltype(v)>

template<class T>
struct GetReference
{
    using Type = typename T::reference;
};
namespace arma{ template<typename eT> class Col; }
template<class T>
struct GetReference<arma::Col<T>>
{
    using Type = T&;
};
}}
