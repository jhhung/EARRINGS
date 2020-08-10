#pragma once
#include <type_traits>
namespace nucleona { namespace type_traits {
struct DummyType{};
template<class T, bool scalar_check = std::is_scalar<T>::value>
struct ToClass
{
    using Type = T;
};
template<class T>
struct ToClass<T, true>
{
    using Type = DummyType;
};
}}
