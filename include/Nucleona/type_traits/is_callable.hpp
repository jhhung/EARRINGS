#pragma once
#include <type_traits>
namespace nucleona { namespace type_traits {

template<typename T>
struct is_callable_impl {
  private:
    typedef char(&yes)[1];
    typedef char(&no)[2];

    struct Fallback { void operator()(); };
    struct Derived : T, Fallback { };

    template<typename U, U> struct Check;

    template<typename>
    static yes test(...);

    template<typename C>
    static no test(Check<void (Fallback::*)(), &C::operator()>*);

  public:
    static const bool value = sizeof(test<Derived>(0)) == sizeof(yes);
};
template<class T, class E = void>
struct is_callable_base 
{ 
};
template<class T>
struct is_callable_base<
    T
    , std::enable_if_t<
        !std::is_function<T>::value
        && !std::is_class<std::decay_t<T>>::value
    >
> : public std::false_type
{ };
template<class T>
struct is_callable_base<
    T
    , std::enable_if_t<
        !std::is_function<T>::value
        && std::is_class<std::decay_t<T>>::value
    >
> : public is_callable_impl<T>
{ };
template<class T>
struct is_callable_base<
    T
    , std::enable_if_t<
        std::is_function<T>::value
    >
> : public std::true_type
{ };

template<class T>
using is_callable = is_callable_base<T>;


}}
