#ifndef GUNUNU_STATIC_IF_HPP
#define GUNUNU_STATIC_IF_HPP
// TODO namespace fix
#include <type_traits>

namespace gununu {
    template<bool b>
    struct bool_ { static const bool value = false; };
    template<> struct bool_<true> { static const bool value = true; };

    struct caller_false {
        template <bool E, typename F> caller_false call(bool_<E>, F&&){return caller_false();}
    };
    struct caller_true {
        template <typename F> caller_false call(bool_<true>, F&& f){ f(std::true_type{}); return caller_false(); };
        template <typename F> caller_true  call(bool_<false>, F&&) { return caller_true(); }
    };

    struct avoid_user_return {};
} //namespace gununu

#define STATIC_BREAKIF   do { return gununu::avoid_user_return(); (void)_gununu_arg_; } while(0); 

#define STATIC_IF(e)     do { gununu::caller_true().call(gununu::bool_<(e)>(),  [&](auto _gununu_arg_)mutable->gununu::avoid_user_return{
#define STATIC_ELSEIF(e) STATIC_BREAKIF }).call(gununu::bool_<(e)>(),  [&](auto _gununu_arg_)mutable->gununu::avoid_user_return{
#define STATIC_ELSE      STATIC_BREAKIF }).call(gununu::bool_<true>(), [&](auto _gununu_arg_)mutable->gununu::avoid_user_return{
#define STATIC_ENDIF     STATIC_BREAKIF }); } while(0);

#define LAZY_VALUE(v)    ((void)(_gununu_arg_), (v)) 
#define LAZY_TYPE(...)   typename std::conditional<_gununu_arg_ != std::false_type{}, __VA_ARGS__, void>::type

#endif //GUNUNU_STATIC_IF_HPP
