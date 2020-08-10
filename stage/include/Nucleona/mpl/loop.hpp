#pragma once
#include <Nucleona/tuple.hpp>
namespace nucleona::mpl {

constexpr struct TransformImpl {
    template<class Func, class... Args>
    inline auto operator()(Func&& func, Args&&... args) const {
        return nucleona::make_tuple(func(FWD(args))...);
    }
} transform_impl;
constexpr struct Transform {
    template<class Func, class... Args>
    inline auto operator()(Func&& func, std::tuple<Args...>& args) const {
        return std::apply([&func](auto&&... _args){
            return tranform_impl(func, FWD(_args)...);
        }, args);
    }
    template<class Func, class... Args>
    inline auto operator()(Func&& func, std::tuple<Args...>&& args) const {
        return std::apply([&func](auto&&... _args){
            return tranform_impl(func, FWD(_args)...);
        }, std::move(args));
    }
} transform;

constexpr struct Foreach {
    template<class Func, class... Args>
    inline void operator()(Func&& func, Args&&... args) const {
        return (void)nucleona::make_tuple(func(FWD(args))...);
    }
    template<class Func, class Int, Int... n>
    inline void operator()(
        Func&& func, 
        std::integer_sequence<Int, n...> idx
    ) const {
        return (void)nucleona::make_tuple(func(n)...);
    }
} foreach;

}
