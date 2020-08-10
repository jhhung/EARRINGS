#pragma once
#include <utility>
#include <Nucleona/language.hpp>
namespace nucleona{
namespace functional {
template<class T>
struct MutableStorage
{
    mutable T storage;
};
template<class T>
struct MutableStorage<T&>
{
    T& storage;
};
template<class T>
auto make_mutable_storage(T&& o)
{
    return MutableStorage<T>{std::forward<T>(o)};
}
FUNC_ALIAS(mms, make_mutable_storage);
// template<class T>
// constexpr auto mms = make_mutable_storage<T>;

template<class F>
struct Mutable : public MutableStorage<F>
{
    using Base = MutableStorage<F>;
    Mutable(F&& f)
    : Base{ std::forward<F>(f) }
    {}
    template<class... ARGS>
    decltype(auto) operator()(ARGS&&... args) const
    {
        return this->storage(std::forward<ARGS>(args)...);
    }
};

template<class F>
auto mutable_(F&& f)
{
    return Mutable<F>(std::forward<F>(f));
}
}}
