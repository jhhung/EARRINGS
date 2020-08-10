#pragma once
namespace nucleona{ namespace mpl{

template<class T, T... args>
struct TPack {
    using ElemType = T;
};

template<class... T>
struct TypeList{};

template<class T> struct Type{};

}}
