#pragma once
#include <utility>
namespace nucleona::language {
template<class T>
struct Make {
    template<class... Args>
    T operator()( Args&&... args) const {
        return T(std::forward<Args>(args)...);
    }
};

}