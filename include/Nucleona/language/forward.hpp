#pragma once
#include <utility>
#include <Nucleona/type_traits/core.hpp>

#define FWD(v) \
std::forward<nucleona::type_traits::ReverseURefDeduceT<decltype(v)>>(v)

