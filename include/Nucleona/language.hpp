#pragma once
#include <Nucleona/language/dir.hpp>
#include <Nucleona/language/forward.hpp>
#include <Nucleona/language/static_if.hpp>
#include <Nucleona/language/copy.hpp>
#include <Nucleona/util.hpp>
#include <Nucleona/language/getter.hpp>
#include <Nucleona/language/make.hpp>
#define DISABLE_COPY(T) \
    T(T const &) = delete; \
    T& operator=(T const &t) = delete;

#define DISABLE_MOVE(T) \
    T(T&&) = delete; \
    T& operator=(T&& t) = delete;

#define DEFAULT_COPY(T) \
    T(T const &) = default; \
    T& operator=(T const &t) = default;

#define DEFAULT_MOVE(T) \
    T(T&&) = default; \
    T& operator=(T&& t) = default;

#define CREATE_DERIVED_TYPE_BODY(TYPE, PARENT_TYPE) \
    using Base = PARENT_TYPE; \
    using Base::Base; \
    TYPE() = default; \
    DEFAULT_COPY(TYPE) \
    DEFAULT_MOVE(TYPE)

#define CREATE_DERIVED_TYPE(TYPE, PARENT_TYPE) \
    struct TYPE : public PARENT_TYPE \
    { \
        CREATE_DERIVED_TYPE_BODY(TYPE, PARENT_TYPE) \
    };

#ifdef _MSC_VER
    #define APPLY(callable, ...) \
        EX(callable(__VA_ARGS__))
#else
    #define APPLY(callable, ...) \
        callable(__VA_ARGS__)
#endif

#define FOREACH_BASIC_TYPE(callable) \
callable(std::string) \
callable(float) \
callable(double) \
callable(int8_t) \
callable(int16_t) \
callable(int32_t) \
callable(int64_t) \
callable(uint8_t) \
callable(uint16_t) \
callable(uint32_t) \
callable(uint64_t) \

#define __SDIR__ sdir(__FILE__)
#define COMMA ,

namespace nucleona
{
    using language::copy;
} /*  nucleona */ 

#define FUNC_ALIAS(new_f, old_f) \
template<class...T> \
decltype(auto) new_f(T&&... o) \
{ return old_f(std::forward<T>(o)...); }


#define VDUMP( exp ) #exp << ": " <<  exp

#define STR(sym) #sym
#define UNWRAP_SYM_STR(sym) STR(sym)