#pragma once
#include <Nucleona/exception.hpp>

#if defined( __GNUC__) || defined(__clang__)
    #define NUCLEONA_DEBUG_THROW(e) \
        throw nucleona::make_debug_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__, e)
#elif defined(_MSC_VER)
    #define NUCLEONA_DEBUG_THROW(e) \
        throw nucleona::make_debug_exception(__FILE__, __LINE__, __FUNCION__)
#endif

#if NUCLEONA_USE_DEBUG_EXCEPTION
#   define NUCLEONA_THROW(e) NUCLEONA_DEBUG_THROW(e)
#else 
#   define NUCLEONA_THROW(e) throw e
#endif