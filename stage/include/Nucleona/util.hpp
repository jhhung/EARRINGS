#pragma once
#include <vector>
#include <string>
#include <Nucleona/util/core.hpp>
#include <Nucleona/util/remove_const.hpp>
namespace nucleona {
using namespace nucleona::util;
template<class T>
void reset_vector( 
      std::vector<T>&   vec
    , const size_t&     n_cluster
    , const T&          val 
)
{
    vec.resize(n_cluster);
    std::fill( vec.begin(), vec.end(), val );
}
std::string indent(int il);
}
#define NOINLINE __attribute__ ((noinline))
#ifdef _MSC_VER
    #define EX(x) x
#endif

#ifdef SINGLE_CPP
    #include <Nucleona/util_impl.hpp>
#endif