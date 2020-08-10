#pragma once
#include <vector>
#include <string>
namespace nucleona {
namespace util{
template<class T>
void clear( T& v ) {}

template<class T>
void clear( std::vector<T>& v )
{
    std::vector<T> blank;
    v.swap(blank);
}
}
}
