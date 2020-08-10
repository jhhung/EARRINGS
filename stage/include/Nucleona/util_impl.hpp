#pragma once
#ifndef SINGLE_CPP
    #include <Nucleona/util.hpp>
#endif
#include <string>
namespace nucleona {
std::string indent(int il)
{
    static const std::string ele = "    ";
    std::string res;
    for ( int i = 0; i < il; i ++ )
        res += ele;
    return res;
}
}