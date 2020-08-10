#pragma once
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

namespace nucleona{ namespace algo{

template<class STR, class T>
auto split( STR&& str, const T& deli ) {
    std::vector<std::string> res;
    boost::split(res, std::forward<STR>(str), boost::is_any_of(deli));
    return res;
}

}}