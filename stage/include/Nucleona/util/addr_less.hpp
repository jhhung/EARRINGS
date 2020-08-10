#pragma once

namespace nucleona { namespace util {
struct AddrLess
{
    template<class T, class Y>
    bool operator()( T&& t, Y&& y ) const
    {
        return &t < &y;
    }
};
} /* util */ 
} /* nucleona */ 
