#pragma once
#include <Nucleona/memory/cdelete.hpp>
namespace nucleona { namespace util { 

class Cdeleter
{
public:
    template<class T>
    void operator()( T* ptr )
    {
        nucleona::memory::cdelete( ptr );
    }
};


} /* util */ 
} /* nucleona */ 
