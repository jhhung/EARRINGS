#pragma once
#include <memory>
namespace nucleona{

template<class T>
std::unique_ptr<T> mkuniq( T* ptr )
{
    return std::unique_ptr<T>( ptr );
}
template<class T, class Deleter>
std::unique_ptr<T, Deleter> mkuniq( T* ptr, Deleter&& del )
{
    return std::unique_ptr<T, Deleter>( ptr,  std::forward<Deleter>(del));
}
template<class T>
std::unique_ptr<T[]> mkuniq_arr( T* ptr )
{
    return std::unique_ptr<T[]>( ptr );
}
template<class T, class Deleter>
std::unique_ptr<T[], Deleter> mkuniq_arr( T* ptr, Deleter&& del )
{
    return std::unique_ptr<T[], Deleter>( ptr,  std::forward<Deleter>(del));
}


}
