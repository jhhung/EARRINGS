#pragma once
#include <future>
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{ namespace traits{

template<class TASK> struct GetFutureType{};

template<class R, class... P>
struct GetFutureType<std::packaged_task<R(P...)>>
{ 
    using Result = std::future<R>;

};

template<class T>
using get_future_type = typename GetFutureType<T>::Result;

template<class CORE>
struct CoreHelper
{
    static auto get_id()
    {}
};
template<>
struct CoreHelper<std::thread>
{
    static auto get_id()
    {
        return std::this_thread::get_id();
    }
};
}}}}}
