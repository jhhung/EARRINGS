#pragma once
#include <thread>
#include <future>
#include <boost/lockfree/spsc_queue.hpp>
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{ namespace traits{
template<class CORE>
struct BasicTrait
{};
template<>
struct BasicTrait<std::thread>
{
    typedef std::thread                     Core;
    template<class T>
    using Queue = boost::lockfree::spsc_queue<T>;
};
// template<>
// struct BasicTrait<boost::thread>
// {
// 	typedef std::packaged_task<void(void)> _TaskStore;
// 	typedef std::size_t _TaskId;
// 	typedef boost::thread _ParallelType;
// 	
// 	static void worker_yield()
// 	{
// 		boost::this_thread::yield();
// 	}
// };

}}}}}
