#pragma once
#include <thread>
#include <future>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{ namespace traits{
template<class CORE>
struct DependencyTrait
{};
template<>
struct DependencyTrait<std::thread>
{
	typedef std::packaged_task<void(void)>  TaskStore;
	typedef std::size_t                     TaskId;
    typedef std::thread                     Core;

    template<class T>
    using Queue = boost::lockfree::spsc_queue<T>;
};
}}}}}
