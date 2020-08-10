/**
 * @file para_thread_pool.hpp
 * @brief The definition of thread pool implementation by parallel pool. 
 */
#pragma once
#include <thread>
#include <Nucleona/parallel/model/distsys/service.hpp>
#include <Nucleona/language.hpp>
#include <Nucleona/language/dll_config.hpp>
/**
 * @brief The ParaThreadPool is a thread pool which can parallize the job post by worker.
 * @details The ParaThreadPool is a implementation of ParallelPool, which worker is std::thread and use the basic trait.
 * The job function signature is void(void).
 * This implementation has the compatible interface of original version ThreadPool.
 * 
 */
namespace pokemon{ class ThreadPool; };
namespace nucleona{ namespace parallel{

using ThreadPoolProto = nucleona::parallel::model::distsys::Service< 
    nucleona::parallel::model::distsys::traits::BasicTrait< std::thread > 
>;
API std::atomic<int>& available_thread_num();
class ThreadPool : public ThreadPoolProto
{
API friend ThreadPool make_thread_pool( const std::size_t& n );
API friend ThreadPool make_thread_pool();
API friend bool set_max_thread_num( const int& n ); // not thread-safe
API friend ThreadPool make_empty_thread_pool();
friend pokemon::ThreadPool;
    using Base = ThreadPoolProto;
    static int& max_thread_num()
    {
        static int n(10);
        return n;
    };
protected:
    ThreadPool( const std::size_t n )
    : Base( n )
    {}
    ThreadPool()
    : Base()
    {}
public:
    template<class T>
    using Future = model::distsys::Future<T, nucleona::parallel::model::distsys::traits::BasicTrait< std::thread > >;
    DEFAULT_MOVE( ThreadPool );
    DISABLE_COPY( ThreadPool );
    
    template<class...ARGS>
    auto submit(ARGS... args)
    {
        return Base::job_post( FWD(args)... );
    }
    
    auto flush()
    {
        Base::flush_tasks();
    }
    bool resize_pool( std::size_t pool_size )
    {
        auto& available_thread_num_ = available_thread_num();
        auto tmp = Base::get_parallization();
        if ( pool_size == tmp ) ;
        else if ( pool_size < tmp )
        {
            Base::resize_pool( pool_size );
            available_thread_num_ += ( tmp - pool_size );
        }
        else
        {
            auto thread_space = available_thread_num_.load( std::memory_order_acquire );
            auto required_threads_num = pool_size - Base::get_parallization();
            do{
                if ( thread_space < required_threads_num ) 
                    return false;
            }
            while( !available_thread_num_.compare_exchange_weak(
                  thread_space, thread_space - required_threads_num
                , std::memory_order_release
                , std::memory_order_relaxed
            ));
            Base::resize_pool(pool_size);
        }
        return true;
    }
    ~ThreadPool()
    {
        try {
            flush();
            available_thread_num() += Base::get_parallization();
        } catch( const std::exception& e ) {

        }
    }
};


ThreadPool make_thread_pool( const std::size_t& n );
ThreadPool make_thread_pool();
ThreadPool make_empty_thread_pool();
bool set_max_thread_num( const int& n ); // not thread-safe

}}
