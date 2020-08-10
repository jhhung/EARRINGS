#pragma once
#include <Nucleona/parallel/thread_pool.hpp>
#include <Nucleona/language/dll_config.hpp>
namespace nucleona{ namespace parallel{ 
// int ThreadPool::max_thread_num(10);
// std::atomic<int> ThreadPool::available_thread_num(10);

// struct ThreadPoolImplInit
// {
//     ThreadPoolImplInit()
//     {
//         ThreadPool::max_thread_num = 10;
//         ThreadPool::available_thread_num = 10;
//     }
// };
// ThreadPoolImplInit initiallib;
API ThreadPool make_thread_pool( const std::size_t& n )
{
    auto& available_thread_num_ = available_thread_num();
    auto thread_space = available_thread_num_.load( std::memory_order_acquire );
    do
    {
        if ( thread_space < n ) 
            throw std::runtime_error("thread create fail, global thread space not enough.");
    }
    while( !available_thread_num_.compare_exchange_weak(
          thread_space, thread_space - n
        , std::memory_order_release
        , std::memory_order_relaxed
    ));
    return ThreadPool( n );
}
API ThreadPool make_empty_thread_pool()
{
    return ThreadPool();
}
API ThreadPool make_thread_pool()
{
    auto& available_thread_num_ = available_thread_num();
    auto thread_space = available_thread_num_.load( std::memory_order_acquire );
    while( !available_thread_num_.compare_exchange_weak(
          thread_space, 0
        , std::memory_order_release
        , std::memory_order_relaxed
    ));
    if ( thread_space > 0 )
        return ThreadPool( thread_space );
    else throw std::runtime_error("thread create fail, global thread space not enough.");
}
API bool set_max_thread_num( const int& n ) // not thread-safe
{
    auto& available_thread_num_ = available_thread_num();
    if ( n == ThreadPool::max_thread_num() ) ;
    else if ( n > ThreadPool::max_thread_num() )
    {
        auto old_max = ThreadPool::max_thread_num();
        ThreadPool::max_thread_num() = n;
        available_thread_num_ += (n - old_max);
    }
    else
    {
        auto old_ava = available_thread_num_.load( std::memory_order_acquire );
        do{
            auto use = ThreadPool::max_thread_num() - old_ava;
            auto new_ava = n - use;
            if ( new_ava < 0 ) 
                return false;
        }
        while( !available_thread_num_.compare_exchange_weak(
              old_ava
            , n - ( ThreadPool::max_thread_num() - old_ava )
            , std::memory_order_release 
            , std::memory_order_relaxed
        ));
        ThreadPool::max_thread_num() = n;
    }
    return true;
}
API std::atomic<int>& available_thread_num()
{
    static std::atomic<int> n(std::thread::hardware_concurrency());
    return n;
};

}}
