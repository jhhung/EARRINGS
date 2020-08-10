#pragma once
#include <boost/asio.hpp>
#include <thread>
#include <Nucleona/language.hpp>
#include <Nucleona/tuple.hpp>
#include <future>
#include "asio_pool/future.hpp"
#include <iostream>
namespace nucleona::parallel {
    
template<class IOS>
struct BasicAsioPool {
    template<class Ret>
    using Future = asio_pool::Future<Ret>;
    BasicAsioPool(IOS&& ios, std::size_t thread_num)
    : ios_ (std::forward<IOS>(ios))
    , work_guard_(ios_)
    , workers_mux_()
    , worker_num_(thread_num)
    , workers_ ()
    {
        init_all_workers_(worker_num_);
    }
    BasicAsioPool(std::size_t thread_num)
    : ios_ ()
    , work_guard_(ios_)
    , workers_mux_()
    , worker_num_(thread_num)
    , workers_ ()
    {
        init_all_workers_(worker_num_);
    }

    template<class Func>
    auto submit(Func&& func) {
        using Ret = decltype(func());
        auto pm = std::make_shared<std::promise<Ret>>();
        asio_pool::Future<Ret> fut(pm->get_future(), ios_);
        ios_.post([ 
            pm_ = std::move(pm),
            func_ = FWD(func)
        ](){
            if constexpr(std::is_void_v<Ret>) {
                func_();
                pm_->set_value();
            } else {
                pm_->set_value(func_());
            }
        });
        return fut;
    }

    void flush() {
        while(ios_.poll() > 0) {}
        ios_.stop();
    }

    template<class... Args>
    auto poll_one(Args&&... args) {
        return ios_.poll_one(FWD(args)...);
    }
    template<class... Args>
    auto poll(Args&&... args) {
        return ios_.poll(FWD(args)...);
    }
    inline bool support_owner_as_worker() const {
        return true;
    }
    ~BasicAsioPool() {
        if(!ios_.stopped()) {
            ios_.stop();
        }
        for(auto&& t : workers_) {
            if(t->joinable()) {
                t->join();
            }
        }
    }
    auto get_parallization() const {
        return worker_num_;
    }
private:
    void init_all_workers_(std::size_t n ) {
        for(std::size_t i = 0; i < n; i ++) {
            workers_.emplace_back(new std::thread([this](){
                ios_.run();
            }));
        }

    }
    IOS ios_;
    typename std::decay_t<IOS>::work work_guard_;
    std::mutex workers_mux_;
    std::size_t worker_num_;
    std::vector<std::unique_ptr<std::thread>> workers_;
};

auto make_asio_pool(std::size_t thread_num, boost::asio::io_service& ios) {
    return BasicAsioPool<boost::asio::io_service&>(ios, thread_num);
}
auto make_asio_pool(std::size_t thread_num) {
    return BasicAsioPool<boost::asio::io_service>(thread_num);
}

}
