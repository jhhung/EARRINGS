#pragma once
#include <future>
#include <Nucleona/language.hpp>
#include <boost/asio.hpp>
namespace nucleona::parallel::asio_pool {

template<class T>
struct Future {

    DEFAULT_MOVE(Future)
    DISABLE_COPY(Future)

    Future(std::future<T>&& core, boost::asio::io_service& ios)
    : core_(std::move(core)) 
    , ios_(&ios)
    {}
    bool valid() const noexcept {
        return core_.valid();
    }
    decltype(auto) sync() {
        while(
            core_.wait_for(std::chrono::seconds(0)) != 
            std::future_status::ready
        ) {
            ios_->poll_one();
        }
        return core_.get();
    }
private:
    std::future<T> core_;
    boost::asio::io_service* ios_;
};

}
