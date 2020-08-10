#pragma once
#include <future>
#include <Nucleona/language.hpp>
#include <Nucleona/parallel/model/distsys/service_identifier.hpp>
#include <Nucleona/parallel/model/distsys/traits/helpers.hpp>
#include <map>
#include <Nucleona/parallel/model/distsys/worker.hpp> 
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{

// template<class TRAIT> class Worker; 
template<class T>
using FutureProto = std::future<T>;
template<class T, class TRAIT>
struct Future : public FutureProto<T>
{
    using Base = FutureProto<T>;
    using WorkerType = Worker<TRAIT>;
    using CoreMapWorker = std::map<typename TRAIT::Core::id, Worker<TRAIT>*>;
    Future( 
          std::future<T>&& f
        , CoreMapWorker& cmw
    )
    : Base( std::move( f ) )
    , core_map_worker_( cmw )
    {}
    Future( 
        CoreMapWorker& cmw
    )
    : Base()
    , core_map_worker_( cmw )
    {}
    DISABLE_COPY( Future );

    Future( Future&& fut )
    : Base( std::move( (Base&) fut ) )
    , core_map_worker_( fut.core_map_worker_ )
    {}

    Future& operator=( Future&& fut )
    {
        if( &core_map_worker_ == &(fut.core_map_worker_) )
            (Base&)*this = std::move( (Base&)fut );
        else
            throw std::runtime_error("cannot move bettween 2 future generate by different core_map_worker_");
    }
    
    T sync()
    {
        auto this_core_id = traits::CoreHelper<typename TRAIT::Core>::get_id();
        auto itr = core_map_worker_.find( this_core_id );
        if( itr == core_map_worker_.end() )
        {
            Base::wait();
        }
        else
        {
            auto& worker = *(itr->second);
            while( Base::wait_for(std::chrono::milliseconds(1) ) != std::future_status::ready )
            {
                // Task<TRAIT> job;
                bool b;
                WorkerType::run_next_task( &worker, b);
            }
        }
        return Base::get();
    }
  protected:
    void wait() const
    {
        return Base::wait();
    }
    CoreMapWorker& core_map_worker_;
};

template<class T, class TRAIT>
auto make_future( std::future<T>&& f, std::map<typename TRAIT::Core::id, Worker<TRAIT>*>& core_map_worker)
{
    return Future<T, TRAIT>( std::move(f), core_map_worker );
}



}}}}
