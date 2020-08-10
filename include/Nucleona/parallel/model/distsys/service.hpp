#pragma once
#include <Nucleona/parallel/model/distsys/traits/basic_trait.hpp>
#include <Nucleona/parallel/model/distsys/worker.hpp>
#include <Nucleona/language.hpp>
#include <Nucleona/parallel/model/distsys/service_identifier.hpp>
#include <Nucleona/parallel/model/distsys/future.hpp>
// using namespace std::chrono_literals;
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{
struct IdRoundIterator
{
    IdRoundIterator(size_t max)
    : max_(max)
    {}

    IdRoundIterator() = default;

    size_t mod_max( size_t& n )
    {
        return n % max_;
    }
    size_t get_round_next()
    {
        std::size_t tmp = p_next_id_->load(std::memory_order_acquire );
        while(!p_next_id_->compare_exchange_weak(
            tmp, tmp + 1, std::memory_order_release, std::memory_order_relaxed ) );
        return mod_max(tmp);
    }
    std::tuple<size_t, size_t> get_next()
    {
        std::size_t tmp = p_next_id_->load(std::memory_order_acquire );
        while(!p_next_id_->compare_exchange_weak(
            tmp, tmp + 1, std::memory_order_release, std::memory_order_relaxed ) );
        return std::make_tuple(tmp, mod_max(tmp));
    }
    void reset(size_t n)
    {
        max_ = n;
    }
    size_t max_;
    std::unique_ptr<std::atomic<std::size_t>> p_next_id_ { new std::atomic<std::size_t>(0) };
};
template<class TRAIT>
class ServiceProto
{
    typedef TRAIT Trait;
public:
    typedef Worker<Trait>                   WorkerType;
    typedef Task<Trait>                     TaskType;
protected:
    std::vector<WorkerType> worker_group_;
    std::size_t pool_size_;

    IdRoundIterator counter_;
    // std::unique_ptr<std::mutex> p_counter_mux_ { new std::mutex };

    std::map<typename Trait::Core::id, WorkerType*> core_map_worker_;
    // std::unique_ptr<std::mutex> p_core_map_worker_mux_ { new std::mutex };


    bool try_job_post ( TaskType& task, WorkerType& worker )// true for success else fail
    {
        if(pool_size_ > 1)
        {
            if(!worker.is_worker_loaded())
            {
                worker.load_worker();
            }
            return worker.try_push_job(task);
        }
        else
        {
            WorkerType::inner_worker_content(task, &worker);
            return true;
        }
    }
    template<class R>
    auto pack_task_post ( std::packaged_task<R(void)>& task_pack, const size_t& worker_id )
    {
        auto future = task_pack.get_future();
        TaskType task ( std::move( task_pack ) );
        bool post_fail = true;
        auto& worker = worker_group_[worker_id];
        while(post_fail)
        {
            if(try_job_post(task, worker_group_[worker_id]))
            {
                post_fail = false;
            }
            std::this_thread::yield();
        }
        return make_future(std::move(future), core_map_worker_ );
    }
    void join_pool()
    {
        if(pool_size_ > 1 )
        {
            for( size_t i = 0; i < worker_group_.size(); i++ )
            {
                worker_group_[i].join();
            }
        }
    }
    template<class F>
    auto for_each_worker_post( F&& task )
    {
        using Ret = decltype(task());
        using Fut = Future<Ret, Trait>;
        std::vector<Fut> check_points;
        for ( size_t i (0); i < pool_size_; i++ )
        {
            std::packaged_task<Ret(void)> tmp( nucleona::language::copy( task ) );
            check_points.emplace_back( pack_task_post( tmp, i ) );
        }
        return check_points;
    }

public:
    template<class... W_ARGS>
    ServiceProto(std::size_t pool_size, W_ARGS&&... w_args )
    : pool_size_        ( pool_size )
    , counter_          ( pool_size )
    {
        for (std::size_t i (0); i < pool_size; i ++ )
        {
            worker_group_.emplace_back(FWD(w_args)...);
        }
        for( size_t i (0); i < pool_size_; i++ )
        {
            std::packaged_task<void(void)> task_pack( 
                [this,i]()
                {
                    // std::lock_guard<std::mutex>(*p_core_map_worker_mux_);
                    core_map_worker_.emplace( std::this_thread::get_id(), &worker_group_.at(i) );
                }
            );
            auto fut = pack_task_post( task_pack, i );
            fut.sync();
        };
    }
    ServiceProto()
    {}
    DEFAULT_MOVE( ServiceProto )

    size_t get_parallization()
    {
        return pool_size_;
    }
    template<class F> 
    auto job_post( F&& func, const size_t& worker_id )
    {
        std::packaged_task<decltype(func())()> task_pack( FWD(func) );
        return pack_task_post(task_pack, worker_id );
    }

    template<class F> 
    auto job_post( F&& func )
    {
        std::packaged_task<decltype(func())()> task_pack( FWD(func) );
        auto id ( counter_.get_next() );
        size_t worker_id = std::get<1>(id);
        return pack_task_post(task_pack, worker_id );
    }
    void flush_tasks()
    {
        if(pool_size_ > 1 )
        {
            for( size_t i (pool_size_); i > 0; i -- )
            {
                worker_group_.at(i - 1).flush();
            }
        }
    }

    void resize_pool(std::size_t pool_size)
    {
        join_pool();
        pool_size_ = pool_size;
        counter_.reset(pool_size);
        std::vector<WorkerType> tmp(pool_size);
        worker_group_.swap(tmp);
    }
    inline bool support_owner_as_worker() const {
        return false;
    }
    ~ServiceProto()
    {
        join_pool();
    }

};    
template<class TRAIT>
class Service
: public ServiceProto<TRAIT>
{
  public:
    using Base = ServiceProto<TRAIT>;
    using Base::Base;
    DEFAULT_COPY( Service );
    DEFAULT_MOVE( Service );
    Service() = default;
};

template<class CORE>
class Service<traits::DependencyTrait<CORE>> 
: public ServiceProto<traits::DependencyTrait<CORE>>
{
  public:
    using Base = ServiceProto<traits::DependencyTrait<CORE>>;
    using Base::Base;
    Service( std::size_t n )
    : Base( n, Base::job_tree_ )
    {}
    DEFAULT_COPY( Service )
    DEFAULT_MOVE( Service )
    template<class F> 
    decltype(auto) job_post( F&& task_content, const std::vector<std::size_t>& dep )
    {
        auto task = [dep, job = std::forward<F>(task_content), this]()
        {
            for( auto& i : dep )
            {
                this->flush_one(i);
            }
            job();
        };
        return Base::job_post( std::move(task) );
    }
    template<class T>
    decltype(auto) job_post( T&& arg )
    {
        return Base::job_post( FWD(arg) );
    }
    // void resize_pool( std::size_t pool_size )
    // {
    //     Base::join_pool();
    //     Base::pool_size_ = pool_size;
    //     Base::counter_.reset(pool_size);
    //     std::vector<typename Base::WorkerType> tmp;
    //     for ( std::size_t i(0); i < pool_size; i ++ )
    //     {
    //         tmp.emplace_back( Base::job_tree_ );
    //     }
    //     Base::worker_group_.swap(tmp);
    // }
    
};
}}}}
