#pragma once
#include <Nucleona/parallel/model/distsys/traits/basic_trait.hpp>
#include <Nucleona/parallel/model/distsys/task.hpp>
#include <Nucleona/parallel/model/distsys/config.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <Nucleona/language.hpp>
#include <memory>
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{

template<class TRAIT, class DERIVE>
class WorkerProto
{
protected:
    typedef TRAIT Trait;
    typedef WorkerProto<Trait, DERIVE> This;

public:
    typedef typename Trait::Core Core;
    typedef Task<Trait> TaskType;

protected:
    static const std::size_t BUFFER_SIZE = Config::worker_buffer_size;
    bool leave;
    bool running;
    // bool empty;
    std::unique_ptr<typename Trait:: template Queue<TaskType>> p_task_queue;

    std::unique_ptr<Core> worker_;
    static void do_job(TaskType& job )
    {
        job();
    }
    static void push_job( typename Trait:: template Queue<TaskType>& queue, TaskType& job )
    {
        int try_times = 0;
        while(!queue.push(job))   
        {
            if(try_times > 10000)
            {
                throw std::runtime_error( "worker crash : push job fail" );
            }
            try_times ++;
            std::this_thread::yield();
        }
    }    
public:
    inline static void inner_worker_content(TaskType& job, This* worker)
    {
        job();        
    }
    static std::function<void(void)> on_worker_load_;
    static std::function<void(void)> on_worker_leave_;
    static bool run_next_task( DERIVE* worker, bool& pop_state )
    {
        TaskType job;
        pop_state = worker->p_task_queue->pop(job);
        // TODO if another bug happened then add empty state check, 
        // here we assume boost spcs queue's empty is set after pop success
        if(pop_state) //non-block if fail then go else.
        {
            job();
        }
        else
        {
            std::this_thread::sleep_for(std::chrono::microseconds(3000));
        }
        return pop_state;

    }
    static void worker_content(DERIVE* worker)
    {
        // worker->empty = false;
        if(on_worker_load_)
        {
            on_worker_load_();
        }
        worker->running = true;
        while(!worker->p_task_queue->empty() || !worker->leave)
        {
            run_next_task( worker, worker->running );
        }
        if(on_worker_leave_)
        {
            on_worker_leave_();
        }
    }
    WorkerProto()
    : worker_       (nullptr)
    , p_task_queue  (new typename Trait:: template Queue<TaskType>( BUFFER_SIZE ))
    , leave         (false)
    {}
    DISABLE_COPY( WorkerProto );
    DEFAULT_MOVE( WorkerProto );
    bool is_worker_loaded()
    {
        return worker_ != nullptr;
    }
    void load_worker()
    {
        leave = false;
        worker_.reset(new Core(DERIVE::worker_content, (DERIVE*)this));
    }
    void push_job( TaskType& job )
    {
        push_job( *p_task_queue, job );
    }
    bool try_push_job( TaskType& job ) // true for success else fail
    {
        return p_task_queue->push(job);
    }
    void flush()
    {
        if(is_worker_loaded())
        {
            while(!p_task_queue->empty() || running)
            {
                std::this_thread::yield();
            }
        }
        else
        {
            return;
        }
    }
    void join()
    {
        flush();
        if(is_worker_loaded())
        {
            TaskType final_job(
                [this](){ leave = true; }
            );
            push_job(final_job);
            worker_->join();
            worker_.reset(nullptr);
        }
        else
        {
            return;
        }
    }
};
template<class TRAIT, class DERIVE>
std::function<void(void)> WorkerProto<TRAIT, DERIVE>::on_worker_load_;
template<class TRAIT, class DERIVE>
std::function<void(void)> WorkerProto<TRAIT, DERIVE>::on_worker_leave_;

template< class TRAIT>
class Worker : public WorkerProto<TRAIT, Worker<TRAIT> >
{
  public:
    CREATE_DERIVED_TYPE_BODY(Worker, WorkerProto<TRAIT COMMA Worker<TRAIT>> );
};


template<class CORE>
class Worker<traits::DependencyTrait<CORE>> 
: public WorkerProto<traits::DependencyTrait<CORE>, Worker<traits::DependencyTrait<CORE>>>
{
  public:
    using Base = WorkerProto<
          traits::DependencyTrait<CORE>
        , Worker<traits::DependencyTrait<CORE>>
    >;
    using Base::Base;
    using This = Worker<traits::DependencyTrait<CORE>>;
    typedef std::map<
          typename Base::Trait::TaskId
        , traits::get_future_type<typename Base::TaskType::TaskStore>
    > FutureMap;
    DISABLE_COPY( Worker );
    DEFAULT_MOVE( Worker );
    
    Worker(const FutureMap& fm)
    : Base()
    , p_repost_queue( new typename Base::Trait::template Queue<typename Base::TaskType> ( 
        Base::BUFFER_SIZE 
    ) )
    , task_tree( fm )
    {}
    static void worker_content(This* worker)
    {
        if(Base::on_worker_load_)
        {
            Base::on_worker_load_();
        }
        while(true && !worker->leave)
        {
            typename Base::TaskType job;
            if(worker->p_task_queue->pop(job)) //non-block if fail then go else.
            {
                // if ( job.has_dep( worker->task_tree ) )
                // {
                //     Base::push_job( *(worker->p_repost_queue), job );
                //     continue;
                // }
                // else
                // {
                //     job();
                // }
                job();

            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
        }
        if(Base::on_worker_leave_)
        {
            Base::on_worker_leave_();
        }
    }
    
    std::unique_ptr<typename Base::Trait::template Queue<typename Base::TaskType>> p_repost_queue;
    const FutureMap& task_tree;
};



//===============boost thread version===========//
//#include <boost/chrono.hpp>
//template<>
//class Worker<BasicTrait<boost::thread>> 
//    : public WorkerProto<BasicTrait<boost::thread>> 
//{
//protected:
//    typedef WorkerProto<BasicTrait<boost::thread>> _Base;
//    
//public:
//    //using _Base::_Base;
//    bool try_join()
//    {
////        return worker_->try_join_for(boost::chrono::milliseconds(ParallelPoolConstDef::boost_join_time_out));
//    }
//
//    void join()
//    {
//    }
//};
}}}}
