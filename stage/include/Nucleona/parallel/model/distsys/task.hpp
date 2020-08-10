#pragma once
#include <set>
#include <map>
#include <Nucleona/parallel/model/distsys/traits/basic_trait.hpp>
#include <Nucleona/parallel/model/distsys/traits/dependency_trait.hpp>
#include <Nucleona/language.hpp>
#include <Nucleona/parallel/model/distsys/traits/helpers.hpp>
#include <Nucleona/parallel/model/distsys/service_identifier.hpp>
#include <Nucleona/pointer.hpp>
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{
template<class TRAIT>
class Task
{};
template<class TRAIT>
class TaskProto
{
    struct VCall
    {
        virtual void call() = 0;
    };
    template<class T>
    struct CallWrapped : public VCall
    {
        T func;
        CallWrapped( T&& o )
        : func ( FWD(o) )
        {}
        virtual void call() override
        {
            func();
        }
    };
protected:
    typedef TRAIT Trait;
protected:
    std::unique_ptr<VCall> storage_;
    // std::function<void(void)> storage_;
public:
    TaskProto() = default;
    DEFAULT_MOVE( TaskProto )
    // DISABLE_COPY( Task ), workaround for boost spsc_queue
    
    TaskProto( const TaskProto& task )
    : storage_( std::move ( const_cast<TaskProto&>( task ).storage_ ) )
    {}
    TaskProto& operator= ( const TaskProto& task )
    {
        *this = std::move( const_cast<TaskProto&>(task) );
        return *this;
    }
    template<class FUNC>
    TaskProto(FUNC&& packf )
    : storage_( nucleona::mkuniq( new CallWrapped<FUNC>(FWD(packf))) )
    // : storage_( [ func = std::move(packf) ]() mutable { func(); } )
    {}

    void operator()()
    {
        storage_->call();
        // storage_();
    }
};
template<class CORE>
class Task<traits::BasicTrait<CORE>>
: public TaskProto<traits::BasicTrait<CORE>>
{ 
  public:
    CREATE_DERIVED_TYPE_BODY( Task, TaskProto<traits::BasicTrait<CORE>> );
};

template<class CORE>
class Task<traits::DependencyTrait<CORE>>
: public TaskProto<traits::DependencyTrait<CORE>>
{ 
  public:
    CREATE_DERIVED_TYPE_BODY( Task, TaskProto<traits::DependencyTrait<CORE>> );
  private:
    typedef Task<typename Base::Trait> This;
  public:
    typedef typename Base::Trait::TaskStore TaskStore;
    typedef typename Base::Trait::TaskId    TaskId;
    bool has_dep( const std::map<TaskId, traits::get_future_type<TaskStore>>& history_tasks )
    {
        update_dep( history_tasks );
        if ( upstream_deps_.size() == 0 )
            return false;
        else return true;
    }
  private:
    void update_dep( const std::map<TaskId, traits::get_future_type<TaskStore>>& history_tasks )
    {
        std::vector<TaskId> tasks_not_found;
        for ( auto&& up_task_id : upstream_deps_ )
        {
            auto itr = history_tasks.find( up_task_id );
            if ( itr == history_tasks.end() )
            {
                tasks_not_found.push_back( up_task_id );
            }
            else
            {
                if( itr->second.valid() )
                {
                    if ( itr->wait_for( std::chrono::milliseconds(0) ) == std::future_status::ready )
                    {
                        tasks_not_found.push_back( up_task_id );
                    }
                }
                else
                {
                    tasks_not_found.push_back( up_task_id );
                }
            }
        }
        if ( tasks_not_found.size() != upstream_deps_.size() )
        {
            for ( auto&& miss_task : tasks_not_found )
            {
                upstream_deps_.erase( miss_task );
            }
        }
        else
        {
            upstream_deps_.clear();
        }
    }
    std::set< TaskId > upstream_deps_;
};
/* This is a compiler test 
typedef Task<BasicTrait<boost::thread>> BoostPj;
void foo()
{}
void compile_test()
{
    BoostPj boost_pj([](){}); // pass in g++ fail for clang++
}
*/
}}}}
