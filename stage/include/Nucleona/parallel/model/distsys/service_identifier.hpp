#pragma once
#include <Nucleona/language.hpp>
#include <thread>
#include <Nucleona/parallel/model/distsys/traits/helpers.hpp>
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{

template<class CORE>
struct ServerIdentifierProto
{
  public:
    ServerIdentifierProto()
    : master_id_( traits::CoreHelper<CORE>::get_id() )
    {}

    bool is_master() const
    {
        return traits::CoreHelper<CORE>::get_id() == master_id_;
    }

  protected:
    const typename CORE::id master_id_;
};

template<class TRAIT>
struct ServiceIdentifier 
: public ServerIdentifierProto< typename TRAIT::Core >
{
    CREATE_DERIVED_TYPE_BODY(
          ServiceIdentifier
        , ServerIdentifierProto< 
            typename TRAIT::Core 
        > 
    );

};

}}}}
