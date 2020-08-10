#pragma once

#include <range/v3/all.hpp>
#include <boost/hana/equal.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/string.hpp>
#include <boost/hana/for_each.hpp>

namespace pipeline::range {
    
  constexpr auto filter_list = boost::hana::make_tuple(
    boost::hana::string_c<'m', 'i', 'R', 'N', 'A'>
  , boost::hana::string_c<'m', 'i', 'r', 't', 'r', 'o', 'n'>
  , boost::hana::string_c<'t', 'R', 'N', 'A'>
  , boost::hana::string_c<'u','n','p','r','o','c','e','s','s','e','d','_',
        'p','s','e','u','d','o','g','e','n','e'>
  );

  template <class AnnoBed>
  constexpr void match(AnnoBed&& input)
  {
    bool res = false;
    boost::hana::for_each(filter_list, [&](auto&& str) {
        for (auto&& anno_type : input.anno_type_) {
            res |= (anno_type == str.c_str());
        }
    });
    if (res) 
    {
        input.is_filtered = 1;
    }
    return;
  }

  auto filter()
  {
    return ranges::view::transform (
      [=](auto&& anno_bed_adaptor) {
        //auto& anno_bed = *anno_bed_adaptor.get_data_ptr();
        match(*anno_bed_adaptor.ptr_);
        return anno_bed_adaptor.get_data_ptr();
    });
  }  
    
}
