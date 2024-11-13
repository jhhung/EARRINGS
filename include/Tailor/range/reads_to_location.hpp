#pragma once

#include <range/v3/all.hpp>
#include <OldBiovoltron/annotator/annotator.hpp>

namespace pipeline::range {

auto convert_to_location()
{
  return ranges::view::transform (
    [](auto&& reads) {
      if (reads.flag_ != 0 && reads.flag_ != 16)
        return std::vector<biovoltron::annotator::Annotation>{};
      
      std::vector<biovoltron::annotator::Annotation> res;
      for (int i = 0; i < reads.real_pos_.size(); ++i)
      {
        res.emplace_back (
          biovoltron::annotator::Annotation {
            reads.chrs_[i], 
            reads.real_pos_[i] - 1,
            reads.real_pos_[i] + static_cast<std::uint32_t>(reads.fq_.seq.size()) - 1
        });
      }
      return res;
  });
}
}
