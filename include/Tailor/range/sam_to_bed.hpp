#pragma once

#include <range/v3/all.hpp>
#include <mutex>

namespace pipeline::range {

std::mutex push_mutex;

template <class AnnoBed_v>
auto sam_to_bed(AnnoBed_v&& anno_bed_v)
{
  return ranges::view::transform (
    [&anno_bed_v](auto&& AlignedReads) {
      using AnnoBed = 
        typename std::remove_reference_t<decltype(anno_bed_v)>::value_type;
      using SEQ = typename std::remove_reference_t<decltype(AlignedReads)>::seq_type;
      using IndexType = typename std::remove_reference_t<decltype(AlignedReads)>::index_type;

      const auto& AR = AlignedReads;

      std::scoped_lock<std::mutex> sl(push_mutex);
      for (int i = 0; i < AR.chrs_.size(); ++i) 
      {
        anno_bed_v.emplace_back(AnnoBed {
            std::move(AR.chrs_[i])
          , std::move(AR.real_pos_[i])
          , AR.real_pos_[i] + AR.fq_.seq.size() - 1
          , AR.flag_ == 0 ? "+" : "-"
          , std::move(AR.NH_num_)
          , std::move(AR.reads_count_)
          , AR.reads_count_ / AR.NH_num_ // watch out the type conversion
          , 0.0
          , AR.fq_.seq.size()
          , std::move(AR.tail_pos_)
          , std::move(AR.fq_.seq)
          , std::move(AR.tail_seq_)
          , std::move(AR.MM_tag_)
          , ""
          , ""
        });
      }
      return 0;
  });
} 

template <class AnnoBed_v>
auto sam_to_bed_set(AnnoBed_v&& anno_bed_v)
{
  return ranges::view::transform (
    [&anno_bed_v](auto&& AlignedReads_v) {
      using AnnoBed = 
        typename std::remove_reference_t<decltype(anno_bed_v)>::value_type;
      using AlignedReads = typename std::remove_reference_t<decltype(AlignedReads_v)>::value_type;
      using SEQ = typename AlignedReads::seq_type;
      using IndexType = typename AlignedReads::index_type;

      for (auto&& AR : AlignedReads_v)
      {
        std::scoped_lock<std::mutex> sl(push_mutex);
        for (int i = 0; i < AR.chrs_.size(); ++i) 
        {
          anno_bed_v.emplace_back(AnnoBed {
              std::move(AR.chrs_[i])
            , std::move(AR.real_pos_[i])
            , AR.real_pos_[i] + (std::uint32_t)AR.fq_.seq.size() - 1
            , AR.flag_ == 0 ? "+" : "-"
            , std::move(AR.NH_num_)
            , std::move(AR.reads_count_)
            , (double)AR.reads_count_ / (double)AR.NH_num_ // watch out the type conversion
            , 0.0
            , AR.fq_.seq.size()
            , std::move(AR.tail_pos_)
            , std::move(AR.fq_.seq)
            , std::move(AR.tail_seq_)
            , std::move(AR.MM_tag_)
            , std::vector<std::string> {}
            , std::vector<std::string> {}
            , false
          });
        }
      }
      return 0;
  });   
}

template <class AnnoBed_v>
auto sam_to_anno_bed_set(AnnoBed_v&& anno_bed_pv)
{
  return ranges::view::transform (
    [&anno_bed_pv](auto&& AlignedReads_v) {
      using AnnoBedAdaptor = 
        typename std::remove_reference_t<decltype(anno_bed_pv)>::value_type;
      using AnnoBed = typename AnnoBedAdaptor::AnnoBedType;
      using AlignedReads = typename std::remove_reference_t<decltype(AlignedReads_v)>::value_type;
      using SEQ = typename AlignedReads::seq_type;
      using IndexType = typename AlignedReads::index_type;

      for (auto&& AR : AlignedReads_v)
      {
        std::scoped_lock<std::mutex> sl(push_mutex);
        for (int i = 0; i < AR.chrs_.size(); ++i) 
        {
          anno_bed_pv.emplace_back( AnnoBedAdaptor ( AnnoBed {
              std::move(AR.chrs_[i])
            , std::move(AR.real_pos_[i])
            , AR.real_pos_[i] + (std::uint32_t)AR.fq_.seq.size() - 1
            , AR.flag_ == 0 ? "+" : "-"
            , std::move(AR.NH_num_)
            , std::move(AR.reads_count_)
            , (double)AR.reads_count_ / (double)AR.NH_num_ // watch out the type conversion
            , 0.0
            , AR.fq_.seq.size()
            , std::move(AR.tail_pos_)
            , std::move(AR.fq_.seq)
            , std::move(AR.tail_seq_)
            , std::move(AR.MM_tag_)
            , std::vector<std::string> {}
            , std::vector<std::string> {}
            , false
          })
          );
        }
      }
      return 0;
  });   
}


}
