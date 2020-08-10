#pragma once

#include <Biovoltron/annotator/annotator.hpp>

namespace tailor {
  
  template <class Query, class Seq, class IndexType>
  void annotate(const std::string& db_path, Query&& anno_bed_pv)
  {
    std::ifstream anno_input(db_path);
    biovoltron::annotator::BedSixAdaptor bed_file_adaptor(anno_input);
    auto db = std::move(bed_file_adaptor());
    
    biovoltron::annotator::Points_Sweep<
      biovoltron::annotator::Adaptor<biovoltron::format::BedSix>,
      biovoltron::annotator::AnnoBedAdaptor<Seq, IndexType>
      > points_sweep;
    
    points_sweep.set_db(db);
    points_sweep.set_query(anno_bed_pv);

    auto anno_results = points_sweep.annotate();
    for (size_t idx = 0; idx < anno_results.size(); ++idx)
    {
      auto bed_adaptor = std::get<0>(anno_results[idx]);
      auto& bed_anno = *anno_bed_pv[idx].ptr_;
      for (auto&& i : bed_adaptor)
      {
        bed_anno.anno_type_.emplace_back(i->get_gene_type());
        bed_anno.anno_seed_.emplace_back(i->get_anno_seed());
      }
    }      
  }
}
