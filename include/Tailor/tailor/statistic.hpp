#pragma once

namespace tailor {
  
  template <class Query, class OStream>
  void anno_stat(Query&& anno_bed_pv, OStream&& out)
  {
    double ppm_rate = 0.0;
    std::map<std::string, double> annoTypeStat;
    double total = 0.0;
    for (size_t idx = 0; idx < anno_bed_pv.size(); ++idx)
    {
      auto& bed_anno = *anno_bed_pv[idx].ptr_;
      double reads_count = bed_anno.reads_count_;
      auto& annoType = bed_anno.anno_type_;
      ppm_rate += reads_count;
      
      auto anno_count = annoType.size();
      
      for (auto&& i : annoType)
      {
        if (annoTypeStat.find(i) == annoTypeStat.end()) 
        {
          annoTypeStat[i] = 0.0;
        }
        annoTypeStat[i] += reads_count / anno_count; 
        total += reads_count / anno_count; 
      }
    }

    ppm_rate = 1000000.0 / ppm_rate;

    for (size_t idx = 0; idx < anno_bed_pv.size(); ++idx)
    {
      auto& bed_anno = *anno_bed_pv[idx].ptr_;
      bed_anno.PPM_ = bed_anno.reads_count_ * ppm_rate;
    }

    // output
    for (auto&& i : annoTypeStat)
    {
        out << i.first << "\t" << std::fixed << std::setprecision(0) << i.second << "\n";
    }
    out << "Total" << "\t" << std::fixed << std::setprecision(0) << total << "\n";
  }
}
