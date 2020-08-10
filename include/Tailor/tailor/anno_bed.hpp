#pragma once
#include <iostream>

namespace tailor {

template <class SEQ, class IndexType>
struct AnnoBed
{
  std::string chr_;
  IndexType start_;
  IndexType end_;
  std::string strand_;
  IndexType aligned_count_;
  IndexType raw_count_;
  double reads_count_;
  double PPM_;
  std::size_t length_;
  int tail_len_;
  SEQ seq_;
  SEQ tail_;
  std::string MM_tag_;
  std::vector<std::string> anno_type_;
  std::vector<std::string> anno_seed_;
  bool is_filtered;
};

}
