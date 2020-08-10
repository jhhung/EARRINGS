#pragma once

namespace tailor {
struct TailParas
{
  const std::size_t n_thread;
  const uint32_t min_prefix_len;
  const uint32_t min_multi;
  const bool allow_mm;
  //const bool allow_t2c;
  
  explicit TailParas(const std::size_t n_thread, const uint32_t min_prefix_len, const uint32_t min_multi, const bool allow_mm)
    : n_thread(n_thread), min_prefix_len(min_prefix_len), min_multi(min_multi), allow_mm(allow_mm)
  {}

  void print() const
  {
    // std::cout << "nthread: " << n_thread << std::endl;
    // std::cout << "minimal prefix length: " << min_prefix_len << std::endl;
    // std::cout << "multi: " << min_multi << std::endl;
    // std::cout << "allow_mm: " << allow_mm << std::endl;
    //std::cout << "allow_t2c: " << allow_t2c << std::endl;
  }
};
}
