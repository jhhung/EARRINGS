#pragma once

#include <iostream>
#include <vector>
#include <range/v3/all.hpp>

template <class DB>
class Stage
{
public:
    using OutputType = typename DB::Output;
    using InputType = typename DB::Input;
    Stage() = default;
    Stage(Stage&& stage) {}
    Stage(const Stage& stage) {}
    virtual OutputType run(InputType&&) const {}
};

template <class DB, class C_Table, class OCC_Table>
class Searcher
  : public Stage<DB>
{
    const C_Table& c_t_;
    const OCC_Table& occ_t_;
    using OutputType = typename DB::Output;
    using InputType = typename DB::Input;
public:
    explicit Searcher(
      const C_Table& c_t, 
      const OCC_Table& occ_t
    )
      : c_t_(c_t)
      , occ_t_(occ_t)
    {}

    OutputType run(InputType&& input) const override
    //std::pair<int, int> run(const std::string& s) override
    {
      const auto& seq = input.seq;
      //const auto& seq = input;
      auto iter = seq.rbegin();
      //std::pair<int, int> reg {0, input.size()};
      OutputType reg {0, seq.size()};
      while (iter != seq.rend())
      {
          int encode;
          switch(*iter)
          {
              case 'A':
                  encode = 0; break;
              case 'C':
                  encode = 1; break;
              case 'G':
                  encode = 2; break;
              case 'T':
                  encode = 3; break;
          }
        reg.first = c_t_[encode] + occ_t_[encode][reg.first];
        reg.second = c_t_[encode] + occ_t_[encode][reg.second];
        std::cout << *iter << " ";
        std::cout << "first: " << reg.first << ", second: " << reg.second << "\n";
        ++iter;
      }
      return reg;
	  }
};

/*template <Fastq>
class tailor_like
  : public ranges::view_facade<tailor_like<Fastq>>
{
  friend ranges::range_access;
  std::vector<Fastq> v_;
  Fastq& read() const { return *v_; }
  bool equal(range::default_sentinel) const { return }
public:
  tailor_like() = default;
  explicit tailor_like(std::vector<Fastq>&)

};
*/

/*template <class Rng, class IStream>
class Reader
  : public ranges::view_adaptor<transform_view<Rng, IStream>, Rng>
{
    friend ranges::range_access;
    IStream&& is_;
    class adaptor : public ranges_adaptor_base
    {
       IStream&& is_;
    public:
       adaptor() = default;
       adaptor(IStream&& is) : is_(is) {}
       void read(ranges::iterator_t<Rng> it) const
       {
           is_ >> *it;
       }
    };
    adaptor begin_adaptor() const ()
};*/

