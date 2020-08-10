#pragma once

#include <istream>
#include <range/v3/all.hpp>

namespace pipeline::range {

template <class Rng>
struct get_fasta_view
  : ranges::view_facade<get_fasta_view<Rng>, ranges::unknown>
{
private:
    friend ranges::range_access;
    std::istream* sin_;
    Rng* data_;
    
    struct cursor
    {
    private:
        get_fasta_view* rng_;
    public:
        cursor() = default;
        explicit cursor(get_fasta_view& rng)
          : rng_(&rng)
        {}
        void next()
        {
            rng_->next();
        }
        void read() const noexcept
        {
            return;
        }
        bool equal(ranges::default_sentinel) const
        {
            if constexpr (ranges::RandomAccessRange<decltype(*(rng_->data_))>())
            {
              //std::cout << "It's random range type!!!!!!!!!" << std::endl; 
              return !*rng_->sin_ || rng_->count_ > rng_->pool_size_;
            }
            else
              return !*rng_->sin_;
        }
    };

    void next()
    {
      typename std::remove_reference_t<decltype(*data_)>::key_type fa_;
      *sin_ >> fa_;
      
      if (!*sin_)
        return;
      if ((*data_).find(fa_) != (*data_).end())
        ++(*data_)[fa_];
      else
        (*data_)[fa_] = 1;
    }
    
    cursor begin_cursor()
    {
        return cursor{*this};
    }

public:
    get_fasta_view() = default;
    get_fasta_view(std::istream& sin, Rng& data)
      : sin_(&sin), data_(&data)
    {
        this->next();
    }
    void cached() noexcept
    {
        return;
    }
};

template <class Rng>
get_fasta_view(std::istream&, Rng&&) -> get_fasta_view<Rng>;

struct get_fasta_fn
{

    template <class Rng>
    auto operator()(std::istream& sin, Rng&& data) const
    {
        return get_fasta_view{sin, data};
    }
    
    template <class Rng>
    auto operator()(Rng&& data) const
    {
        return ranges::make_pipeable([&data](std::istream& sin){
            return get_fasta_view{sin, data};
        });
    }
};


RANGES_INLINE_VARIABLE(get_fasta_fn, get_fasta)

}
