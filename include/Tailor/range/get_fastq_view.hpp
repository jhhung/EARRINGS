#pragma once

#include <istream>
#include <range/v3/all.hpp>

namespace pipeline::range {

template <class Rng>
struct get_fastq_view
  : ranges::view_facade<get_fastq_view<Rng>, ranges::unknown>
{
private:
    friend ranges::range_access;
    static constexpr auto pool_size_ = 10000;
    std::istream* sin_;
    Rng* data_;
    int count_;
    
    struct cursor
    {
    private:
        get_fastq_view* rng_;
    public:
        cursor() = default;
        explicit cursor(get_fastq_view& rng)
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
        if constexpr (ranges::RandomAccessRange<Rng>())
        {
          typename std::remove_reference_t<decltype(*data_)>::value_type fq_;
          *sin_ >> fq_;
          ++count_;
          
          //if (!*sin_)
          //std::cout << "count: " << count_ << std::endl;
          //std::cout << "if true ? " << (count_ > pool_size_) << std::endl;
          if (!*sin_ || count_ > pool_size_)
            return;

          (*data_).emplace_back(fq_);  
        }
        else
        {
          typename std::remove_reference_t<decltype(*data_)>::key_type fq_;
          *sin_ >> fq_;
          
          if (!*sin_)
            return;
  
          if ((*data_).find(fq_) != (*data_).end())
            ++(*data_)[fq_];
          else
            (*data_)[fq_] = 1;
        }
    }
    
    cursor begin_cursor()
    {
        return cursor{*this};
    }

public:
    get_fastq_view() = default;
    get_fastq_view(std::istream& sin, Rng& data)
      : sin_(&sin), data_(&data), count_(0)
    {
        this->next();
    }
    void cached() noexcept
    {
        return;
    }
};

template <class Rng>
get_fastq_view(std::istream&, Rng&&) -> get_fastq_view<Rng>;

struct get_fastq_fn
{

    template <class Rng>
    auto operator()(std::istream& sin, Rng&& data) const
    {
        return get_fastq_view{sin, data};
    }
    
    template <class Rng>
    auto operator()(Rng&& data) const
    {
        return ranges::make_pipeable([&data](std::istream& sin){
            return get_fastq_view{sin, data};
        });
    }
};


RANGES_INLINE_VARIABLE(get_fastq_fn, get_fastq)

}
