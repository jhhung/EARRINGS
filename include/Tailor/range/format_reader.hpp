// format_reader.hpp
#pragma once

#include <istream>
#include <mutex>
#include <Biovoltron/format/fastq.hpp>
#include <range/v3/all.hpp>

namespace pipeline::range {

template <class Val = biovoltron::format::FASTQ<std::string, std::string>>
struct format_reader
  : ranges::view_facade<format_reader<Val>, ranges::unknown>
{
private:
    static std::mutex input_mutex;
    friend ranges::range_access;
    static constexpr auto pool_size_ = 10001;
    std::size_t count_;
    std::istream* sin_;
    //ranges::semiregular_box_t<Val> obj_;
    Val obj_;
    
    struct cursor
    {
    private:
        format_reader* rng_ = nullptr;
    public:
        cursor() = default;
        explicit cursor(format_reader* rng)
          : rng_(rng)
        {}
        void next()
        {
            rng_->next();
        }
        Val& read() const noexcept
        {
            return rng_->cached();
        }
        bool equal(ranges::default_sentinel) const
        {
            return !rng_->sin_ || rng_->count_ == format_reader::pool_size_;
        }
        bool equal(cursor that) const
        {
            return !rng_->sin_ == !that.rng_->sin_;
        }
    };

    void next()
    {
      std::scoped_lock<std::mutex> sl(this->input_mutex);      
      if(!(*sin_ >> cached()))
        sin_ = nullptr;
      ++count_;
    }
    
    cursor begin_cursor()
    {
        return cursor{this};
    }

public:
    format_reader() = default;
    format_reader(std::istream& sin)
      : sin_(&sin), obj_{}, count_{0}
    {
        this->next();
    }
    Val& cached() noexcept
    {
        return obj_;
    }
};

/*
template <class Rng>
read_fastq_view(std::istream&, Rng&&) -> read_fastq_view<Rng>;
*/

template<class Val>
struct format_reader_fn
{
    auto operator()(std::istream& sin) const
    {
        return format_reader<Val>{sin};
    }
    
    auto operator()() const
    {
        return ranges::make_pipeable([](std::istream& sin){
            return format_reader<Val>{sin};
        });
    }
};

template <class Val>
std::mutex format_reader<Val>::input_mutex;

//RANGES_INLINE_VARIABLE(format_reader_fn, format_reader)

}
