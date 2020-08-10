#pragma once

#include <istream>
#include <range/v3/all.hpp>

namespace pipeline::range {

template <class Stream, class Fastq>
auto reader()
{
    return ranges::view::transform (
        [](auto&& rng) {
            Stream buf;
            for (auto&& str : rng)
            {
                buf << str << "\n";
            }
            Fastq res;
            Fastq::get_obj(buf, res);
            return res;
        }
    );
}

template <class IStream, class Stream, class Fastq>
auto fastq_reader(IStream&& input)
{
     return ranges::getlines_range(input) 
          | ranges::view::chunk(4) 
          | reader<Stream, Fastq>();
}

template <class Fastq>
struct get_fastq_view
  : ranges::view_facade<get_fastq_view<Fastq>, ranges::unknown>
{
private:
    friend ranges::range_access;
    std::istream* sin_;
    Fastq fq_;
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
        Fastq& read() const noexcept
        {
            return rng_->fq_;
        }
        bool equal(ranges::default_sentinel) const
        {
            return !*rng_->sin_;
        }
    };
    void next()
    {
        *sin_ >> fq_;
    }
    cursor begin_cursor()
    {
        return cursor{*this};
    }
public:
    get_fastq_view() = default;
    get_fastq_view(std::istream& sin)
      : sin_(&sin), fq_{}
    {
        this->next();
    }
    Fastq& cached() noexcept
    {
        return fq_;
    }
};

template <class Fastq>
struct get_fastq_fn
{
    //template <class Fastq>
    get_fastq_view<Fastq> operator()(std::istream& sin) const
    {
        return get_fastq_view<Fastq>{sin};
    }
    
    //template <class Fastq>
    auto operator()() const
    {
        return ranges::make_pipeable([](std::istream& sin){
            //return (*__this)<Fastq>(sin);
            return get_fastq_view<Fastq>{sin};
        });
    }
};

//RANGES_INLINE_VARIABLE(get_fastq_fn, get_fastq)

}
