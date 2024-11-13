#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <type_traits>
#include <charconv>
#include <OldBiovoltron/format/gff.hpp>
#include <Tailor/tailor/anno_bed.hpp>
#include <charconv>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>


namespace biovoltron::format
{
    class BedSix
    {
    public:
      using BedType = 
        std::tuple<
          std::string
        , std::uint32_t
        , std::uint32_t
        , std::string
        , std::string
        , std::string
      >;
      
      BedType data;
      
      template<int N>
      auto& get() const
      {
				return std::get<N>(data);
      }
      
      template <int N>
      static void getbed(const std::vector<std::string>& content, BedType& b)
      {
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(std::get<N>(b))>, std::uint32_t>)
          std::from_chars(content[N].data(), content[N].data() + content[N].size(), std::get<N>(b));
        else
          std::get<N>(b) = content[N];
        if constexpr (N != 5)
          getbed<N+1>(content, b);
      }
       
      static std::istream& get_obj(std::istream& is, BedSix& bed)
      {
        std::string str;
        std::getline(is, str);
        if (!is.good())
          return is;
        std::vector<std::string> content;
        boost::iter_split(content, str, boost::algorithm::first_finder("\t"));
        getbed<0>(content, bed.data);

        return is;
      }
    };
}

namespace biovoltron::annotator 
{
    template<class FORMAT>
    class Adaptor {};
    /*
    template<>
    class Adaptor<biovoltron::format::VCF>
    {
    public:
        Adaptor(biovoltron::format::VCF&& vcf)
            : ptr_(std::make_unique<biovoltron::format::VCF>(vcf)) {}

        std::unique_ptr<biovoltron::format::VCF> 
        get_data_ptr() { return std::move(ptr_); }
        std::string get_chrom() { return ptr_->chrom; }
        // vcf starts from 1
        uint32_t get_start() { return ptr_->pos - 1; }
        uint32_t get_end() { return ptr_->pos + ptr_->ref.length() - 2; }

    private:
        std::unique_ptr<biovoltron::format::VCF> ptr_;
    };*/

    template<>
    class Adaptor<biovoltron::format::GFF>
    {
    public:
        // chromosome names can be given with or without the 'chr' prefix
        Adaptor(biovoltron::format::GFF&& gff)
            : ptr_(std::make_unique<biovoltron::format::GFF>(gff)) {}

        std::unique_ptr<biovoltron::format::GFF> 
        get_data_ptr() { return std::move(ptr_); }

        std::string get_chrom() 
        { 
            if (ptr_->get<0>().substr(0, 3) == "chr") 
                return ptr_->get<0>().substr(3);
            else 
                return ptr_->get<0>(); 
        }  // format: chr1, chr2...
        
        // gff sequence numbering starting at 1
        uint32_t get_start() { return std::stoi(ptr_->get<3>()) - 1; }
        uint32_t get_end() { return std::stoi(ptr_->get<4>()) - 1; }

    private:
        std::unique_ptr<biovoltron::format::GFF> ptr_;
    };

    template<>
    class Adaptor<biovoltron::format::BedSix>
    {
    public:
        Adaptor(biovoltron::format::BedSix&& bed_six)
            : ptr_(std::make_unique<biovoltron::format::BedSix>(bed_six)) {}

        std::unique_ptr<biovoltron::format::BedSix> get_data_ptr() { return std::move(ptr_); }

        std::string get_chrom()
        {
            if (ptr_->get<0>().substr(0, 3) == "chr") 
                return ptr_->get<0>().substr(3, ptr_->get<0>().find(" ") - 3);
            else 
                return ptr_->get<0>();  
        }

        uint32_t get_start() { return ptr_->get<1>() - 1; }
        uint32_t get_end() { return ptr_->get<2>() - 1; }
        std::string get_strand() { return ptr_->get<3>(); }
        std::string get_gene_type() { return std::move(ptr_->get<4>()); }
        std::string get_anno_seed() { return std::move(ptr_->get<5>()); }
    
    private:
        std::unique_ptr<biovoltron::format::BedSix> ptr_;
    };

    template<class SEQ, class IndexType>
    class AnnoBedAdaptor
    {
    public:
        using AnnoBedType = tailor::AnnoBed<SEQ, IndexType>;
        AnnoBedAdaptor(AnnoBedType&& anno_bed)
            : ptr_(std::make_unique<AnnoBedType>(anno_bed)) {}

        std::unique_ptr<AnnoBedType> get_data_ptr() 
        { 
            return std::move(ptr_); 
        }

        std::string get_chrom() 
        { 
            if (ptr_->chr_.substr(0, 3) == "chr")
              return ptr_->chr_.substr(3, ptr_->chr_.find(" ") - 3);
            else
              return ptr_->chr_; 
        }

        uint32_t get_start() { return ptr_->start_ - 1; }
        uint32_t get_end() { return ptr_->end_ - 1; }

        //std::string& get_gene_type() { return ptr_->anno_type_; }
        //std::string& get_anno_seed() { return ptr_->anno_seed_; }
    
    //private:
        std::unique_ptr<AnnoBedType> ptr_;
    };

    /*
    template<>
    class Adaptor<biovoltron::format::WIG>
    {
    public:
        Adaptor(biovoltron::format::WIG&& wig, size_t idx)
            : ptr_(std::make_unique<biovoltron::format::WIG>(wig))
            , idx_(idx) {}

        std::unique_ptr<biovoltron::format::WIG> 
        get_data_ptr() { return std::move(ptr_); }

        std::string get_chrom() 
        {
            // format: chr1, chr2... -> triming string to 1, 2, 3...
            // eg. chr13 -> 13
            if (ptr_->header.chrom.substr(0, 3) == "chr") 
                return ptr_->header.chrom.substr(3);
            else 
                return ptr_->header.chrom; 
        }  

        // wig sequence starts from 1 (1-start, fully-closed)
        uint32_t get_start() 
        {   
            if (ptr_->header.is_var_step) 
                return ptr_->chrom_pos - 1;
            else
                return ptr_->header.start + (idx_ * ptr_->header.step) - 1;
        }

        uint32_t get_end()
        {
            return get_start() + ptr_->header.span - 1;
        }

    private:
        std::unique_ptr<biovoltron::format::WIG> ptr_;
        size_t idx_ = 0;
    };*/

}
/*
namespace biovoltron::annotator 
{
    
    class VCFFileAdaptor
    {
    public:
        using Adaptor_VCF = Adaptor<biovoltron::format::VCF>;

        template<class IFS>
        VCFFileAdaptor(IFS&& ifs) 
        {
            biovoltron::format::VCF::VCFHeader vcf_h(ifs);
            biovoltron::format::VCF vcf_d(vcf_h);
            
            while (biovoltron::format::VCF::get_obj(ifs, vcf_d)) 
                adaptor_ptrs.emplace_back(Adaptor_VCF(std::move(vcf_d)));
            ifs.close();
        }

        VCFFileAdaptor& operator=(
                const VCFFileAdaptor& vcf_file_adaptor) = delete;

        VCFFileAdaptor& operator=(
                VCFFileAdaptor&& vcf_file_adaptor) = default;

        std::vector<Adaptor_VCF>& operator()()
        {
            return adaptor_ptrs;
        }

    private:
        std::vector<Adaptor_VCF> adaptor_ptrs;
    };

}
*/

namespace biovoltron::annotator 
{
    class GFFFileAdaptor 
    {
    public:
        using Adaptor_GFF = Adaptor<biovoltron::format::GFF>;

        template<class IFS>
        GFFFileAdaptor(IFS &&ifs
                     , biovoltron::format::GFF::GFFFasta &gff_fasta) 
        {
            biovoltron::format::GFF gff_d(gff_fasta);

            while(biovoltron::format::GFF::get_obj( ifs, gff_d ) )
            {  
                adaptor_ptrs.emplace_back(Adaptor_GFF(std::move(gff_d)));
            }
            ifs.close();
        }

        GFFFileAdaptor& operator=(
                const GFFFileAdaptor& gff_file_adaptor) = delete;

        GFFFileAdaptor& operator=(
                GFFFileAdaptor&& gff_file_adaptor) = default;

        std::vector<Adaptor_GFF>& operator()()
        {
            return adaptor_ptrs;
        }

    private:
        std::vector<Adaptor_GFF> adaptor_ptrs;
    };
}

namespace biovoltron::annotator 
{
    class BedSixAdaptor 
    {
    public:
        using Adaptor_BedSix = Adaptor<biovoltron::format::BedSix>;

        template<class IFS>
        BedSixAdaptor(IFS &&ifs) 
        {
            biovoltron::format::BedSix bed_six;
            while(biovoltron::format::BedSix::get_obj( ifs, bed_six ))
            {  
                adaptor_ptrs.emplace_back(Adaptor_BedSix(std::move(bed_six)));
            }
            ifs.close();
        }

        BedSixAdaptor& operator=(
                const BedSixAdaptor& bed_six_file_adaptor) = delete;

        BedSixAdaptor& operator=(
                BedSixAdaptor&& bed_six_file_adaptor) = default;

        std::vector<Adaptor_BedSix>& operator()()
        {
            return adaptor_ptrs;
        }

    private:
        std::vector<Adaptor_BedSix> adaptor_ptrs;
    };
}

/*
namespace biovoltron::annotator
{
    class WIGFileAdaptor
    {
    public:
        using Adaptor_WIG = Adaptor<biovoltron::format::WIG>;

        template<class IFS>
        WIGFileAdaptor(IFS &&ifs, biovoltron::format::WIGHeader& wig_h)
        {
            biovoltron::format::WIG wig(wig_h);

            if (wig_h.is_var_step) 
            {
                while(biovoltron::format::WIG::get_obj( ifs, wig )) {
                    adaptor_ptrs.emplace_back(
                                    Adaptor_WIG(std::move(wig), 0));
                }    
            } 
            else 
            {
                size_t counter = 0;
                while(biovoltron::format::WIG::get_obj( ifs, wig )) {
                    adaptor_ptrs.emplace_back(
                                Adaptor_WIG(std::move(wig), counter++));
                }
            }
            
            ifs.close();
        }

        WIGFileAdaptor& operator=(
                const WIGFileAdaptor& wig_file_adaptor) = delete;

        WIGFileAdaptor& operator=(
                WIGFileAdaptor&& wig_file_adaptor) = default;

        std::vector<Adaptor_WIG>& operator()()
        {
            return adaptor_ptrs;
        }

    private:
        std::vector<Adaptor_WIG> adaptor_ptrs;
    };
}
*/
