/**
 * @file vcf.hpp
 * @brief a parser of vcf format file
 *
 * @author JHH corp
 */

#pragma once
#include <cstring>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/mpl/string.hpp>
#include <OldBiovoltron/utils/string_view.hpp>
#include <OldBiovoltron/format/vcf/header.hpp>
#include <OldBiovoltron/format/vcf/utils.hpp>


namespace biovoltron::format{


/**
 * @class VCF
 * @brief A class which store vcf entry and header, class VCFHeader
 * is nested inside
 *
 * Member variable:<br>
 *      VCFHeader headers_<br>
 *      bool is_no_data<br>
 *      std::string chrom_, id_, ref_, alt_, filter_<br>
 *      int32_t pos_<br>
 *      float qual_<br>
 *      std::vector<InfoContentType> info_<br>
 *      std::vector<SampleContentType> fmt_data_<br>
 *
 * This class store VCF entry and headers of this entry, class 
 * VCFHeaders are also define nested, and also provides entry parsing
 *  , access, modify and entry access functions
 */
class VCF
{
  public:

    ///VCFHeader is the VCF header type
    using VCFHeader = vcf::Header;

    ///VecStrType is just for shorten code, actually a vector<string>
    using VecStrType = typename VCFHeader::VecStrType;

    /**
     * @brief Indicates a line of header
     * 
     * HeaderType is vector<pair>, indicates contents of a header
     *
     * The pair is (string, VecStrType), indicates sub-field of 
     * content
     *
     * The first of pair is string, indicates name of sub-field<br>
     * The second of pair is VecStrType, indicates contents of 
     * sub-field
     */
    using HeaderType = typename VCFHeader::HeaderType;

    /** 
     * @brief Indicates hole headers of a vcf file
     *
     * HeaderFmtType is vector of pair with size 8, first for 
     * user-defind, others are default<br>
     *
     * The pair is (string, vector<HeaderType>), indicate a 
     * category of headers<br>
     *
     * The first of pair is string, indicates name of category
     *
     * The second of pair is vector of HeaderType, indicate 
     * headers belongs to this category
     */    
    using HeaderFmtType = typename VCFHeader::HeaderFmtType;

    /**
     * @brief Indicates contents of sub-field of INFO field
     *
     * InfoContentType is pair of two string, the first is name of 
     * sub-field, and the second is value of sub-field
     */
    using InfoContentType = std::pair<std::string, std::string>;

    /**
     * @brief Indicates contents of sub-field of FORMAT and samples 
     * fields
     *
     * SampleContentType is pair of string and vector<string>, the 
     * first is name of sub-field of format, and the second is values 
     * of samples of sub-field
     */
    using SampleContentType = 
        std::pair<std::string, std::vector<std::string>>;

  private:
    ///A constant size for size of string reserve
    static const size_t buf_default_size_ = 256;
    
    ///A VCFHeader reference contains headers
    VCFHeader& headers_;

  public:
    ///A flag to represent if has no data
    bool is_no_data;

    ///Some strings which store data of CHROM, ID, REF, ALT, FILTER
    std::string chrom, id, ref, alt, filter;

    ///A 32bit-int to store data of POS
    int32_t pos;

    ///A float to store data of QUAL
    float qual;

    ///A vector of InfoContentType to store data of INFO
    std::vector<InfoContentType> info;
    
    ///A vector of SampleContentType to store data of FORMAT and samples
    std::vector<SampleContentType> fmt_data;
    
    /**
     * @brief A constructor copy given VCF::VCFHeader reference to 
     * member variable VCF::headers_
     *
     * @param vcf_h A VCF::VCFHeader object used to initialize this 
     * VCF by using copy
     *
     * Copy the given reference and initialize some member variable
     *
     * Time Complexity: O(1)
     *
     * @sa vcf::Header
     */
    VCF(VCF::VCFHeader& vcf_h) noexcept
    : headers_      ( vcf_h )
    , is_no_data    (  true )
    , chrom         (   "." )
    , id            (   "." )
    , ref           (   "." )
    , alt           (   "." )
    , filter        (   "." )
    , pos           (   -1  )
    , qual          (   0.0 )
    {
    }

    VCF(const VCF& vcf_d) = default;
    VCF(VCF&& vcf_d) = default;

    /**
     * @brief = copy assign operator which can assign data to this 
     * VCF object from other VCF
     *
     * @param vcf_d Another VCF variable which whose data will be  
     * copy to this VCF variable
     *
     * @return Reference of this object
     *
     * A copy assign = operator overloading which copy data from VCF 
     * another VCF object. Containing headers_, is_no_data, chrom, 
     * id, alt, filter, pos, qual, info, fmt_data
     *
     * Time Complexity: O(1)
     */
    VCF& operator=(const VCF& vcf_d) noexcept
    {
	this->headers_ = vcf_d.headers_;
	this->is_no_data = vcf_d.is_no_data;
	this->chrom = vcf_d.chrom;
	this->id = vcf_d.id;
	this->ref = vcf_d.ref;
	this->alt = vcf_d.alt;
	this->filter = vcf_d.filter;
	this->pos = vcf_d.pos;
	this->qual = vcf_d.qual;
	this->info = vcf_d.info;
	this->fmt_data = vcf_d.fmt_data;

	return *this;
    }

    /**
     * @brief = move assign operator which can move data to this VCF 
     * object from other VCF
     *
     * @param vcf_d Another VCF variable which whose data will be  
     * move to this VCF variable
     *
     * @return Reference of this object
     *
     * A move assign = operator overloading which move data from VCF 
     * another VCF object by using std::move. Containing headers_, 
     * is_no_data, chrom, id, alt, filter, pos, qual, info, fmt_data
     *
     * Time Complexity: O(1)
     */
    VCF& operator=(VCF&& vcf_d) noexcept
    {
	this->headers_ = vcf_d.headers_;
	this->is_no_data = std::move(vcf_d.is_no_data);
	this->chrom = std::move(vcf_d.chrom);
	this->id = std::move(vcf_d.id);
	this->ref = std::move(vcf_d.ref);
	this->alt = std::move(vcf_d.alt);
	this->filter = std::move(vcf_d.filter);
	this->pos = std::move(vcf_d.pos);
	this->qual = std::move(vcf_d.qual);
	this->info = std::move(vcf_d.info);
	this->fmt_data = std::move(vcf_d.fmt_data);

	return *this;
    }

    /**
     * @brief < operator which can help to compare VCF
     *
     * @param vcf_d Another VCF variable which gonna compare VCF::pos 
     * with this VCF variable
     *
     * @return if VCF::pos of this VCF variable is less than VCF::pos 
     * of vcf_d
     *
     * A simple < operator overloading which compare VCF by their 
     * VCF::chrom and VCF::pos. First compare VCF::chrom, if they're 
     * the same then compare VCF::pos
     *
     * Time Complexity: O(1)
     */
    bool operator<(const VCF& vcf_d) const noexcept
    {
	if (chrom < vcf_d.chrom)
	    return true;
	else if (chrom == vcf_d.chrom)
	    if (pos < vcf_d.pos)
		return true;
	    else
		return false;
	else
	    return false;
    }

    /**
     * @brief > operator which can help to compare of VCF
     *
     * @param vcf_d Another VCF variable which gonna compare VCF::pos 
     * with this VCF variable
     *
     * @return if VCF::pos of this VCF variable is great than 
     * VCF::pos of vcf_d
     *
     * A simple > operator overloading which compare VCF by their 
     * CHROM and POS field. First compare CHROM, if they're the same 
     * then compare POS
     *
     * Time Complexity: O(1)
     */
    bool operator>(const VCF& vcf_d) const noexcept
    {
	if (chrom > vcf_d.chrom)
	    return true;
	else if (chrom == vcf_d.chrom)
	    if (pos > vcf_d.pos)
		return true;
	    else
		return false;
	else
	    return false;
    }

    /**
     * @brief >> operator which get entry from an istream, then 
     * store read vcf format data into VCF
     *
     * @param is An istream where to get entry
     * @param vcf A VCF object to store gotten entry
     *
     * @return Is identical to parameter is
     *
     * Just call VCF::get_obj() to help parse given istream, and 
     * store them in given vcf object.<br>
     * Should have the some property to get_obj function
     *
     * @sa VCF::get_obj()
     */
    friend std::istream& operator>>(std::istream& is, 
                    VCF& vcf)
    {
        get_obj(is, vcf);
        return is;
    }

    /**
     * @brief << operator which get entry from an ostream
     *
     * @param os An ostream where to dump entry
     * @param vcf A VCF object which want to dump
     *
     * @return Is identical to parameter os
     *
     * Just call VCF::to_string() to help change data of given vcf 
     * object to string then output to given ostream, and store them 
     * in given vcf object.<br>
     * Should have the some property to VCF::to_string()
     *
     * @sa VCF::to_string()
     */
    friend std::ostream& operator<<(std::ostream& os, 
                    const VCF& vcf)
    {
        os << vcf.to_string();
        return os;
    }

    /**
     * @brief can get VCF::headers_ of this VCF object
     *
     * @return A VCF::VCFHeader reference which is private member 
     * variable VCF::headers_
     */
    const VCF::VCFHeader& get_header() const noexcept
    {
        return headers_;
    }

    /**
     * @brief Can set VCF::data_ of this VCF object by given string
     *
     * @param buf A string used to set VCF::data_ of this VCF object
     *
     * Store given string in an istringstream, then use 
     * VCF::get_obj() to parse this istringstream<br>
     * Should have almost the same property of VCF::get_obj()
     *
     * @sa VCF::get_obj()
     */
    void set_vcfdata(const std::string& buf)
    {
        std::stringstream iss;
        iss.str(buf);
        get_obj(iss, *this);
    }

    /**
     * @brief Can get VCF::data_ of this VCF object
     *
     * @return A string of data in vcf format
     *
     * Concatenate string of all data Field by order defined in 
     * VCF::header_, and make this string fit vcf format, then 
     * return it
     *
     * Time Complexity: O(n)<br>
     *	    \e n is the number of Field
     */
    std::string to_string() const;

    /**
     * @brief Parse a line of vcf formatted input from specified 
     * istream
     *
     * @param is An istream contains entry data
     * @param vcf_d A VCF object used to store parsed entry data
     * @return An istream identical to parameter is
     *
     * Use std::getline from an istream then parse and store in 
     * vcf_d. Parsing a line by the order defined in VCF::headers_ , 
     * and store them into corresponding member variable<br>
     * Since the format of storing is related to headers stored in 
     * vcf_d, so should use vcf::Header::preparse() to parse headers first, 
     * or it will parse headers automatically
     *
     * Time Complexity: O(i+j*k)<br>
     *	    \e i: the length of str<br>
     *	    \e j: the size of OTHER FIELD<br>
     *	    \e k: number of sub-fields of samples<br>
     * Time Complexity: O(O(m*n)+i+j*k), if need to preparse<br>
     *	    \e m: lines of headers<br>
     *	    \e n: average length of each header<br>
     *
     * @sa VCF::split_with_quote_handle()
     */
    static std::istream& get_obj(std::istream& is, VCF& vcf_d);


    /**
     * @brief Convert vector of vcf to string then dump to an ostream
     *
     * @param os An ostream to dump vcf format data
     * @param v_vcf A vector<VCF>, which we want to dump them to vcf 
     * format
     * @param is_print_header A bool as a flag to decide print 
     * headers before data entrys or not, default value is true
     *
     * Convert vector<VCF> to string by first use 
     * <em>os << vcf::Header::to_string</em>  to output formatted 
     * headers, if \e is_print_header is false, skip this step. Then 
     * use VCF::to_string() to output VCF in while loop
     *
     * Time Complexity: O(m*n*(i+j*k) + (p*q))<br>
     *	    \e m: size of v_vcf<br>
     *	    \e n: number of FIELD<br>
     *	    \e i: number of default FIELD<br>
     *	    \e j: size of OTHER FIELD<br>
     *	    \e k: sub-field number of FORMAT FIELD<br>
     *	    \e p: number of headers<br>
     *	    \e q: average length of each header
     * 
     * @sa vcf::Header::to_string(), VCF::to_string()
     */
    static void dump(std::ostream& os, const std::vector<VCF>& v_vcf, 
	    bool is_print_header = true);
};

}
