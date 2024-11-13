#pragma once
#include <iostream>
#include <cstdint>
#include <vector>
#include <string>
#include <OldBiovoltron/utils/string_view.hpp>
#include <OldBiovoltron/format/vcf/exception.hpp>
#include <OldBiovoltron/format/vcf/utils.hpp>
#include <functional>

namespace biovoltron::format {

class VCF;

}

namespace biovoltron::format::vcf{

/**
 * @class Header
 * @brief Headers of a vcf file, nested in VCF
 *
 * Member variable<br>
 *          HeaderFmtType headers_<br>
 *
 * This class stores the headers and order of data fields
 * Also provide header pre-parsing, access, modify functions
 */
class Header
{
friend ::biovoltron::format::VCF;
  public:
    ///Exactly vector<string>, just for shorten code
    using VecStrType = std::vector<std::string>;

    /**
     * @brief Indicates a line of header
     * 
     * HeaderType vector<pair>, indicates contents of a header
     *
     * The pair is (string, VecStrType), indicates sub-field of 
     * content
     *
     * The first of pair is string, indicates name of sub-field
     *
     * The second of pair is VecStrType, indicates contents of 
     * sub-field
     */
    using HeaderType =
        std::vector<
        std::pair<
            std::string, VecStrType
        >>;
    
    /** 
     * @brief Indicates hole headers of a vcf file
     *
     * HeaderFmtType is vector<pair> with size 8, first for 
     * user-defined, others are default
     *
     * The pair is (string, vector<HeaderType>), indicate a 
     * category of headers
     *
     * The first of pair is string, indicates name of category
     *
     * The second of pair is vector of HeaderType, indicate 
     * headers belongs to this category
     */    
    using HeaderFmtType =
        std::vector<    
        std::pair<
            std::string, std::vector<HeaderType>
        >>;

    /**
     * @brief The function type of each field parser.
     * 
     * FieldParser is a function type which expected to do 
     * the string to integer/float/string vector conversion.
     * The target data type is defined by VCF member variable.
     */
    using FieldParser = 
	std::function<void(const utils::StringView&, VCF&)>;

    /**
     * @brief The function type of each field writer.
     * 
     * FieldWriter is a function type which expected to do 
     * the integer/float/string vector to conversion.
     * The source data type is defined by VCF member variable.
     */
    using FieldWriter = 
	std::function<bool(std::string&, const VCF&)>;
    
  private:
    ///A constant size for size of string reserve
    static const size_t buf_default_size_ = 8192;
    
    ///A constant size for size of field(vector) reserve
    static const size_t field_order_default_size_ = 16;
    
    ///A private member variable which stores headers 
    HeaderFmtType headers_;

    ///The parser function for each field, the order is same as the field_order_.
    std::vector<FieldParser> field_parser_;

    ///The writer function for each field, the order is same as the field_order_.
    std::vector<FieldWriter> field_writer_;

    ///A public member variable which stores order of data fields
    std::vector<std::string> field_order_;

  public:
    /**
     * @brief A default constructor initialize Field names and 
     * reserve space of variable
     *
     * Time Complexity: O(1)
     */
    Header() noexcept
    {
        headers_.resize(8);
        headers_[0].first = "";
        headers_[1].first = "INFO";
        headers_[2].first = "FILTER";
        headers_[3].first = "FORMAT";
        headers_[4].first = "ALT";
        headers_[5].first = "contig";
        headers_[6].first = "SAMPLE";
        headers_[7].first = "PEDIGREE";

        field_order_.reserve(field_order_default_size_);
    }

    /**
     * @brief A constructor which will preparse headers from given istream
     *
     * @param is A istream which contains headers can be preparse
     *
     * Do the same things like default constructor, and use 
     * Header::preparse() to preparse the header by given is
     *
     * Time Complexity: O(i*j)<br>
     *	    \e i: the number of headers<br>
     *	    \e j: the average length of headers<br>
     *
     * @sa Header::preparse()
     */
    Header(std::istream& is)
	: Header()
    {
        preparse(is);
    }

    /**
     * @brief istream operator >> overload
     *
     * @param is An istream which contains headers can preparse
     * @param vcf_h A Header to store headers preparsed
     * @return Is identical to parameter \e is
     *
     * Use Header::preparse() to parse vcf formatted headers in \e is
     *
     * @sa Header::preparse()
     */
    friend std::istream& operator>>(std::istream& is, 
                    Header& vcf_h)
    {
        vcf_h.preparse(is);
        return is;
    }

    /**
     * @brief ostream operator << overload
     *
     * @param os An ostream where want to dump to
     * @param vcf_h A Header want to dump
     * @return Is identical to parameter \e os
     *
     * Use Header::to_string to generate string of vcf formatted 
     * headers, and output it to \e os
     *
     * Time Complexity: O(i*j)<br>
     *	    \e i: the number of headers<br>
     *	    \e j: the average length of headers
     *
     * @sa Header::to_string()
     */
    friend std::ostream& operator<<(std::ostream& os, 
                    Header& vcf_h)
    {
        os << vcf_h.to_string();
        return os;
    }

    /**
     * @brief Add a header from given string
     *
     * @param buf A string which stores headers
     * @exception Thrown when \e buf can't be parse to vcf formatted 
     * header
     *
     * First split the string by Header::parse_one_line_header(), 
     * and get a vector<string> which contains string of each FIELD, 
     * then pass contents to Header::parse_one_line_header() by loop. 
     * And decide to throw exception or not by flag returned
     *
     * Time Complexity: O(i)<br>
     *	    \e i: the length of the string
     *
     * @sa VCF::split_with_quote_handle(), 
     * Header::parse_one_line_header()
     */
    void add_header(std::string& buf)
    {
        std::vector<utils::StringView> order_split;
        utils::StringView sv_line(buf), delim("\n");
        int flag;

        split_with_quote_handle(order_split, sv_line, delim);
        if (order_split.back() == "")
        {
            order_split.pop_back();
        }

        for (auto& i: order_split)
        {
            flag = parse_one_line_header(i);

            if (flag < 0)
            {
                throw VCFHeaderException(
                    "ERROR: add_header(): add given header fail\n"
                );
            }
        }
    }

    /**
     * @brief Convert Header::headers_ to vcf formatted string, and 
     * return it
     *
     * @return A string of headers which is vcf formatted
     *
     * Use Header::get_one_line_header() to get a string of a vcf 
     * formatted header, and append all headers of string to a 
     * string \e buf by loop, then return \e buf
     *
     * Time Complexity: O(i*j)<br>
     *	    \e i: the number of headers<br>
     *	    \e j: the average length of headers
     *
     * @sa Header::get_one_line_header()
     */
    std::string to_string() const noexcept
    {
        std::string buf = "";

        buf.reserve(buf_default_size_);
        for (const auto& category_pair: headers_)
        {
            for (const auto& one_header: category_pair.second)
            {
                if (one_header.size() != 1)
                {
                    buf.append("##");
                    buf.append(category_pair.first);
                    buf.append("=");
                }
                get_one_line_header(buf, one_header);
            }
        }

        buf.append("#");
        for (auto& field: field_order_)
        {
            buf.append(field);
            buf.append("\t");
        }

        if (buf.size() == 1)
        {
            buf = "";
        }
        else
        {
            buf.pop_back();
        }

        return buf;
    }

    /**
     * @brief To preparse headers from given istream
     *
     * @param is A istream contains headers can be preparse
     * @exception Thrown VCFHeaderException when fail to parse 
     * headers
     *
     * Use std::getline() to get a line of header, then use 
     * Header::parse_one_line_header() to extract and store data in 
     * Header::headers_ repeatly.<br>
     * If the return of Header::parse_one_line_header() is 0, means 
     * header preparse is over; If it is 1, means should keep 
     * parsing; And if it's less than 0, means parsing fail, so will 
     * throw an exception
     *
     * Time Complexity: O(i*j)<br>
     *	    \e i: the number of headers<br>
     *	    \e j: the average length of headers
     *
     * \sa Header::parse_one_line_header, VCFHeaderException
     */
    void preparse(std::istream& is)
    {
        std::string line;
        int flag;

        while (std::getline(is, line))
        {
            utils::StringView sv_line(line);
            flag = parse_one_line_header(sv_line);
            if (flag == 0)
            {
                break;
            }
            else if (flag == 1)
            {
                continue;
            }
            else
            {
                throw VCFHeaderException(
                    "ERROR: preparse(): parse headers fail\n"
                );
            }
        }
    }

  private:
    /**
     * @brief Convert a Header::HeaderType to string which meets 
     * vcf format
     *
     * @param buf A string to store converted header
     * @param one_header A HeaderType object which want to 
     * convert to string
     *
     * Convert a line of header to string by concatenate data in 
     * \e one_header
     *
     * Time Complexity: O(n)<br>
     *	    \e n: the length of this header
     */
    void get_one_line_header(std::string& buf, const HeaderType& one_header) const noexcept
    {
        if (one_header.size() == 1)
        {
            buf.append("##");
            buf.append(one_header[0].first);
            buf.append("=");
            buf.append(one_header[0].second[0]);
        }
        else
        {
            buf.append("<");
            for (const auto& key_pair: one_header)
            {
                buf.append(key_pair.first);
                buf.append("=");
                for (const auto& content: key_pair.second)
                {
                    buf.append(content);
                    buf.append(";");
                }
                buf.back() = ',';
            }
            buf.back() = '>';
        }

        buf.append("\n");
    }

    /**
     * @brief Parse a line of header from given string_view
     *
     * @param line A string_view which want to parse to header
     * @return A int {-1,0,1} to indicate fail,header,field 
     * order respectively
     *
     * Parse a line of header. Since header contains two categories, 
     * pure header and Field order, we distinguish it by their first 
     * two char. Pure header is \e ##, and Field order is \e #*(* is 
     * wild card).<br>
     * If it's pure header then we divide it into name and data 
     * parts, and use Header::store_header() to store it. And if it's 
     * Field order, we use Header::field_parser_ and 
     * Header::field_writer_ to record the order<br>
     * If this header is pure header, return 1; If it's Field order, 
     * return 0. If it has wrong format, return -1
     *
     * Time Complexity: O(n)<br>
     *	    \e n: the length of header
     *
     * @sa Header::field_parser_, Header::field_writer_, 
     * Header::store_header()
     */
    int parse_one_line_header(const utils::StringView& line);
   

    /**
     * @brief store header with given header name and content
     *
     * @param name Header name of header
     * @param content Header content of header
     *
     * First find the Field index specified by <em>name</em>, then 
     * split content by comma with Header::split_with_quote_handle(). 
     * And for each segment, split them by colon, then put them in 
     * the correct place of Header::headers_<br>
     * Or if we can't find Field specific by <em>name</em>, means 
     * it's not a default Field, so we put them in index 0, which is 
     * a space prepared for non-default Field
     *
     * Time Complexity: O(n)<br>
     *	    \e n: the length of header
     *
     * @sa Header::headers_, Header::split_with_quote_handle()
     */
    void store_header(const utils::StringView& name, 
            utils::StringView& content)
    {
        size_t field_idx(1);

        for (; field_idx < headers_.size(); field_idx++)
        {
            if (headers_[field_idx].first == name)
            {
                break;
            }
        }

        HeaderType h_tmp; 

        if (field_idx != headers_.size())
        {
            std::vector<utils::StringView> comma_split, colon_split;
            utils::StringView delim(",");

            remove_angle_brackets(content);
            split_with_quote_handle(comma_split, content, delim);

            delim = ";";
            for (size_t i(0); i < comma_split.size(); i++)
            {
                size_t pos(comma_split[i].find("="));
                split_with_quote_handle(
                    colon_split,
                    comma_split[i].substr(pos + 1),
                    delim
                );

                VecStrType vs_tmp;
                for (auto& i : colon_split)
                {
                    vs_tmp.emplace_back(i);
                }

                h_tmp.emplace_back(
                    std::make_pair(
                        comma_split[i].substr(0, pos), vs_tmp
                    )
                );
            }

            headers_[field_idx].second.emplace_back(h_tmp);
        }
        else
        {
            h_tmp.emplace_back(
                std::make_pair(
                    name, VecStrType(1, std::string(content))
                )
            );
            headers_[0].second.emplace_back(h_tmp);
        }
    }

    /**
     * @brief Remove angle brackets for given string_view
     *
     * @param line A string_view may or may not contains angle 
     * bracket
     *
     * Since every content of pure header must be surrounded by angle 
     * brackets, so we have to remove them before parsing contents of 
     * headers. In this method we use std::find_first_of() and 
     * std::find_last_of to find angle brackets and remove them
     *
     * Time Complexity: O(n)<br>
     *	    \e n: the length of header
     */
    inline void remove_angle_brackets(utils::StringView& line) noexcept
    {
        size_t pos1, pos2;

        if ((pos1 = line.find_first_of('<')) != std::string::npos)
        {
            line.remove_prefix(pos1 + 1);
        }
        if ((pos2 = line.find_last_of('>')) != std::string::npos)
        {
            line.remove_suffix(line.size() - pos2);
        }
    }
};
/*Header end*/
}
