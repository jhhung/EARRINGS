#pragma once

#include <vector>
#include <string>
#include <utility>
#include <boost/tokenizer.hpp>

namespace biovoltron::format {

class BEDGRAPH;

}

/**
 * @class Header
 * @brief Header of a bedgraph file, and is nested in BEDGRAPH
 * 
 *
 */
namespace biovoltron::format::bedgraph{

class Header
{
friend ::biovoltron::format::BEDGRAPH;
  public:
    using HeaderValueType = std::pair<std::string, std::string>;
    using HeaderType = std::vector<HeaderValueType>;
    
  private:
    ///A private member variable storing a header
    HeaderType header;
    
  public:
    ///A boolean to determine whether to print out quotation marks
    bool double_quoted = false;

    ///@brief A default constructor
    Header(){};

    /**
     * @brief A constructor contains a header from given istream
     * 
     * @param is: A istream which contains headers
     */
    Header(std::istream& is)
    {
        std::string str;
        std::getline(is, str);
        set_header(str);
    }
    /**
     * @brief copy constructor
     * 
     * @param _header 
     */
    Header(const Header& _header) = default;

    /**
     * @brief move constructor
     * 
     * @param _header 
     */
    Header(Header&& _header) = default;

    /**
     * @brief operator = overload
     * 
     * @param _header 
     * @return Header& 
     */
    Header& operator= (const Header& _header) = default;
    Header& operator= (Header&& _header) = default;

    /**
     * @brief return a pointer pointing to this header
     * 
     * @return const std::shared_ptr<Header> 
     */
    const std::shared_ptr<Header> to_ptr()
    {
        return std::make_shared<Header>(*this);
    }

    /**
     * @brief Get the header object in string format
     * 
     * @return std::string 
     */    
    std::string to_string()
    {
        std::string buf = "";
        for(auto& header_type_pair : header)
        {
            buf.append(" ");
            buf.append(header_type_pair.first);
            buf.append("=");
            if(double_quoted &&
                (header_type_pair.first == "name" || header_type_pair.first == "description"))
            {
                buf.append("\"");
                buf.append(header_type_pair.second);
                buf.append("\"");
            }
            else buf.append(header_type_pair.second);
        }
        return "track"+buf;
    }

    /**
     * @brief Set the header object 
     * 
     * @param str: A string which we want to parse
     */
    void set_header(std::string str)
    {
        if(str.find("\"")) double_quoted = true;
        boost::escaped_list_separator<char> els("", " ", "\"");
        boost::tokenizer<boost::escaped_list_separator<char>> 
            tok(str, els);
        for(auto&t : tok)
        {
            std::size_t found_equal = t.find("=");
            if(found_equal != std::string::npos)
            {
                HeaderValueType value_type (
                    std::move(t.substr(0, found_equal)), 
                    std::move(t.substr(found_equal+1))
                    );
                header.emplace_back(value_type);
            }
        }
    }

    /**
     * @brief istream operator >> overload
     * 
     * @param is: An istream which contains the objects of header
     * @param header: A header to store parsed objects
     * @return std::istream&: Is identical to parameter is
     */
    friend std::istream& operator>>(std::istream& is, Header& header)
    {
        std::string str;
        std::getline(is, str);
        header.set_header(str);
        return is;
    }

    /**
     * @brief ostream operator << overload
     * 
     * @param os: A ostream where to dump to
     * @param header: A header being dump 
     * @return std::ostream&: Is identical to parameter os
     */
    friend std::ostream& operator<<(std::ostream& os, Header& header)
    {
        os << header.to_string();
        return os;
    }
};
}