/**
 * @file bedgraph.hpp
 * @author Nathalie
 * @date 2018-09-28
 */

#pragma once

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <Biovoltron/format/bedgraph/header.hpp>
namespace biovoltron::format{


/**
 * @class BEDGRAPH
 * @brief A class which store bedgraph entry and header, 
 *        Class Header is nested inside
 * 
 */

class BEDGRAPH
{
  public:
    using BGHeader = bedgraph::Header;
    using BGHeaderValueType = typename BGHeader::HeaderValueType;
    using BGHeaderType = typename BGHeader::HeaderType;
    
    
  public:
    ///A shared pointer pointing to a header
    std::shared_ptr<BGHeader> ptr_BGheader;
    ///A string to store chromosome name
    std::string name;
    ///two unsigned int to store start and end position of chromosome respectively
    uint32_t pos_start, pos_end;
    ///A float to store score
    float score;
    ///A unsigned int to store float precision
    uint32_t float_precision = 0;

    ///A default constructor
    BEDGRAPH(){};

    /**
     * @brief Construct a new BEDGRAPH object 
     *        with current header pointer 
     * 
     * @param ptr_header: a pointer pointing to a header
     * @param precision: an unsigned integer to set the precision
     */
    BEDGRAPH(const std::shared_ptr<BGHeader>& ptr_header, 
            const uint32_t precision)
        : ptr_BGheader(ptr_header), float_precision(precision){};

    /**
     * @brief default copy constructor
     * 
     * @param bdgraph: the bedgraph that the user want to copy
     */
    BEDGRAPH(const BEDGRAPH& bdgraph) = default;
    BEDGRAPH(BEDGRAPH&& bdgraph) = default;


    /**
     * @brief operator = overload
     * 
     * @param bedgraph: the bedgraph that the user want to be assigned
     * @return BEDGRAPH&
     */
    BEDGRAPH& operator= (const BEDGRAPH& bedgraph) = default;
    BEDGRAPH& operator= (BEDGRAPH&& bedgraph) = default;

    /**
     * @brief Set all data objects
     * 
     * @param data: BedGraph data. 
     *              4 columns: {chrom chromStart chromEnd dataValue}
     */
    void set_bgdata(const std::string& data)
    {
        std::vector<std::string> vec_str_data;
        boost::split(vec_str_data, data, boost::is_any_of(" \t"), boost::token_compress_on);
    
        name = vec_str_data[0];
        pos_start = std::stoi(vec_str_data[1]);
        pos_end   = std::stoi(vec_str_data[2]);
        score = std::stof(vec_str_data[3]);
    }

    /**
     * @brief Get the bedgraph data object in string format
     * 
     * @return std::string 
     */
    std::string to_string()
    {
        std::stringstream str_score;
        str_score << std::fixed 
                  << std::setprecision(float_precision) 
                  << score;

        std::string data(name);
        data.append(" ");
        data.append(std::to_string(pos_start));
        data.append(" ");
        data.append(std::to_string(pos_end));
        data.append(" ");
        data.append(str_score.str());
        
        return data;
    }

    /**
     * @brief Get bedgraph header object in string format
     * 
     * @return std::string 
     */
    std::string get_header()
    {
        return (*ptr_BGheader).to_string();
    }


    /**
     * @brief >> operator which gets entry from an istream
     * 
     * @param is: An istream where to get entry
     * @param bedgraph: A bedgraph object to store entry
     * @return std::istream&: Is identical to parameter is 
     */
    friend std::istream& operator>>
            (std::istream& is, BEDGRAPH& bedgraph)
    {
        std::string str;
        std::getline(is, str);
        bedgraph.set_bgdata(str);
        return is;
    }

    /**
     * @brief << operator which gets entry from an ostream
     * 
     * @param os: An ostream where to dump entry
     * @param bedgraph: A bedgraph object which want to dump
     * @return std::ostream&: Is identical to parameter os
     */
    friend std::ostream& operator<<
            (std::ostream& os, BEDGRAPH& bedgraph)
    {
        os << bedgraph.to_string();
        return os;
    }

    /**
     * @brief Parse bedgraph entry from specific istream
     * 
     * @param is: An istream contains entry data
     * @param ptr_header: A pointer pointing to the header
     * @param precision: An unsigned integer to set the precision
     * @return BEDGRAPH: The bedgraph gained from the istream
     */
    
    static BEDGRAPH get_obj(std::istream& is, 
            const std::shared_ptr<BGHeader>& ptr_header, 
            const uint32_t precision)
    {
        
        std::string str;
        getline(is, str);

        ///Construct the first bedgraph and build its shared_ptr for the header
        BEDGRAPH bedgraph_(ptr_header, precision);
        bedgraph_.set_bgdata(str);

        return bedgraph_;
    }
    
    /**
     * @brief A template function 
     *        to dump a BEDGRAPH object from a container to an ostream
     * 
     * @param os: An ostream to dump bedgraph format data
     * @param container_bedgraph: A container of bedgraph
     */
    template<typename T>
    static void dump
            (std::ostream& os, T& container_bedgraph)
    {
        os << container_bedgraph[0].get_header() << std::endl;
        for(auto& _bedgraph : container_bedgraph)
        {
            os << _bedgraph << std::endl;
        }
    } 
};
}