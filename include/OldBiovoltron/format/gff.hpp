/**
 * @file gff.hpp
 * @brief a parser of gff format file
 *
 * @author JHH corp
 */
#pragma once

#include <cstring>
#include <exception>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <memory>
#include <OldBiovoltron/utils/string_view.hpp>
#include <OldBiovoltron/base_vector.hpp>
#include <OldBiovoltron/format/fasta.hpp>
#include <boost/algorithm/string.hpp>

namespace biovoltron::format{

/**
 * @class GFF
 * @brief A class which store gff entry 
 *
 * Member variable:<br>
 *      GFFFasta fasta<br>
 *      MemType member<br>
 *
 *	This class stores GFF's columns as well as fasta data (if exist), 
 *  and it also provides entry parsing and accessing functions.
 *
 *  GFF contains the following 9 columns:
 *  seqid, source, type, start, end, score, phase and attributes.
 *  Please visit http://gmod.org/wiki/GFF3 for more information.
 */

class GFF
{
  public:
    
    ///UInt is 64 bits unsigned int
    using UInt = uint64_t;

    ///MemType stores gff's columns using vector of strings
    using MemType = std::vector<std::string>;

    ///GFFFasta uses std::shared_ptr<std::string> to store fasta's content
    using GFFFasta = std::shared_ptr<std::string>;

    ///GFFFasta is a type that stores reference to shared_ptr of fasta data
    GFFFasta fasta;

    ///MemType is a type used stores gff's columns.
    MemType member;

  public:
    /**
     * @brief A getter function that returns the 
     * corresponding gff member. "pos" should not
     * be greater than 8, since there're only 9 columns.
     *
     * @return A string which stores the member.
     */
    template<size_t pos>
    std::string get() const
    {
        return member[pos];
    }

    /**
     * @brief One can use a GFF object to call this 
     * function to get a string which contains all the 
     * columns of the GFF object.
     *
     * @return A string which stores the GFF content.
     */
    std::string to_string() const
    {
        std::string data_str("");
        for (auto& str: member)
        {
            data_str += str;
            
            if (str != member[member.size() - 1])
            {
                data_str += '\t';
            }
            else
            ;  // the last string, don't concatenate '\t'
        }

        return data_str;
    }

    /**
     * @brief If a GFF file contains fasta data, then one
     * can use this function to get all the fasta sequence 
     * srting.
     *
     * @return A string which stores fasta sequence in GFF file.
     */
    std::string fasta_to_string() const
    {   
        return *fasta;
    }

    /**
     * @brief A constrctor with fasta sequence pointer.
     * See GFFFasta to learn more about input type.
     *
     * @param fasta A GFFFasta object that is used 
     * to store the fasta sequence pointer.
     */
//default constructor with table
    GFF(const GFFFasta& fasta)
    : fasta( fasta )
    , member(  )
    {
        member.reserve(9);
    }

    /**
     * @brief A copy constructor which copies another 
     * GFF object.
     *
     * Assume A, B are two different GFF object. Then 
     * A(B) implies copying B's content to A's.
     * 
     * @param gff A GFF object which one intends to copy.
     */
//copy constructor
    GFF( const GFF& gff ) = default;

    /**
     * @brief A move constructor which moves another 
     * GFF object.
     *
     * Assume A, B are two different GFF object. Then 
     * A(std::move(B)) implies moving B's content to A's.
     * 
     * @param gff A GFF object which one intends to move.
     */
//move constructor
    GFF( GFF&& gff )
    {
        fasta = gff.fasta;
        member = std::move(gff.member);
    }

    /**
     * @brief An copy assignment copies from another GFF object.
     *
     * Assume A, B are two different GFF object. Then A = B implies
     * copying B's content to A's.
     *
     * @param gff A GFF which one intends to copy.
     */
//copy assignment
    GFF& operator = ( const GFF& gff )
    {
        fasta = gff.fasta;
        member = gff.member;

        return *this;
    }
      
    /**
     * @brief A assignment move another GFF object.
     *
     * Assume A, B are two different GFF object. Then 
     * A = std::move(B) implies moving B's content to A's.
     *
     * @param gff A GFF object which one intends to move.
     */
//move assignment
    GFF& operator = ( GFF&& gff )
    {
        fasta = gff.fasta;
        member = std::move(gff.member);

        return *this;
    }

    /**
     * @brief A >> operator that copies entries from an istream.
     *
     * Assume A is an istream object and B is a GFF object.
     * Doing A >> B means entries are copied from A to B.
     *
     * @param is An istream which one intends to get entries.
     * @param gff A GFF object used to store the input entries.
     *
     * @return Is This is identical to the input parameter "is".
     */
//istream >> operator
friend std::istream& operator >> (std::istream& is, GFF& gff)
{
    get_obj(is, gff);
    return is;
}

    /**
     * @brief A << operator which gets all the GFF entries 
     * from an ostream.
     *
     * @param os An ostream used dump all the entries.
     * @param gff A GFF object which one intends to dump.
     *
     * @return Is This is identical to the input paramter "os".
     */
//ostream << operator
friend std::ostream& operator << (std::ostream& os, GFF& gff)
{
    os << gff.to_string();
    return os;
}

    /**
     * @brief A operator< which compares two GFF object. The 
     * first element to compare is chromosome then the start 
     * position and lastly, the end position.
     *
     * @param gff_0 A GFF object.
     * @param gff_1 A GFF object.
     *
     * @return A bool that indicates whether gff_0 is small 
     * than gff_1.
     */
//operator <
friend bool operator < ( const GFF& gff_0, const GFF& gff_1 )
{
    if ( gff_0.get<0>() > gff_1.get<0>() )
    {
	    return false;
    }
    else if ( std::stoi(gff_0.get<3>()) > std::stoi(gff_1.get<3>()) )
    {
	    return false;
    }
    else if ( 
            std::stoi(gff_0.get<3>()) == std::stoi(gff_1.get<3>())
            && std::stoi(gff_0.get<4>()) < std::stoi(gff_1.get<4>())
            )
    {
        return false;
    }
    else
    {
        return true;
    }
}

    /**
     * @brief A function used to parse the input stream and stores
     * the parsed content into a GFF object.
     *
     * @param is An istream contains entry data.
     * @param vcf_d A GFF object used to store parsed entry data.
     * 
     * @return An istream indicates the current stream pointers.
     */
//get object
    static std::istream& get_obj( std::istream& is, GFF& gff)
    {
        std::string str;
        while (std::getline(is, str)) 
        {
            if(str[0] != '#' && str[0] != '>')
            {
                // read chromosome data
                boost::split( gff.member, str, boost::is_any_of("\t"));
            
                return is;
            } 
            else if (str == "##FASTA")
            {
                // Fasta sequences are preceded by a ##FASTA line.
                std::string fstr("");

                // read to the end of the data
                while(std::getline( is, str ))
                {
                    fstr += str + '\n';
                }
                *(gff.fasta) = fstr;
                is.seekg( 0, std::ios::end );

                return is;
            }
            else
            ;  // continue reading data 
            
        } 
        
        return is;
    }

    /**
     * @brief Convert the content of vector of all GFFs to string 
     * then dump the string to an ostream.
     *
     * @param os An ostream used to dump GFF format data.
     * @param v_gff A vector of GFF objects which one intends to dump.
     */

/// dumpFunc
    static void dump(std::ostream& os, std::vector<GFF>& v_gff)
    {
        std::string str;
        for ( auto& i: v_gff )
        {
            str += i.to_string() + '\n'; 
        }
        
        str.pop_back();  // remove the last '\n'
        os << str;
    }
};

}

