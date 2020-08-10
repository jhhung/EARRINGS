/// @file header.hpp
/// @brief Header of wig object 

#pragma once
#include <sstream>
#include <vector>
#include <Biovoltron/format/wig/exception.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/constants.hpp>

namespace biovoltron::format::wig {

/// @class Header
/// @brief A header class that store header infomation for 
/// wig object
///
/// Member variable<br>
///         string definition<br>
///         string wig_type<br>
///         string chrom<br>
///         uint32_t start<br>
///         uint32_t step<br>
///         uint32_t span<br>
///         bool is_var_step<br>
///
/// Member functions<br>
///         Constructors<br>
///         operator=<br>
///         to_string<br>
///         set_header<br>
///         preparse<br>
///
/// Non-member functions<br>
///         operator>>
///         operator<<
///     
/// This class stores the header of a wig object, also provide
/// pre-parsing, access, modify functions.
class Header
{
  public:
    /// @brief Default constructor, default is variableStep
    Header()
    : definition(""),
      wig_type("variableStep"),
      chrom(""),
      start(0),
      step(0),
      span(1),
      is_var_step(true)
    {}

    /// @brief Istream constructor
    ///
    /// @param Istream object where the input come from
    Header(std::istream& is)
    {
        preparse(is);
    }

    /// @brief Copy constructor
    ///
    /// @param wig header to be used as source to initialize the 
    /// wig header with
    Header(const Header& other) = default;

    /// @brief Move constructor
    ///
    /// @param wig header to be used as source to initialize the 
    /// wig header with
    Header(Header&& other) = default;

    /// @brief Copy assignment operator
    ///
    /// @param wig header to be used as source to initialize the 
    /// wig header with
    /// @return *this
    Header& operator=(const Header& other) = default;

    /// @brief Move assignment operator
    ///
    /// @param wig header to be used as source to initialize the 
    /// wig header with
    /// @return *this
    Header& operator=(Header&& other) = default;

    /// @brief Print header data as string
    std::string to_string() const
    {
        // numeric to string
        std::string str_start(std::to_string(start));
        std::string str_step(std::to_string(step));
        std::string str_span(std::to_string(span));

        if(is_var_step)
            return definition + "\n" +
                   wig_type + " " +
                   "chrom=" + chrom + " " +
                   "span=" + str_span;
        else
            return definition + "\n" +
                   wig_type + " " +
                   "chrom=" + chrom + " " +
                   "start=" + str_start + " " +
                   "step=" + str_step + " " +
                   "span=" + str_span;
    }

    /// @brief Set header data with string
    ///
    /// @param String to parse
    void set_header(const std::string& str)
    {
        // check header format validability
        size_t def_end;
        if(str.substr(0, 19) == "track type=wiggle_0")
        {
            // get definition line
            def_end = str.find_first_of('\n');
            definition = str.substr(0, def_end);
        }
        else
            throw WIGHeaderException(
                "ERROR: set_header(): parse headers fail\n"
            );

        // store rest of the str
        std::string decl_str = 
            str.substr(def_end+1, str.size()-(def_end+1));

        // get rest of the header member
        std::vector<std::string> v_decl_str;
        boost::algorithm::split(v_decl_str,
                                decl_str,
                                boost::algorithm::is_any_of(" ="),
                                boost::algorithm::token_compress_on);

        // set rest of the header member
        wig_type = v_decl_str[0];
        chrom = v_decl_str[2];
        is_var_step = wig_type == "variableStep" ?true :false;
        // only take odd element in vector
        if(is_var_step)
        {
            if(v_decl_str.size() > 3)
                span = std::stol(v_decl_str[4]);
        }
        else
        {
            start = std::stol(v_decl_str[4]);
            step = std::stol(v_decl_str[6]);

            if(v_decl_str.size() > 7)
                span = std::stol(v_decl_str[8]);
        }
    }

    /// @brief Preparse istream to header data
    ///
    /// @param Istream object to parse
    void preparse(std::istream& is)
    {
        // set definition line
        std::getline(is, definition);

        // check header format validability
        if(definition.substr(0, 19) != "track type=wiggle_0")
            throw WIGHeaderException(
                "ERROR: preparse(): parse headers fail\n"
            );


        // get rest of the header member
        std::string str_decl;
        std::vector<std::string> v_str;
        std::getline(is, str_decl);
        boost::algorithm::split(v_str,
                                str_decl,
                                boost::algorithm::is_any_of(" ="),
                                boost::algorithm::token_compress_on);

        // set rest of the header member
        wig_type = v_str[0];
        chrom = v_str[2];
        is_var_step = wig_type == "variableStep" ?true :false;
        // only take even element in vector
        if(is_var_step)
        {
            if(v_str.size() > 3)
                span = std::stol(v_str[4]);
        }
        else
        {
            start = std::stol(v_str[4]);
            step = std::stol(v_str[6]);
            if(v_str.size() > 7)
                span = std::stol(v_str[8]);
        }
    }

    /// @brief Stream extraction operator which extract header
    /// from istream
    ///
    /// @param An istream where header come from
    /// @param A header object which store the header info
    friend std::istream& operator>>(std::istream& is, Header& header)
    {
        header.preparse(is);
        return is;
    }

    /// @brief Stream insertion operator which inserts header
    /// from ostream
    ///
    /// @param An ostream where header info will be dumped 
    /// @param A header object which store the header info
    friend std::ostream& operator<<(std::ostream& os, const Header& header)
    {
        os << header.to_string();
        return os;
    }

    /// @brief Store wig format definition line
    std::string definition;

    /// @brief Store wig format type: variableStep or fixStep
    std::string wig_type;

    /// @brief Store chromosome type
    std::string chrom;

    /// @brief Store start position (fixStep only)
    uint32_t start;

    /// @brief Store steps between position (fixStep only)
    uint32_t step;

    /// @brief Store span window size (optional)
    uint32_t span;

    /// @brief A flag to tell varStep or fixStep
    bool is_var_step;
};

}
