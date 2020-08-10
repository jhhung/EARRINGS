/// @file wig.hpp
/// @brief A parser of wig format file 

#pragma once
#include <Biovoltron/format/wig/header.hpp>

namespace biovoltron::format{

using WIGHeader = wig::Header;

/// @class WIG
/// @brief Wig class which store wig data and wig header.
///
/// Member variable<br>
///         WIGHeader& header<br>
///         uint32_t chrom_pos<br>   
///         double data_value<br>
///
/// Member functions<br>
///         Constructors<br>
///         operator=<br>
///         data_to_string<br>
///         get_wig_hader<br>
///         set_wig_data<br>
///         dump<br>
///         get_obj<br>
///
/// Non-member functions<br>
///         operator>>
///         operator<<
///     
/// This class store wig data and header from an wig format file.
/// The wig class only stores one data entry from the wig file. 
///  reference.
class WIG
{
  public:
    /// @brief Copy header contructor
    ///
    /// @param Wig header to be used as source to initialize the
    /// wig object
    WIG(WIGHeader& other_header)
    : header(other_header),
      chrom_pos(0),
      data_value(0.0)
    {}

    /// @brief Copy constructor
    ///
    /// @param Wig object to be used as source to initialize the 
    /// wig object with
    WIG(const WIG& other) = default;

	/// @brief Move constructor
    ///
    /// @param Wig object to be used as source to initialize the 
    /// wig object with
    WIG(WIG&& other) = default;

    /// @brief Copy assignment operator
    ///
    /// @param Wig object to be used as source to initialize the 
    /// wig object with
    /// @return *this
    WIG& operator=(const WIG& other)
    {
        header = other.header;
        chrom_pos = other.chrom_pos;
        data_value = other.data_value;

        return *this;
    }

    /// @brief Move assignment operator
    ///
    /// @param Wig object to be used as source to initialize the 
    /// wig object with
    /// @return *this
    WIG& operator=(WIG&& other)
    {
        header = std::move(other.header);
        chrom_pos = std::move(other.chrom_pos);
        data_value = std::move(other.data_value);

        return *this;
    }

    /// @brief Get the header reference own by the wig object
    /// 
    /// @return The wig header reference
    WIGHeader get_wig_header() const
    {
        return header;
    }

    /// @brief Print wig data to string
    ///
    /// @return A string that contains the data value information
    std::string data_to_string() const
    {
        // numeric to string
        std::string ouput;
        std::string str_chrom_pos = std::to_string(chrom_pos);
        std::string str_data_value = std::to_string(data_value);

        // delete trailing zeros
        size_t point_pos = str_data_value.find_first_of('.');
        size_t last_not_zero = str_data_value.find_last_not_of('0');
        if(str_data_value[point_pos+1] == '0')
            str_data_value.erase(point_pos, std::string::npos);
        else
            str_data_value.erase(last_not_zero+1, std::string::npos);

        // output
        if(header.is_var_step)
            return str_chrom_pos + " " + str_data_value;
        else
            return str_data_value;
    }

    /// @brief Set wig data from string
    ///
    /// @param A reference string that stores data entry
    ///
    /// The expected string is in one of the below two format: <br>
    /// 1. varStep: "   [uint32_t] [double]"
    /// 2. fixStep: "   [double]"
    void set_wig_data(const std::string& str)
    {
        // get first data boudary
        size_t begin = str.find_first_not_of(" \t");
        size_t end = str.find_first_of(" \t");

        if(header.is_var_step)
        {
            chrom_pos = std::stol(str.substr(begin, end-begin+1));
            data_value = std::stod(str.substr(end+1));
        }
        else
            data_value = std::stod(str.substr(begin));
    }

    /// @brief Dump vector of WIG objects (with the same header
    /// reference) to ostream (file), i.e. print a series of wig
    /// object as a wig file.
    /// 
    /// @param An ostream object to dump to, usually file.
    /// @param An vector that stroes a list of wig objects with the
    /// same header reference.
    /// 
    /// The function will first dump header info. After a '\n', dump 
    /// data stored in each wig object and separate them by '\n'.
    static void dump(std::ostream& os, const std::vector<WIG>& v_wig)
    {
        // make sure vector not empty
        if(v_wig.size() == 0)
            return;

        os << v_wig[0].header.to_string() << "\n";

        for(auto item: v_wig)
            os << item.data_to_string() << "\n";
    }

    /// @brief Parse wig data from specific istream
    ///
    /// @param An istream object contains wig data
    /// @param A wig object to store the input wig data
    /// @return The passed in istream object
    ///
    /// Store one data entry from istream.
    /// Upon calling, the function assume that the header has already
    /// been constructed.
    static std::istream& get_obj(std::istream& is, WIG& wig)
    {
        std::string str;
        if(std::getline(is, str))
            wig.set_wig_data(str);

        return is;
    }

    /// @brief Stream extraction operator which extract data
    /// value from istream
    ///
    /// @param An istream where data value come from
    /// @param A wig object which store the data value
    ///
    /// Upon calling, the function assume that the header has already
    /// been constructed.
    friend std::istream& operator>>(std::istream& is, WIG& wig)
    {
        get_obj(is, wig);
        return is;
    }

    /// @brief Stream insertion operator which dump object 
    /// data value to ostream
    ///
    /// @param An ostream object where data value dump to 
    /// @param A wig object which store the data value to dump
    friend std::ostream& operator<<(std::ostream& os, const WIG& wig)
    {
        os << wig.data_to_string();
        return os;
    }

    /// @brief wig header reference, shared by WIG objects 
    /// of the same file
    WIGHeader& header;

    /// @brief chromosome position, only variableStep have it.
    /// If the wig object is fixStep, this member is 0, and should
    /// not be accessed.
    uint32_t chrom_pos; 

    /// @brief data value of specific position
    double data_value;
};

}
