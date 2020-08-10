/// @file exception.hpp
/// @brief Exception wig related class 
#pragma once
#include <stdexcept>

namespace biovoltron::format::wig{

/// @class WIGHeaderException
/// @brief An exception class which inherit std::runtime_error
///
/// This exception will be thrown if action which belongs to 
/// WIGHeader have exceptions.
class WIGHeaderException : public std::runtime_error
{
    /**
     * @brief Message assign constructor
     *
     * @param m Message string wanted to be assigned
     */
  public:
    WIGHeaderException(const std::string& mesg)
    : runtime_error(std::string("WIGHeaderException " + mesg)) {}
};

}
