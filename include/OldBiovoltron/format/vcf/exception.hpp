#pragma once
#include <stdexcept>
namespace biovoltron::format::vcf{
/**
 * @class VCFException
 * @brief An exception class which inherit std::runtime_error
 *
 * This exception will be thrown if actions which belongs to VCF 
 * have exceptions.
 */
class VCFException : public std::runtime_error
{
  public:

    /**
     * @brief Message assign constructor
     *
     * @param m Message string wanted to be assigned
     *
     * This message will be shown if exception thrown
     */
    VCFException(const std::string& m)
    : runtime_error(std::string("VCFException " + m))
    {
    }
};

/**
 * @class VCFHeaderException
 * @brief An exception class which inherit std::runtime_error
 *
 * This exception will be thrown if actions which belongs to 
 * VCF::VCFHeader have exceptions.
 */
class VCFHeaderException : public std::runtime_error
{
  public:
    
    /**
     * @brief Message assign constructor
     *
     * @param m Message string wanted to be assigned
     *
     * This message will be shown if exception thrown
     */
    VCFHeaderException(const std::string& m)
    : runtime_error(std::string("VCFHeaderException " + m))
    {
    }
};
}
