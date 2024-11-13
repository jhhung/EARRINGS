#pragma once

#include <iostream>
#include <stdexcept>
#include <string>

namespace biovoltron::annotator
{
class Overlap_Annotator: public std::runtime_error
{
public:
    Overlap_Annotator(const std::string& msg)
    : runtime_error(std::string(msg)) {}
};
}