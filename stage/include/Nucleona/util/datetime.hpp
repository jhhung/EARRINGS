#pragma once
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
namespace nucleona::util {

/**
 * %Y: 4 digits year
 * %m: [01, 12] month
 * %d: [01, 31] day ( 2 digits prefer )
 * %e: [1, 31] day ( 1 digits prefer )
 * %H: [00, 23] hour
 * %I: [01, 12] hour
 * %M: [00, 59] min
 * %S: [00, 60] sec
 * 
 */
struct DateTimeFn {
    template<class CLOCK, class DU>
    std::string operator()(
        const std::string& format, 
        const std::chrono::time_point<CLOCK, DU>& time
    ) const {
        auto in_time = std::chrono::system_clock::to_time_t( time );
        std::stringstream ss;
        ss << std::put_time( std::localtime(&in_time), format.c_str() );
        return ss.str();
    }
};
constexpr DateTimeFn datetime;
}