#pragma once
#include <cstdint>

using namespace EARRINGS;
namespace EARRINGS 
{
class RWCount
{
private:
    uint32_t _rcount;

public:
    RWCount() : _rcount(0) {}

    uint32_t rend_idx = -1, rend_count = -1;

    uint32_t rcount_fetch_add()
    {
        _rcount += 1;
        return _rcount - 1;
    }

    // Chunks read so far
    uint32_t rcount() const { return _rcount; }
};
}