#pragma once
#include <atomic>

using namespace EARRINGS;
namespace EARRINGS 
{
class RWCount
{
private:
    std::atomic_uint32_t _wcount;
    uint32_t _rcount;

public:
    RWCount() : _wcount(0), _rcount(0) {}

    bool all_finished = false;
    uint32_t rend_idx = -1, rend_count = -1;

    uint32_t rcount_fetch_add()
    {
        _rcount += 1;
        return _rcount - 1;
    }

    uint32_t wcount_fetch_add()
    {
        return _wcount.fetch_add(1);
    }
};
}