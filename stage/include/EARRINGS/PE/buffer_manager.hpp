#pragma once
#include <bitset>
#include <mutex>
#include <vector>
#include <string>

using namespace EARRINGS;
namespace EARRINGS
{
class BufferManager
{
private:
    std::mutex _chunk_mutex;
    std::bitset<32> _chunk_used;
    uint32_t _num_chunks;

public:
    typedef std::unique_lock<std::mutex> Buf_lock;
    std::vector<std::vector<std::string>> buf;
    BufferManager() {}

    void set_chunk_size(uint32_t chunk_size, uint32_t num_chunks)
    {
        _num_chunks = num_chunks;
        _chunk_used.set(); // all set to 1(used)
        _chunk_used <<= _num_chunks;
        buf = std::vector<std::vector<std::string>>{
                          2
                        , std::vector<std::string>(chunk_size * num_chunks)};
    }

    uint32_t inc_chunk_cnt()
    {
        uint32_t buf_idx = -1;
        size_t i(0);
        Buf_lock lock(_chunk_mutex, std::defer_lock);

        while (true)
        {
            if (lock.try_lock())
            {
                if (_chunk_used.all()) // buffer full
                {
                    break;
                }

                for (; i < _num_chunks; ++i)
                {
                    if (_chunk_used[i] == 0)
                    {
                        _chunk_used[i] = 1;
                        buf_idx = i;
                        break;
                    }
                    else
                        ;
                }
                lock.unlock();
                break;
            }
        }

        return buf_idx;
    }

    void dec_chunk_cnt(uint32_t idx)
    {
        Buf_lock lock(_chunk_mutex, std::defer_lock);
        while (true)
        {
            if (lock.try_lock())
            {
                _chunk_used[idx] = 0;
                lock.unlock();
                break;
            }
        }
    }
};
}