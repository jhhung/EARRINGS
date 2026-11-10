#pragma once

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <map>
#include <mutex>
#include <optional>
#include <utility>

namespace EARRINGS
{
// Hands items that were finished in any order to one consumer in the order of their index
// (0, 1, 2, ...). push() never waits, so a worker is never held up by a slow consumer; the
// producer limits how far ahead of the consumer it runs by looking at consumed().
template<typename T>
class OrderedQueue
{
private:
    std::mutex _mutex;
    std::condition_variable _cv;
    std::map<uint32_t, T> _items;
    uint32_t _next = 0;
    std::optional<uint32_t> _total;  // number of items there will ever be, once known
    bool _closed = false;
    std::atomic<uint32_t> _consumed{0};

public:
    // Add item `idx`; each index is pushed once. Dropped after close().
    void push(uint32_t idx, T item)
    {
        {
            std::lock_guard<std::mutex> lock(_mutex);
            if (_closed) return;
            _items.emplace(idx, std::move(item));
        }
        _cv.notify_one();
    }

    // No index >= total will be pushed.
    void set_total(uint32_t total)
    {
        {
            std::lock_guard<std::mutex> lock(_mutex);
            _total = total;
        }
        _cv.notify_one();
    }

    // The consumer is gone (or stops caring): drop what is queued and what is pushed later.
    void close()
    {
        {
            std::lock_guard<std::mutex> lock(_mutex);
            _closed = true;
            _items.clear();
        }
        _cv.notify_all();
    }

    // Wait for the next item in order. Returns false once all `total` items have been taken,
    // or after close().
    bool pop_next(T &out)
    {
        std::unique_lock<std::mutex> lock(_mutex);
        _cv.wait(lock, [this] {
            return _closed || _items.count(_next) || (_total && _next >= *_total);
        });
        if (_closed) return false;
        auto node = _items.extract(_next);
        if (node.empty()) return false;  // all taken
        out = std::move(node.mapped());
        ++_next;
        _consumed.store(_next, std::memory_order_release);
        return true;
    }

    // Number of items taken so far.
    uint32_t consumed() const { return _consumed.load(std::memory_order_acquire); }
};
}
