#pragma once
#include <iostream>
#include <string>

namespace biovoltron::annotator{
    class Annotation {
    public:
        Annotation(std::string chrom, uint32_t start, uint32_t end):
                chrom_(chrom), start_(start), end_(end) {}
        std::string get_chrom() const { return chrom_; }
        uint32_t get_start() const { return start_; }
        uint32_t get_end() const { return end_; }
    private:
      std::string chrom_;
      uint32_t start_;
      uint32_t end_;
    };

};
