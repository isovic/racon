/*!
 * @file bed.hpp
 *
 * @brief BED file containers and parser.
 */

#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace racon {

class BedRecord;

class BedFile {
public:
    static bool Deserialize(const std::string& line, BedRecord& record);
    static void Serialize(std::ostream& os, const BedRecord& record);
    static std::string Serialize(const BedRecord& record);
};

/*
 * \brief BedRecord container.
 * Note: BED records have 0-based coordinates, and the end coordinate is non-inclusive.
*/
class BedRecord {
public:
    ~BedRecord() = default;

    BedRecord() = default;

    BedRecord(std::string _chrom, int64_t _chrom_start, int64_t _chrom_end)
        : chrom_(std::move(_chrom))
        , chrom_start_(_chrom_start)
        , chrom_end_(_chrom_end) {}

    const std::string& chrom() const {
        return chrom_;
    }
    int64_t chrom_start() const {
        return chrom_start_;
    }
    int64_t chrom_end() const {
        return chrom_end_;
    }

    void chrom(const std::string& val) {
        chrom_ = val;
    }
    void chrom_start(int64_t val) {
        chrom_start_ = val;
    }
    void chrom_end(int64_t val) {
        chrom_end_ = val;
    }

    bool operator==(const BedRecord& rhs) const
    {
        return chrom_ == rhs.chrom_ && chrom_start_ == rhs.chrom_start_
                    && chrom_end_ == rhs.chrom_end_;
    }

    std::ostream& operator<<(std::ostream& os) const {
        BedFile::Serialize(os, *this);
        return os;
    }

private:
    std::string chrom_;
    int64_t chrom_start_ = 0;
    int64_t chrom_end_ = 0;
};

class BedReader {
public:
    BedReader(const std::string& in_fn);
    BedReader(std::istream& in);

    static std::vector<BedRecord> ReadAll(const std::string& fn);
    static std::vector<BedRecord> ReadAll(std::istream& in);
    bool GetNext(BedRecord& record);

private:
    std::unique_ptr<std::ifstream> file_;
    std::istream& in_;
    std::string line_;
};

}
