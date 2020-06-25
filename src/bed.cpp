/*!
 * @file bed.cpp
 *
 * @brief BED reader source file
 */

#include <iostream>
#include <memory>
#include <sstream>

#include "bed.hpp"

namespace racon {

bool BedFile::Deserialize(const std::string& line, BedRecord& record) {
    if (line.empty()) {
        return false;
    }
    char name_buff[1024];
    int64_t chrom_start = 0, chrom_end = 0;
    int32_t n = sscanf(line.c_str(), "%s %ld %ld", name_buff, &chrom_start, &chrom_end);
    if (n < 3 || chrom_end <= chrom_start) {
        throw std::runtime_error("Invalid BED line: '" + line + "'");
    }
    record = BedRecord(name_buff, chrom_start, chrom_end);
    return true;
}

void BedFile::Serialize(std::ostream& os, const BedRecord& record) {
    os << record.chrom() << " " << record.chrom_start() << " " << record.chrom_end();
}

std::string BedFile::Serialize(const BedRecord& record) {
    std::ostringstream oss;
    Serialize(oss, record);
    return oss.str();
}

BedReader::BedReader(const std::string& in_fn)
    : file_{std::unique_ptr<std::ifstream>(new std::ifstream(in_fn))}
    , in_{*file_.get()}
{
}

BedReader::BedReader(std::istream& in)
    : in_{in}
{
}

bool BedReader::GetNext(BedRecord& record) {
    const bool rv1 = !std::getline(in_, line_).fail();
    if (!rv1)
        return false;
    const bool rv2 = BedFile::Deserialize(line_, record);
    return rv2;
}

std::vector<BedRecord> BedReader::ReadAll(const std::string& fn)
{
    std::ifstream in{fn};
    return ReadAll(in);
}

std::vector<BedRecord> BedReader::ReadAll(std::istream& in)
{
    std::vector<BedRecord> records;
    BedReader reader{in};
    BedRecord record;
    while (reader.GetNext(record)) {
        records.emplace_back(std::move(record));
    }
    return records;
}

}
