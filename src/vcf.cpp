/*!
 * @file vcf.cpp
 *
 * @brief VCF data structures and tools.
 */

#include <sstream>
#include "vcf.hpp"

namespace racon {

std::vector<VcfDiff> ExtractVCFEventsFromCigarString(const racon::Cigar& cigar, const std::string& qseq, const std::string& tseq) {

    std::vector<VcfDiff> vcf_diffs;

    const int64_t qlen = qseq.size();
    const int64_t tlen = tseq.size();

    int64_t qpos = 0, tpos = 0;
    int64_t first_diff = -1;
    int64_t first_qpos = -1, first_tpos = -1;

    for (int64_t i = 0; i < static_cast<int64_t>(cigar.size()); ++i) {
        const auto& c = cigar[i];

        if (c.op != '=' && first_diff < 0) {
            first_diff = i;
            first_qpos = qpos;
            first_tpos = tpos;
        }
        // Check if we found a stretch of diffs.
        if (c.op == '=' && first_diff >= 0) {
            const int64_t qspan = qpos - first_qpos;
            const int64_t tspan = tpos - first_tpos;

            VcfDiff v;
            // VCF format is poorly specified. It has a special case when an indel event happens at the first base.
            if (first_tpos == 0) {
                if ((first_tpos + tspan + 1) > tlen || (first_qpos + qspan + 1) > qlen) {
                    std::ostringstream oss;
                    oss << "The entire contig alignment is a diff?"
                            << " first_tpos = " << first_tpos << ", tspan = "
                            << tspan << ", tlen = " << tlen
                            << ", first_qpos = " << first_qpos << ", qspan = "
                            << qspan << ", qlen = " << qlen;
                    throw std::runtime_error(oss.str());
                }
                v.pos = first_tpos;
                v.ref = tseq.substr(first_tpos, tspan + 1);
                v.alt = qseq.substr(first_qpos, qspan + 1);
            } else {
                v.pos = first_tpos - 1;
                v.ref = tseq.substr(first_tpos - 1, tspan + 1);
                v.alt = qseq.substr(first_qpos - 1, qspan + 1);
            }
            vcf_diffs.emplace_back(v);

            first_diff = -1;
            first_qpos = -1;
            first_tpos = -1;
        }
        if (c.op == '=' || c.op == 'X') {
            qpos += c.count;
            tpos += c.count;
        } else if (c.op == 'I') {
            qpos += c.count;
        } else if (c.op == 'D') {
            tpos += c.count;
        }
    }

    return vcf_diffs;
}

}
