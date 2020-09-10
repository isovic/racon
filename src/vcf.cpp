/*!
 * @file vcf.cpp
 *
 * @brief VCF data structures and tools.
 */

#include <stdexcept>
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
    int64_t n_ops = cigar.size();

    for (int64_t i = 0; i <= n_ops; ++i) {

        if (i < n_ops && cigar[i].op != '=' && first_diff < 0) {
            first_diff = i;
            first_qpos = qpos;
            first_tpos = tpos;
        }

        // Check if we found a stretch of diffs.
        if ((i == n_ops || cigar[i].op == '=') && first_diff >= 0) {
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
        if (i == n_ops) {
            break;
        }
        if (cigar[i].op == '=' || cigar[i].op == 'X') {
            qpos += cigar[i].count;
            tpos += cigar[i].count;
        } else if (cigar[i].op == 'I') {
            qpos += cigar[i].count;
        } else if (cigar[i].op == 'D') {
            tpos += cigar[i].count;
        }

        if (qpos > qlen || tpos > tlen) {
            std::ostringstream oss;
            oss << "Invalid CIGAR vector, the operations do not match the provided sequence "
                << "lengths. qpos = " << qpos << ", tpos = " << tpos << ", qlen = " << qlen << ", tlen = " << tlen;
            throw std::runtime_error(oss.str());
        }
    }

    return vcf_diffs;
}

}
