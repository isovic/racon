/*
 * consensus.h
 *
 *  Created on: Feb 15, 2016
 *      Author: isovic
 */

#ifndef SRC_CONSENSUS_CONSENSUS_H_
#define SRC_CONSENSUS_CONSENSUS_H_

#include "sequences/single_sequence.h"
#include "sequences/sequence_file.h"
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include "parameters.h"
#include "mhap.h"

#include "intervaltree/IntervalTree.h"
typedef Interval<const SingleSequence *> IntervalSS;
typedef IntervalTree<const SingleSequence *> IntervalTreeSS;

// Data for a sampled sequence alignment. Sampling is performed in reference coordinates.
class SampledAlignment {
 public:
  SampledAlignment(const SingleSequence *_seq, OverlapLine &_mhap) : seq(_seq), mhap(_mhap) { };
//  const SingleSequence *seq;
//  std::shared_ptr<SingleSequence> seq;
  const SingleSequence *seq;
  std::unordered_map<int32_t, int32_t> qpos;  // key is the ref pos, and value is the query pos for the ref pos.
  OverlapLine mhap;

  // Find a query pos for a given ref pos, or return an invalid value (-1).
  int32_t find(int32_t rpos) {
    auto it = qpos.find(rpos);
    if (it == qpos.end()) { return -1; }
    return it->second;
  }
};

typedef Interval<std::shared_ptr<SampledAlignment>> IntervalSampled;
typedef IntervalTree<std::shared_ptr<SampledAlignment>> IntervalTreeSampled;

struct seqaln_sort_key {
  inline bool operator() (const SingleSequence* op1, const SingleSequence* op2) {
    return ((op1)->get_aln().get_pos() < (op2)->get_aln().get_pos());
  }
};

struct ContigOverlapLocation {
  int64_t start = 0, end = 0, ctg_id = 0;
};

int GroupAlignmentsToContigs(const SequenceFile &alns, double qv_threshold, std::vector<std::string> &ctg_names, std::map<std::string, std::vector<const SingleSequence *> > &ctg_alns);
int ConsensusFromAln(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &alns);

int GroupOverlapsToContigs(const std::vector<OverlapLine> &sorted_overlaps, std::map<int64_t, ContigOverlapLocation> &map_ctg_to_overlaps);
int ConsensusFromOverlaps(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &reads, const std::map<std::string, int64_t> &qname_to_ids, const std::vector<OverlapLine> &overlaps);

void ExtractWindowFromAlns(const SingleSequence *contig, const std::vector<SingleSequence *> &alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens,
                           IntervalTreeSS &aln_interval_tree, int64_t window_start, int64_t window_end, double qv_threshold, bool use_contig_qvs,
                           std::vector<std::string> &window_seqs, std::vector<std::string> &window_qv, std::vector<const SingleSequence *> &window_refs,
                           std::vector<uint32_t> &window_starts, std::vector<uint32_t> &window_ends,
                           std::vector<uint32_t> &starts_on_read, std::vector<uint32_t> &ends_on_read, FILE *fp_window);

void CreateConsensusAln(const ProgramParameters &parameters, int32_t num_window_threads, const SingleSequence *contig, const std::vector<SingleSequence *> &ctg_alns, std::map<const SingleSequence *, int64_t> &aln_lens_on_ref, std::string &ret_consensus, FILE *fp_out_cons);


void ExtractWindowFromSampledSeqs(const SingleSequence *contig, IntervalTreeSampled &interval_tree,
                                  int64_t window_start, int64_t window_end, double qv_threshold, bool use_contig_qvs,
                                  std::vector<std::string> &window_seqs, std::vector<std::string> &window_qv, std::vector<const SingleSequence *> &window_refs,
                                  std::vector<uint32_t> &window_starts, std::vector<uint32_t> &window_ends,
                                  std::vector<uint32_t> &starts_on_read, std::vector<uint32_t> &ends_on_read, FILE *fp_window);
void PrepareAndSampleOverlaps(const SequenceFile &refs, const SequenceFile &reads,
                             const std::vector<OverlapLine> &sorted_overlaps, int64_t start_overlap, int64_t end_overlap, int32_t num_threads,
                             std::vector<std::shared_ptr<SampledAlignment>> &extracted_overlaps, int32_t window_len, bool verbose_debug);
void CreateConsensusSampling(const ProgramParameters &parameters, int32_t num_window_threads, const SingleSequence *contig,
                             IntervalTreeSampled &interval_tree, std::string &ret_consensus, FILE *fp_out_cons);
void PerformSampling(std::shared_ptr<SampledAlignment> sampling_ovl, const SingleSequence* ref, int64_t window_len);
void PerformDummySampling(std::shared_ptr<SampledAlignment> sampling_ovl, const SingleSequence* ref, int64_t window_len);
int ConsensusFromOverlapsWSampling(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &reads,
                                const std::map<std::string, int64_t> &qname_to_ids, const std::vector<OverlapLine> &sorted_overlaps);
char *Reverse(const char *s, int32_t len);

#endif /* SRC_CONSENSUS_CONSENSUS_H_ */
