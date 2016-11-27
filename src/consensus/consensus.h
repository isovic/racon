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
#include <string>
#include <vector>
#include <map>
#include "parameters.h"
#include "mhap.h"

#include "intervaltree/IntervalTree.h"
typedef Interval<const SingleSequence *> IntervalSS;
typedef IntervalTree<const SingleSequence *> IntervalTreeSS;

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

void CreateConsensus(const ProgramParameters &parameters, int32_t num_window_threads, const SingleSequence *contig, const std::vector<SingleSequence *> &ctg_alns, std::map<const SingleSequence *, int64_t> &aln_lens_on_ref, std::string &ret_consensus, FILE *fp_out_cons);

#endif /* SRC_CONSENSUS_CONSENSUS_H_ */
