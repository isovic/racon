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

#define MAX_COV 100

struct seqaln_sort_key {
//    inline bool operator() (const SingleSequence** op1, const SingleSequence** op2) {
//      return ((*op1)->get_aln().pos < (*op2)->get_aln().pos);
//    }
  inline bool operator() (const SingleSequence* op1, const SingleSequence* op2) {
    return ((op1)->get_aln().pos < (op2)->get_aln().pos);
  }
};

int AlignmentsToContigs(const SequenceFile &alns, std::vector<std::string> &ctg_names, std::map<std::string, std::vector<const SingleSequence *> > &ctg_alns);
int ExtractAltContigs(std::vector<const SingleSequence *> &ctg_alns, int64_t raw_ctg_len, double coverage_threshold, double percent_overlap, double qv_threshold, std::vector<std::vector<const SingleSequence *> *> &ret_alt_contigs, std::vector<const SingleSequence *> &rejected_alns);
int ConstructContigFromAlns(const SingleSequence &orig_contig, const std::vector<const SingleSequence *> *seq_alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens, SingleSequence &new_contig);
int Consensus(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &alns);



#endif /* SRC_CONSENSUS_CONSENSUS_H_ */
