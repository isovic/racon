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
    return ((op1)->get_aln().get_pos() < (op2)->get_aln().get_pos());
  }
};

int GroupAlignmentsToContigs(const SequenceFile &alns, double qv_threshold, std::vector<std::string> &ctg_names, std::map<std::string, std::vector<const SingleSequence *> > &ctg_alns);
//int ExtractAltContigs(std::vector<const SingleSequence *> &ctg_alns, int64_t raw_ctg_len, double coverage_threshold, double percent_overlap, double qv_threshold, std::vector<std::vector<const SingleSequence *> *> &ret_alt_contigs, std::vector<const SingleSequence *> &rejected_alns);
//int ConstructContigFromAlns(const SingleSequence &orig_contig, const std::vector<const SingleSequence *> *seq_alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens, SingleSequence &new_contig);
//int Consensus(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &alns);
//int RunMSAFromSystem(const ProgramParameters &parameters, std::string &cons);

int MajorityVoteFromMSALocal(std::string pir_path, std::string *cons);
int RunMSAFromSystemLocal(const ProgramParameters &parameters, std::string window_path, std::string &cons);
//void ExtractWindowFromAlns(const std::vector<const SingleSequence *> &alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens, int64_t window_start, int64_t window_end, std::vector<std::string> window_seqs, FILE *fp_window);
void ExtractWindowFromAlns(const SingleSequence *contig, const std::vector<const SingleSequence *> &alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens, int64_t window_start, int64_t window_end, std::vector<std::string> &window_seqs, std::vector<std::string> &window_qv, FILE *fp_window);
int ConsensusDirectFromAln(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &alns);
void CreateConsensus(const ProgramParameters &parameters, const SingleSequence *contig, std::vector<const SingleSequence *> &ctg_alns, std::map<const SingleSequence *, int64_t> &aln_lens_on_ref, std::string &ret_consensus, FILE *fp_out_cons);

#endif /* SRC_CONSENSUS_CONSENSUS_H_ */
