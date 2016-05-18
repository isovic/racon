/*
 * consensusdirect.cc
 *
 *  Created on: Feb 29, 2016
 *      Author: isovic
 */

#include "consensus/consensus.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include <stdint.h>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <omp.h>
#include "spoa.hpp"
#include "graph.hpp"
#include "intervaltree/IntervalTree.h"
#include "libs/edlib.h"
#include "pileup.h"

typedef Interval<const SingleSequence *> IntervalSS;
typedef IntervalTree<const SingleSequence *> IntervalTreeSS;

void prepare_indices1(std::vector<uint32_t>& dst, const std::vector<std::string>& sequences, bool sorted) {
    dst.resize(sequences.size());
    std::iota(dst.begin(), dst.end(), static_cast<uint32_t>(0));

    if (sorted) {
        std::sort(dst.begin(), dst.end(),
            [&](uint32_t lhs, uint32_t rhs) {
                return sequences[lhs].size() > sequences[rhs].size();
            }
        );
    }
}

// #define WINDOW_OUTPUT_IN_FASTQ

int GroupAlignmentsToContigs(const SequenceFile &alns, double qv_threshold, std::vector<std::string> &ctg_names, std::map<std::string, std::vector<const SingleSequence *> > &ctg_alns) {
  ctg_names.clear();
  ctg_alns.clear();

  for (int64_t i=0; i<alns.get_sequences().size(); i++) {
    if (alns.get_sequences()[i]->get_aln().IsMapped() == false) continue;

    if (qv_threshold > 0.0) {
      double average_bq = alns.get_sequences()[i]->CalcAverageBQ();
      if (average_bq >= 0 && average_bq < qv_threshold) { continue; }
    }

    auto it = ctg_alns.find(alns.get_sequences()[i]->get_aln().get_rname());
    if (it != ctg_alns.end()) {
      it->second.push_back((const SingleSequence *) (alns.get_sequences()[i]));
    } else {
      ctg_alns[alns.get_sequences()[i]->get_aln().get_rname()] = std::vector<const SingleSequence *> {(const SingleSequence *) alns.get_sequences()[i]};
      ctg_names.push_back(alns.get_sequences()[i]->get_aln().get_rname());
    }
  }

  return 0;
}

void ExtractWindowFromAlns(const SingleSequence *contig, const std::vector<const SingleSequence *> &alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens,
                           IntervalTreeSS &aln_interval_tree, int64_t window_start, int64_t window_end, double qv_threshold,
                           std::vector<std::string> &window_seqs, std::vector<std::string> &window_qv, std::vector<const SingleSequence *> &window_refs,
                           std::vector<uint32_t> &window_starts, std::vector<uint32_t> &window_ends, FILE *fp_window) {
  if (window_start > window_end) {
    return;
  }

  int64_t temp_window_end = std::min((int64_t) window_end, (int64_t) (contig->get_sequence_length()-1));
  window_refs.push_back(contig);
  window_seqs.push_back(GetSubstring((char *) (contig->get_data() + window_start), (temp_window_end - window_start + 1)));
  std::string dummy_quals((temp_window_end - window_start + 1), '!');
  window_qv.push_back(dummy_quals);
  window_starts.push_back(0);
  window_ends.push_back(temp_window_end - window_start);

  // Find seqs which fall into the window region.
  std::vector<IntervalSS> intervals;
  aln_interval_tree.findOverlapping(window_start, temp_window_end, intervals);

  // For each seq, extract its segment which falls into the window.
  for (int64_t i=0; i<intervals.size(); i++) {
    auto seq = intervals[i].value;
    auto aln = seq->get_aln();

    int64_t start_cig_id = 0, end_cig_id = 0;
    int64_t start_seq = aln.FindBasePositionOnRead(window_start, &start_cig_id);
    int64_t end_seq = aln.FindBasePositionOnRead(temp_window_end, &end_cig_id);
    uint32_t seq_start_in_window = 0;
    uint32_t seq_end_in_window = temp_window_end - window_start;

    if (start_seq == -1) {
      start_seq = aln.GetClippedBasesFront();

      seq_start_in_window = aln.get_pos() - 1 - window_start;
      seq_start_in_window = std::max((uint32_t) 0, (uint32_t) ((int32_t) seq_start_in_window - 0));

    } else if (start_seq < 0) {
      fprintf (stderr, "ERROR: start_seq is < 0 and != -1! start_seq = %ld\n", start_seq); exit(1);
    }

    if (end_seq == -2) {
      end_seq = seq->get_data_length() - 1 - aln.GetClippedBasesBack();
      seq_end_in_window = (aln.get_pos() - 1 + aln.GetReferenceLengthFromCigar()) - window_start;
      seq_end_in_window = std::min((uint32_t) (temp_window_end - window_start), (uint32_t) ((int32_t) seq_end_in_window + 0));
    } else if (end_seq < 0) {
      fprintf (stderr, "ERROR: end_seq is < 0 and != -2!\n"); exit(1);
    }

    std::string seq_data = GetSubstring((char *) (seq->get_data() + start_seq), end_seq - start_seq + 1);
    std::string seq_qual = (seq->get_quality() != NULL) ? (GetSubstring((char *) (seq->get_quality() + start_seq), end_seq - start_seq + 1)) : (std::string((end_seq - start_seq + 1), '!' + 0));

    double avg_qual;
    for (int64_t j=0; j<seq_qual.size(); j++) {
      avg_qual += (double) (seq_qual[j] - '!');
    }
    avg_qual /= std::max((double) seq_qual.size(), 1.0);

    if (avg_qual >= qv_threshold) {
      window_refs.push_back(seq);
      window_seqs.push_back(seq_data);
      window_starts.push_back(seq_start_in_window);
      window_ends.push_back(seq_end_in_window);
      window_qv.push_back(seq_qual);
    }

//    if (avg_qual < qv_threshold) {
//      seq_qual = std::string(seq_data.size(), '!');
//    }
//    window_seqs.push_back(seq_data);
//    window_starts.push_back(seq_start_in_window);
//    window_ends.push_back(seq_end_in_window);
//    window_qv.push_back(seq_qual);


    if (fp_window) {
      #ifndef WINDOW_OUTPUT_IN_FASTQ
        fprintf (fp_window, ">%s Window_%d_to_%d\n%s\n", seq->get_header(), window_start, temp_window_end, window_seqs.back().c_str());
      #else
        fprintf (fp_window, "@%s Window_%d_to_%d\n%s\n", seq->get_header(), window_start, temp_window_end, window_seqs.back().c_str());
        fprintf (fp_window, "+\n%s\n", window_qv.back().c_str());
      #endif
    }
  }

}

int ConsensusDirectFromAln(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &alns) {
  LOG_MEDHIGH("Running consensus - directly from alignments.\n");

  std::vector<std::string> ctg_names;
  std::map<std::string, std::vector<const SingleSequence *> > all_ctg_alns;

  // Separate alignments into groups for each contig.
  // Alignments which are called unmapped will be skipped in this step.
  // Also, alignments are filtered by the base quality if available.
  LOG_MEDHIGH("Separating alignments to individual contigs.\n");
  GroupAlignmentsToContigs(alns, -1.0, ctg_names, all_ctg_alns);

  // Verbose.
  LOG_MEDHIGH("In total, there are %ld contigs, each containing:\n", ctg_names.size());
  for (int32_t i=0; i<ctg_names.size(); i++) {
    LOG_MEDHIGH("\t[%ld] %s %ld alignments\n", i, ctg_names[i].c_str(), all_ctg_alns.find(ctg_names[i])->second.size());
  }

  // Hash the sequences by their name.
  std::map<std::string, const SingleSequence *> rname_to_seq;
  for (int32_t i=0; i<contigs.get_sequences().size(); i++) {
    rname_to_seq[contigs.get_sequences()[i]->get_header()] = contigs.get_sequences()[i];
    rname_to_seq[TrimToFirstSpace(contigs.get_sequences()[i]->get_header())] = contigs.get_sequences()[i];
  }

  // Hash all the alignment lengths (which will be used a lot).
  std::map<const SingleSequence *, int64_t> aln_lens_on_ref;
  for (int64_t i=0; i<alns.get_sequences().size(); i++) {
    aln_lens_on_ref[alns.get_sequences()[i]] = alns.get_sequences()[i]->get_aln().GetReferenceLengthFromCigar();
  }

  // Clear the output file for consensus.
  FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "w");
  fclose(fp_out_cons);

  // For each contig (draft can contain multiple contigs), process alignments which only map to that particular contig.
  for (int32_t i=0; i<ctg_names.size(); i++) {
    const SingleSequence *contig = rname_to_seq[ctg_names[i]];
    auto it = all_ctg_alns.find(ctg_names[i]);
    if (it == all_ctg_alns.end()) {
      FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Something strange happened. Contig name, which was extracted from alignments, cannot be found in the std::map containing those same alignments.");
      // Exits.
    }
    // Get alignments for current contig.
    std::vector<const SingleSequence *> &ctg_alns = it->second;

    // This sorts ascending by the pos field.
    std::sort(ctg_alns.begin(), ctg_alns.end(), seqaln_sort_key());

    FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "a");
    std::string consensus;
    if (parameters.do_pileup == false) {
      CreateConsensus(parameters, contig, ctg_alns, aln_lens_on_ref, consensus, fp_out_cons);

    } else {
      Pileup pileup(contig, ctg_alns);
//      pileup.Verbose(stdout);
      pileup.GenerateConsensus(5, consensus);
      fprintf (fp_out_cons, ">Consensus_0\n%s\n", consensus.c_str());
    }
    fclose(fp_out_cons);

    ///////////////////////////////////////
    LOG_MEDHIGH_NOHEADER("\n");
    LOG_ALL("Processed %ld bp of %ld bp (100.00%%)\n", contig->get_data_length(), contig->get_data_length());
    LOG_MEDHIGH_NOHEADER("\n");
  }

  return 0;
}

int ConvertMSAToAln(const std::string &ref_msa, const std::string &seq_msa, int64_t window_start, SequenceAlignment &aln) {
  if (ref_msa.size() != seq_msa.size() || ref_msa.size() == 0 || seq_msa.size() == 0) {
    return 1;
  }

  // The CIGAR vector and the position should not be cleared as they will be filled up accross windows.
  int64_t aln_pos = aln.get_pos();  // aln_pos is 1 based. If it's equal to 0, it's not a valid position.

  int64_t pos_on_ref = 0;
  int64_t pos_on_seq = 0;
  for (int64_t i=0; i<ref_msa.size(); i++) {
    if (seq_msa[i] != '-') {
      if (aln_pos == 0) { aln_pos = window_start + pos_on_ref + 1; }
    }

    char op = 0;
    if (ref_msa[i] == '-' && seq_msa[i] == '-') {
      continue;
    } else if (ref_msa[i] == seq_msa[i]) {
      op = 'M';
    } else if (ref_msa[i] != seq_msa[i] && ref_msa[i] != '-' && seq_msa[i] != '-') {
      op = 'X';
    } else if (ref_msa[i] == '-' && seq_msa[i] != '-') {
      op = 'I';
    } else if (ref_msa[i] != '-' && seq_msa[i] == '-') {
      if (aln_pos != 0) { // In aln_pos == 0 the sequence hasn't already started, which means that the 'D' op should not be counted.
        op = 'D';
      }
    }

    if (aln.get_cigar().size() > 0 && aln.get_cigar().back().op == op) {
      aln.cigar().back().count += 1;
    } else if (op != 0) {
      CigarOp new_cigar_op;
      new_cigar_op.op = op;
      new_cigar_op.count = 1;
      new_cigar_op.pos_ref = -1; // pos_on_ref;     // I can't guarantee the correctness of these positions here, because it depends on the previous windows on the ref.
      new_cigar_op.pos_query = -1; // pos_on_seq;
      aln.cigar().push_back(new_cigar_op);
    }

    if (ref_msa[i] != '-') { pos_on_ref += 1; }
    if (seq_msa[i] != '-') { pos_on_seq += 1; }
  }

  aln.set_pos(aln_pos);

  return 0;
}

int FilterOverhangsFromMsa(const std::vector<std::string> &msa, std::string &consensus) {
  if (msa.size() == 0) { return 1; }

  int64_t seq_len = msa[0].size();
  int32_t first_msa_seq = 1;  // The 0-th sequence is the layout (with qualities set to '!').
  int32_t last_msa_seq = msa.size() - 1;  // The last msa sequence is the consensus sequence, also aligned to the MSA.
  int32_t max_cov = last_msa_seq - first_msa_seq;
  int32_t start_pos = 0, end_pos = seq_len - 1;

  for (start_pos = 0; start_pos < last_msa_seq; start_pos++) {
    int32_t dash_count = 0;
    for (int32_t j=first_msa_seq; j<msa.size(); j++) { if (msa[j][start_pos] == '-') { dash_count += 1; } }
//    if ((max_cov - dash_count) >= overhang_min_threshold) { break; }
    if (dash_count < max_cov/2) { break; }
  }

  for (end_pos = (seq_len - 1); end_pos >= 0; end_pos--) {
    int32_t dash_count = 0;
    for (int32_t j=first_msa_seq; j<last_msa_seq; j++) { if (msa[j][end_pos] == '-') { dash_count += 1; } }
//    if ((max_cov - dash_count) >= overhang_min_threshold) { break; }
    if (dash_count < max_cov/2) { break; }
  }

//  fprintf (stderr, "start_pos = %d, end_pos = %d, seq_len = %d\n", start_pos, end_pos, seq_len);

  std::string new_cons = "";
  std::stringstream ss;
  auto& old_cons = msa.back();
  for (int32_t i=start_pos; i<=end_pos; i++) {
//    if (old_cons[i] != '-') { ss << old_cons[i]; }

//    if (old_cons[i] != '-') { std::string base = MajorityVotePos(msa, first_msa_seq, last_msa_seq-1, i); ss << base; new_cons += (base == "") ? "-" : base; }
//    else {new_cons += " "; }
    if (old_cons[i] != '-') { ss << old_cons[i]; }
//    else {
//      std::string base = MajorityVotePos(msa, first_msa_seq, last_msa_seq-1, i);
//      ss << base;
//      new_cons += (base == "") ? "-" : base;
//    }
  }

//  printf ("%s\n", new_cons.c_str());

//  for (int32_t i=0; i<seq_len; i++) {
//    if (old_cons[i] == '-') { printf ("-"); continue; }
//    int32_t dash_count = 0;
//    for (int32_t j=first_msa_seq; j<last_msa_seq; j++) { if (msa[j][i] == '-') { dash_count += 1; } }
//    if (dash_count < max_cov/2) {
//      ss << old_cons[i];
//      printf ("%c", old_cons[i]);
//    } else {
//      printf (".");
//    }
//  }
//  printf ("\n");

  consensus = ss.str();

  return 0;
}

std::string MajorityVotePos(const std::vector<std::string> &msa, int32_t first_seq, int32_t last_seq, int32_t pos) {
  if (msa.size() == 0) { return std::string(""); }

  int64_t seq_len = msa[0].size();
  std::stringstream ss;

  // Count occurrences for the column.
  int32_t base_counts[256] = {0};
  for (int32_t j=first_seq; j<=last_seq; j++) {
    base_counts[toupper(msa[j][pos])] += 1;
  }

  std::vector<int32_t> filtered_base_counts = {base_counts['A'], base_counts['C'], base_counts['T'], base_counts['G'], base_counts['-']};
  std::vector<size_t> indices;
  ordered_sort_vector(filtered_base_counts, indices);
  if (indices.back() == 4) { return std::string(""); }
  char ret = "ACTG-"[indices.back()];
  return std::string(1, ret);

//  int64_t sum_base_counts = base_counts['A'] + base_counts['C'] + base_counts['T'] + base_counts['G'];
//  int64_t sum_gap_counts = base_counts['-'] + base_counts['.'];
//  if (sum_base_counts > sum_gap_counts) {
//    std::vector<int8_t> bases = {'A', 'C', 'T', 'G'};
//    int8_t max_base = 'A';
//    for (int32_t j=0; j<bases.size(); j++) {
//      if (base_counts[bases[j]] > base_counts[max_base]) max_base = bases[j];
//    }
//    return std::string(1, max_base);
//  }
//
//  return std::string("");
}

int MajorityVoteFromMSA(std::vector<std::string> &msa, std::string &consensus) {
  if (msa.size() == 0) { return 1; }

  int64_t seq_len = msa[0].size();
  std::stringstream ss;

  for (int64_t i=0; i<seq_len; i++) {
    // Count occurrences for the column.
    int32_t base_counts[256] = {0};
    for (int32_t j=1; j<(msa.size()-1); j++) {
      base_counts[toupper(msa[j][i])] += 1;
    }

    int64_t sum_base_counts = base_counts['A'] + base_counts['C'] + base_counts['T'] + base_counts['G'];
    int64_t sum_gap_counts = base_counts['-'] + base_counts['.'];
    if (sum_base_counts > sum_gap_counts) {
      std::vector<int8_t> bases = {'A', 'C', 'T', 'G'};
      int8_t max_base = 'A';
      for (int32_t j=0; j<bases.size(); j++) {
        if (base_counts[bases[j]] > base_counts[max_base]) max_base = bases[j];
      }
      ss << (char) max_base;
    }
  }

  consensus = ss.str();

  return 0;
}

void CreateConsensus(const ProgramParameters &parameters, const SingleSequence *contig, std::vector<const SingleSequence *> &ctg_alns, std::map<const SingleSequence *, int64_t> &aln_lens_on_ref, std::string &ret_consensus, FILE *fp_out_cons) {
  std::stringstream ss_cons;

  if (fp_out_cons) {
    fprintf (fp_out_cons, ">Consensus-%s\n", contig->get_header());
    fflush (fp_out_cons);
  }

  int64_t num_windows = ceil((float) contig->get_sequence_length() / (float) parameters.window_len);
  LOG_DEBUG ("current_contig->get_sequence_length() = %ld, parameters.window_len = %ld, num_windows = %ld\n", contig->get_sequence_length(), parameters.window_len, num_windows);

  // Build the interval tree for fast overlap calculation.
  std::vector<IntervalSS> aln_intervals;
  for (int64_t i=0; i<ctg_alns.size(); i++) {
    int64_t aln_start = ctg_alns[i]->get_aln().get_pos() - 1;
    int64_t aln_end = aln_start + ctg_alns[i]->get_aln().GetReferenceLengthFromCigar() - 1;
    aln_intervals.push_back(IntervalSS(aln_start, aln_end, ctg_alns[i]));
  }
  IntervalTreeSS aln_interval_tree(aln_intervals);

  // For realignment, we need a vector of new alignment objects which will update the existing ones.
  std::map<const SingleSequence *, SequenceAlignment> realigns;
  if (parameters.do_realign) {
//    realigns.resize(ctg_alns.size());
    for (int64_t i=0; i<ctg_alns.size(); i++) {
      realigns[ctg_alns[i]] = SequenceAlignment();
      SequenceAlignment &r = realigns[ctg_alns[i]];
//      r.set_mapq(ctg_alns[i]->get_aln().get_mapq());
//      r.set_flag(ctg_alns[i]->get_aln().get_flag());
//      r.set
      r.CopyFrom(ctg_alns[i]->get_aln());
      r.cigar().clear();
      r.set_pos(0);
      r.set_as(0);
      r.set_evalue(0.0);
      r.optional().clear();
    }
  }

  // Process the genome in windows, but also process windows in batches. Each batch is processed in multiple threads,
  // then the results are collected and output to file. After that, a new batch is loaded.
  for (int64_t window_batch_start = parameters.start_window, num_batches = 0; window_batch_start < num_windows && (parameters.num_batches < 0 || num_batches < parameters.num_batches); window_batch_start += parameters.batch_of_windows, num_batches++) {
    std::vector<std::string> consensus_windows;
    consensus_windows.resize(parameters.batch_of_windows);
    int64_t windows_to_process = std::min(parameters.batch_of_windows, num_windows - window_batch_start);

    #pragma omp parallel for num_threads(parameters.num_threads) schedule(dynamic, 1)
    for (int64_t id_in_batch = 0; id_in_batch < windows_to_process; id_in_batch += 1) {

       int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
       int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin) - 1;
       int32_t thread_id = omp_get_thread_num();

       if (thread_id == 0) { LOG_MEDHIGH("\r(thread_id = %d) Processing window: %ld bp to %ld bp (%.2f%%)", thread_id, window_start, window_end, 100.0 * ((float) window_start / (float) contig->get_data_length())); }

       // Cut a window out of all aligned sequences. This will be fed to an MSA algorithm.
       std::vector<std::string> windows_for_msa;
       std::vector<std::string> quals_for_msa;
       std::vector<uint32_t> starts_for_msa;
       std::vector<uint32_t> ends_for_msa;
       std::vector<const SingleSequence *> refs_for_msa;

       // When realigning reads, we cannot use the QV filtering because chunks of reads would not get realigned.
       double qv_threshold = (parameters.do_realign) ? -1.0 : parameters.qv_threshold;
       ExtractWindowFromAlns(contig, ctg_alns, aln_lens_on_ref, aln_interval_tree, window_start, window_end, qv_threshold, windows_for_msa, quals_for_msa, refs_for_msa, starts_for_msa, ends_for_msa, NULL);
       if (thread_id == 0) { LOG_MEDHIGH_NOHEADER(", coverage: %ldx", windows_for_msa.size()) }
       if (windows_for_msa.size() == 0) {
         consensus_windows[id_in_batch] = "";
         ERROR_REPORT(ERR_UNEXPECTED_VALUE, "windows_for_msa.size() == 0!");
       } else if (quals_for_msa.size() != windows_for_msa.size()) {
         ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Quality values not specified for input reads! Please use FASTQ or another format which provides qualities.");
       }

       // Sort the sequences in the window by their length. This is a temporary workaround to be able to generate the realigned reads.
//       std::vector<uint32_t> indices;
//       if (parameters.do_realign) {
//         prepare_indices1(indices, windows_for_msa, true);
//       }

       auto graph = construct_partial_order_graph(windows_for_msa, quals_for_msa, starts_for_msa, ends_for_msa,
                                                  SPOA::AlignmentParams(parameters.match, parameters.mismatch,
                                                  parameters.gap_open, parameters.gap_ext, (SPOA::AlignmentType) parameters.aln_type));

       if (windows_for_msa.size() <= 2) {  // In case the coverage is too low, just pick the first sequence in the window.
//          int64_t temp_window_end = std::min((int64_t) window_end, (int64_t) (contig->get_sequence_length()-1));
//          consensus_windows[id_in_batch] = GetSubstring((char *) (contig->get_data() + window_start), (temp_window_end - window_start + 1));
           consensus_windows[id_in_batch] = windows_for_msa[0];
       } else {
           consensus_windows[id_in_batch] = graph->generate_consensus();

           std::vector<std::string> msa;
           graph->generate_msa(msa, true);
//           for (int64_t i=0; i<msa.size(); i++) {
//             printf ("%s\n", msa[i].c_str());
//           }

//           int32_t overhang_min_threshold = 3;
           FilterOverhangsFromMsa(msa, consensus_windows[id_in_batch]);
//           MajorityVoteFromMSA(msa, consensus_windows[id_in_batch]);
       }

       if (parameters.do_realign) {
         std::vector<std::string> msa;
         graph->generate_msa(msa, true);
         // Sequence at msa[0] is the reference sequence.
         for (int64_t i=1; i<msa.size(); i++) {
           auto seq = refs_for_msa[i];
           ConvertMSAToAln(msa[0], msa[i], window_start, realigns[seq]);
         }

//         FILE *fp_test = fopen("temp/test.msa", "w");
////         fprintf (fp_test, "%s\n", consensus_windows[id_in_batch].c_str());
//         for (int64_t i1=0; i1<msa.size(); i1++) {
////             fprintf (fp_test, ">ID_%d_MSA_Window_%d\n%s\n", (id_in_batch), i1, msa[i1].c_str());
//           fprintf (fp_test, "%s\n", msa[i1].c_str());
//           fflush(fp_test);
//         }
//         fprintf (fp_test, "\n");
//         fclose(fp_test);
//         exit(1);
       }
    }

     LOG_MEDHIGH_NOHEADER("\n");
     LOG_MEDHIGH("Batch checkpoint: Performed consensus on all windows, joining the windows now.\n");
     for (int64_t id_in_batch = 0; id_in_batch < parameters.batch_of_windows && id_in_batch < num_windows; id_in_batch += 1) {
       int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
       int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin);

       if (id_in_batch == 0) {
         ss_cons << consensus_windows[id_in_batch];
         if (fp_out_cons) {
           fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
           fflush(fp_out_cons);
         }

       } else {
         if (parameters.win_ovl_margin <= 0.0) {  // If overlapping windows is turned off.
           ss_cons << consensus_windows[id_in_batch];
           if (fp_out_cons) {
             fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
             fflush(fp_out_cons);
           }

         } else {     // Otherwise, do the overlap alignment.
           fprintf (stderr, "Overlapping windows.\n");
           fflush(stderr);

           std::string trimmed_window = consensus_windows[id_in_batch-1].substr((1.0 - parameters.win_ovl_margin * 3) * consensus_windows[id_in_batch-1].size());

           std::vector<std::string> windows_for_alignment = {trimmed_window, consensus_windows[id_in_batch]};
           std::vector<std::string> msa;

           SPOA::generate_msa(msa, windows_for_alignment, SPOA::AlignmentParams(1, -1, -1, -1, SPOA::AlignmentType::kOV), false);

           std::stringstream ss_clipped_window;
           int32_t clip_pos = 0;
           int32_t trimmed_id = 0, curr_window_id = 1;
           for (clip_pos=(msa[trimmed_id].size()-1); clip_pos>=0 && msa[trimmed_id][clip_pos] == '-'; clip_pos--);

           for (clip_pos++; clip_pos<msa[curr_window_id].size(); clip_pos++) {
             if (msa[curr_window_id][clip_pos] != '-') { ss_clipped_window << msa[curr_window_id][clip_pos]; }
           }
           std::string clipped_window = ss_clipped_window.str();

           if (clipped_window.size() > 0) {
             ss_cons << clipped_window;
             if (fp_out_cons) {
               fprintf (fp_out_cons, "%s", clipped_window.c_str());
               fflush(fp_out_cons);
             }
//               printf ("[good] window_start = %ld, window_end = %ld, clipped_window.size() = %ld\n", window_start, window_end, clipped_window.size());
//               fflush(stdout);
           } else {
             printf ("\n");
             printf ("[bad] window_start = %ld, window_end = %ld, clipped_window.size() = %ld\n", window_start, window_end, clipped_window.size());
             printf ("\n");
             for (int32_t i2=0; i2<windows_for_alignment.size(); i2++) {
               printf ("%s\n\n", windows_for_alignment[i2].c_str());
             }
             printf ("\n");
             printf ("Alignment:\n\n");
             fflush(stdout);
             for (int32_t i2=0; i2<msa.size(); i2++) {
               printf ("%s\n\n", msa[i2].c_str());
             }
             printf ("\n");

             fflush(stdout);
             exit(1);
           }

//             exit(1);
         }
       }
     }

     LOG_MEDHIGH_NOHEADER("\n");
     LOG_MEDHIGH("Batch checkpoint: Processed %ld windows and exported the consensus.\n", parameters.batch_of_windows);
  }

  if (parameters.do_realign) {
    for (int64_t i=0; i<ctg_alns.size(); i++) {
      const SequenceAlignment &r = realigns[ctg_alns[i]];
      ((SingleSequence *) ctg_alns[i])->aln().CopyFrom(r);
      std::string sam_line = ctg_alns[i]->MakeSAMLine();
      fprintf (stdout, "%s\n", sam_line.c_str());
    }
  }

  ret_consensus = ss_cons.str();
  if (fp_out_cons) {
    fprintf (fp_out_cons, "\n");
  }
}
