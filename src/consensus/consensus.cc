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

typedef Interval<const SingleSequence *> IntervalSS;
typedef IntervalTree<const SingleSequence *> IntervalTreeSS;

// #define WINDOW_OUTPUT_IN_FASTQ

//std::vector<size_t> soort(const std::vector<std::string>& windows) {
//    std::vector<size_t> indices(windows.size());
//    std::iota(begin(indices), end(indices), static_cast<size_t>(0));
//
//    std::sort(
//        begin(indices), end(indices),
//        [&](size_t a, size_t b) { return windows[a].size() > windows[b].size(); }
//    );
//    return indices;
//}

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
                           IntervalTreeSS &aln_interval_tree, int64_t window_start, int64_t window_end, double qv_threshold, std::vector<std::string> &window_seqs, std::vector<std::string> &window_qv,
                           std::vector<uint32_t> &window_starts, std::vector<uint32_t> &window_ends, FILE *fp_window) {
  if (window_start > window_end) {
    return;
  }

  int64_t temp_window_end = std::min((int64_t) window_end, (int64_t) (contig->get_sequence_length()-1));
  window_seqs.push_back(GetSubstring((char *) (contig->get_data() + window_start), (temp_window_end - window_start + 1)));
  std::string dummy_quals((temp_window_end - window_start + 1), '!');
  window_qv.push_back(dummy_quals);
  window_starts.push_back(0);
  window_ends.push_back(temp_window_end - window_start);

//  fflush(stdout);
//  fflush(stderr);
//  printf ("\n");
//  printf ("window_start = %ld, window_end = %ld, temp_window_end = %ld\n", window_start, window_end, temp_window_end);
//  fflush(stdout);
//
//  printf ("seq_start_in_window = %u, seq_end_in_window = %u, seq_len = %ld, qual_len = %ld, %s, %s\n",
//          window_starts.back(), window_ends.back(), window_seqs.back().size(), window_qv.back().size(),
//          window_seqs.back().c_str(), window_qv.back().c_str());
//  fflush(stdout);

  std::vector<IntervalSS> intervals;
  aln_interval_tree.findOverlapping(window_start, temp_window_end, intervals);
  for (int64_t i=0; i<intervals.size(); i++) {
    auto seq = intervals[i].value;
    auto aln = seq->get_aln();

    int64_t start_cig_id = 0, end_cig_id = 0;
//    printf ("Looking for the start: window_start = %ld\n", window_start);
//    fflush(stdout);
    int64_t start_seq = aln.FindBasePositionOnRead(window_start, &start_cig_id);
//    printf ("Looking for the end: temp_window_end = %ld\n", temp_window_end);
//    fflush(stdout);
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
      window_seqs.push_back(seq_data);
      window_starts.push_back(seq_start_in_window);
      window_ends.push_back(seq_end_in_window);
      window_qv.push_back(seq_qual);
    }

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
//  GroupAlignmentsToContigs(alns, parameters.qv_threshold, ctg_names, all_ctg_alns);
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

  // Debug output of alternate contigs, aligned to the raw contig (input sequence), in SAM format.
  // Deprecated.
//  FILE *fp_alt_contig_path = NULL;
//  if (parameters.temp_alt_contig_path != "") {
//    fp_alt_contig_path = fopen(parameters.temp_alt_contig_path.c_str(), "w");
//    fprintf (fp_alt_contig_path, "@HD\tVN:1.0\tSO:unknown\n");
//    for (int32_t i=0; i<ctg_names.size(); i++) {
//      fprintf (fp_alt_contig_path, "@SQ\tSN:%s\tLN:%ld\n", ctg_names[i].c_str(), rname_to_seq[ctg_names[i]]->get_data_length());
//    }
//    fprintf (fp_alt_contig_path, "@PG\tID:consise PN:consise\n");
//  }

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
    CreateConsensus(parameters, contig, ctg_alns, aln_lens_on_ref, consensus, fp_out_cons);
    fclose(fp_out_cons);

    ///////////////////////////////////////
    LOG_MEDHIGH_NOHEADER("\n");
    LOG_ALL("Processed %ld bp of %ld bp (100.00%%)\n", contig->get_data_length(), contig->get_data_length());
    LOG_MEDHIGH_NOHEADER("\n");
  }

//  if (fp_alt_contig_path) { fclose(fp_alt_contig_path); }

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

//  FILE *fp_test = fopen("temp/test.fasta", "w");

  // Process the genome in windows, but also process windows in batches. Each batch is processed in multiple threads,
  // then the results are collected and output to file. After that, a new batch is loaded.
  for (int64_t window_batch_start = parameters.start_window, num_batches = 0; window_batch_start < num_windows && (parameters.num_batches < 0 || num_batches < parameters.num_batches); window_batch_start += parameters.batch_of_windows, num_batches++) {
    std::vector<std::string> consensus_windows;
    consensus_windows.resize(parameters.batch_of_windows);
    int64_t windows_to_process = std::min(parameters.batch_of_windows, num_windows - window_batch_start);

//    windows_to_process = 49;
//    #pragma omp parallel for num_threads(parameters.num_threads) schedule(dynamic, 1)
//     for (int64_t id_in_batch = 48; id_in_batch < windows_to_process; id_in_batch += 1) {
      #pragma omp parallel for num_threads(parameters.num_threads) schedule(dynamic, 1)
       for (int64_t id_in_batch = 0; id_in_batch < windows_to_process; id_in_batch += 1) {

//       if ((id_in_batch + 1) == (windows_to_process)) { break; }

       int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
       int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin) - 1;
       int32_t thread_id = omp_get_thread_num();
//       int32_t thread_id = 0;

       if (thread_id == 0) { LOG_MEDHIGH("\r(thread_id = %d) Processing window: %ld bp to %ld bp (%.2f%%)", thread_id, window_start, window_end, 100.0 * ((float) window_start / (float) contig->get_data_length())); }

       // Cut a window out of all aligned sequences. This will be fed to an MSA algorithm.
       std::vector<std::string> windows_for_msa;
       std::vector<std::string> quals_for_msa;
       std::vector<uint32_t> starts_for_msa;
       std::vector<uint32_t> ends_for_msa;

//       std::vector<std::string> msa;

       // Chosing the MSA algorithm, and running the consensus on the window.
       if (parameters.msa == "poa") {
         ExtractWindowFromAlns(contig, ctg_alns, aln_lens_on_ref, aln_interval_tree, window_start, window_end, parameters.qv_threshold, windows_for_msa, quals_for_msa, starts_for_msa, ends_for_msa, NULL);

         if (thread_id == 0) { LOG_MEDHIGH_NOHEADER(", coverage: %ldx", windows_for_msa.size()) }

         if (windows_for_msa.size() == 0) {
             consensus_windows[id_in_batch] = "";
             ERROR_REPORT(ERR_UNEXPECTED_VALUE, "windows_for_msa.size() == 0!");
         } else if (windows_for_msa.size() <= 2) {
 //          int64_t temp_window_end = std::min((int64_t) window_end, (int64_t) (contig->get_sequence_length()-1));
 //          consensus_windows[id_in_batch] = GetSubstring((char *) (contig->get_data() + window_start), (temp_window_end - window_start + 1));
             consensus_windows[id_in_batch] = windows_for_msa[0];

         } else {
           if (quals_for_msa.size() > 0) {
//             printf ("id_in_batch = %ld\n", id_in_batch);
//             for (int64_t i1=0; i1<windows_for_msa.size(); i1++) {
//               printf ("[i1 = %ld] %u, %u\n", i1, starts_for_msa[i1], ends_for_msa[i1]);
//               fflush(stdout);
//             }
//             printf ("\n");
//             fflush(stdout);

             consensus_windows[id_in_batch] = SPOA::generate_consensus(windows_for_msa, quals_for_msa, starts_for_msa, ends_for_msa,
                                                                       SPOA::AlignmentParams(parameters.match, parameters.mismatch,
                                                                       parameters.gap_open, parameters.gap_ext, (SPOA::AlignmentType) parameters.aln_type));

//             consensus_windows[id_in_batch] = SPOA::generate_consensus(msa, windows_for_msa, quals_for_msa, starts_for_msa, ends_for_msa,
//                                                                       SPOA::AlignmentParams(parameters.match, parameters.mismatch,
//                                                                       parameters.gap_open, parameters.gap_ext, (SPOA::AlignmentType) parameters.aln_type));
//             consensus_windows[id_in_batch] = SPOA::generate_consensus(windows_for_msa, quals_for_msa,
//                                                                       SPOA::AlignmentParams(parameters.match, parameters.mismatch,
//                                                                       parameters.gap_open, parameters.gap_ext, (SPOA::AlignmentType) parameters.aln_type));

           } else {
             consensus_windows[id_in_batch] = SPOA::generate_consensus(windows_for_msa, SPOA::AlignmentParams(parameters.match,
                                                                                                        parameters.mismatch, parameters.gap_open, parameters.gap_ext, (SPOA::AlignmentType) parameters.aln_type), false);
           }

//           fprintf (fp_test, ">ID_%d_Consensus window %ld to %ld\n%s\n", (id_in_batch), window_start, window_end, consensus_windows[id_in_batch].c_str());

////           for (int64_t i1=0; i1<windows_for_msa.size(); i1++) {
////             fprintf (fp_test, ">ID_%d_Window_%d\n%s\n", (id_in_batch), i1, windows_for_msa[i1].c_str());
////             fflush(fp_test);
////           }
////
////           std::vector<std::string> msa;
////           SPOA::generate_msa(msa, windows_for_msa, quals_for_msa, SPOA::AlignmentParams(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_ext, (SPOA::AlignmentType) parameters.aln_type), true);
//           for (int64_t i1=0; i1<msa.size(); i1++) {
////             fprintf (fp_test, ">ID_%d_MSA_Window_%d\n%s\n", (id_in_batch), i1, msa[i1].c_str());
//             fprintf (fp_test, "%s\n", msa[i1].c_str());
//             fflush(fp_test);
//           }
//           fprintf (fp_test, "\n");

         }

//         fprintf (fp_test, ">%d\n%s\n", (id_in_batch), consensus_windows[id_in_batch].c_str());
//         fflush(fp_test);

       } else {
         FILE *fp_window = NULL;
         std::string window_path = FormatString("%s.%ld", parameters.temp_window_path.c_str(), thread_id);
         if (parameters.temp_window_path != "") {
           fp_window = fopen(window_path.c_str(), "w");
         }
         if (fp_window == NULL) {
           ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Window file not opened!\n");
         }
         ExtractWindowFromAlns(contig, ctg_alns, aln_lens_on_ref, aln_interval_tree, window_start, window_end, parameters.qv_threshold, windows_for_msa, quals_for_msa, starts_for_msa, ends_for_msa, fp_window);
         fclose(fp_window);
         RunMSAFromSystemLocal(parameters, window_path, consensus_windows[id_in_batch]);
       }
     }

//     fclose(fp_test);
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

  ret_consensus = ss_cons.str();
  if (fp_out_cons) {
    fprintf (fp_out_cons, "\n");
  }
}



int MajorityVoteFromMSALocal(std::string pir_path, std::string *cons) {
  SequenceFile pir(SEQ_FORMAT_FASTQ, pir_path, false);

  const SequenceVector& seqs = pir.get_sequences();

  if (seqs.size() == 0) { return 1; }


  for (int64_t i=1; i<pir.get_sequences().size(); i++) {
    if (seqs[i]->get_data_length() != seqs[i-1]->get_data_length()) {
      ERROR_REPORT(ERR_FILE_DEFORMED_FORMAT, "MSA sequences not of correct length in file %s!\n", pir_path.c_str());
      return 2;
    }
  }

  int64_t seq_len = seqs[0]->get_data_length();
  *cons = "";
  std::stringstream ss;

  for (int64_t i=0; i<seq_len; i++) {
    // Count occurrences for the column.
    int32_t base_counts[256] = {0};
    for (int32_t j=0; j<seqs.size(); j++) {
      base_counts[toupper(seqs[j]->get_data()[i])] += 1;
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

  *cons = ss.str();

  return 0;
}

int RunMSAFromSystemLocal(const ProgramParameters &parameters, std::string window_path, std::string &cons) {
  std::string msa_path = FormatString("%s.msa", window_path.c_str());
  if (parameters.msa == "mafft") {
    //
//    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --localpair --maxiterate 1000 --quiet %s > %s", // AlignedBases           48306(99.60%)       47762(100.00%)  AvgIdentity                    96.33                96.33
//    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --localpair --op 0 --ep 1 --maxiterate 1000 --quiet %s > %s", // AlignedBases           48306(99.60%)       47353(100.00%)  AvgIdentity                    97.03                97.03
//    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --localpair --lop 0 --lexp 1 --quiet %s > %s", // AlignedBases           48306(99.60%)       47482(100.00%)  AvgIdentity                    96.87                96.87
// Ovaj je dobar    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --localpair --op 0 --ep 1 --quiet %s > %s", // AlignedBases           48306(99.60%)       47482(100.00%)  AvgIdentity                    96.87                96.87
//    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --globalpair --op 0 --ep 1 --quiet %s > %s", // AlignedBases           48306(99.60%)       47482(100.00%)  AvgIdentity                    96.87                96.87

 // Trenutno najbolji rezultat:
//    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --retree 1 --maxiterate 0 --nofft --genafpair --op 0 --ep 1 --quiet %s > %s", // AlignedBases           48306(99.60%)       47482(100.00%)  AvgIdentity                    96.87                96.87

//    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --retree 1 --maxiterate 0 --nofft --op 0 --ep 1 --quiet %s > %s", // AlignedBases           48306(99.60%)       47482(100.00%)  AvgIdentity                    96.87                96.87
//                        parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), window_path.c_str(), msa_path.c_str()).c_str());

    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --op 0 --ep 1 --quiet %s > %s", // AlignedBases           48306(99.60%)       47482(100.00%)  AvgIdentity                    96.87                96.87
                        parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), window_path.c_str(), msa_path.c_str()).c_str());

//    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --op 0 --ep 1 --quiet %s > %s.noleft", // AlignedBases           48306(99.60%)       47482(100.00%)  AvgIdentity                    96.87                96.87
//                        parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), window_path.c_str(), msa_path.c_str()).c_str());
//    int32_t rc1 = system(FormatString("scripts/test-leftalign/convert_msa2sam.py %s.noleft %s",
//                        msa_path.c_str(), msa_path.c_str()).c_str());

  } else if (parameters.msa == "poav2") {
    int32_t rc = system(FormatString("%s/%s/poa -do_local -do_progressive -read_fasta %s -pir %s %s/../settings/all1-poav2.mat",
                        parameters.program_folder.c_str(), parameters.poav2_folder.c_str(), window_path.c_str(), msa_path.c_str(), parameters.program_folder.c_str()).c_str());
  } else {
    ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Unrecognized MSA option '%s'! Exiting.",parameters. msa.c_str());
    return 1;
  }

  MajorityVoteFromMSALocal(msa_path, &cons);

  return 0;
}
