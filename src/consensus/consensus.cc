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
#include "libs/edlib.h"
#include "pileup.h"

// #define WINDOW_OUTPUT_IN_FASTQ

void CreateConsensus(const ProgramParameters &parameters, int32_t num_window_threads, const SingleSequence *contig, const std::vector<SingleSequence *> &ctg_alns, std::map<const SingleSequence *, int64_t> &aln_lens_on_ref, std::string &ret_consensus, FILE *fp_out_cons) {
  std::stringstream ss_cons;

  if (parameters.do_erc == false && fp_out_cons) {
    fprintf (fp_out_cons, ">Consensus_%s\n", contig->get_header());
    fflush (fp_out_cons);
  }

  int64_t num_windows = ceil((float) contig->get_sequence_length() / (float) parameters.window_len);
  if (parameters.do_erc == false) {
    LOG_DEBUG ("current_contig->get_sequence_length() = %ld, parameters.window_len = %ld, num_windows = %ld\n", contig->get_sequence_length(), parameters.window_len, num_windows);
  }

  // Build the interval tree for fast overlap calculation.
  std::vector<IntervalSS> aln_intervals;
  for (int64_t i=0; i<ctg_alns.size(); i++) {
    int64_t aln_start = ctg_alns[i]->get_aln().get_pos() - 1;
    int64_t aln_end = aln_start + ctg_alns[i]->get_aln().GetReferenceLengthFromCigar() - 1;
    aln_intervals.push_back(IntervalSS(aln_start, aln_end, ctg_alns[i]));
  }
  IntervalTreeSS aln_interval_tree(aln_intervals);

  // Process the genome in windows, but also process windows in batches. Each batch is processed in multiple threads,
  // then the results are collected and output to file. After that, a new batch is loaded.
  for (int64_t window_batch_start = parameters.start_window, num_batches = 0; window_batch_start < num_windows && (parameters.num_batches < 0 || num_batches < parameters.num_batches); window_batch_start += parameters.batch_of_windows, num_batches++) {
    std::vector<std::string> consensus_windows;
    consensus_windows.resize(parameters.batch_of_windows);
    int64_t windows_to_process = std::min(parameters.batch_of_windows, num_windows - window_batch_start);

//    #pragma omp parallel for num_threads(parameters.num_threads) schedule(dynamic, 1)
    #pragma omp parallel for num_threads(num_window_threads) schedule(dynamic, 1)
    for (int64_t id_in_batch = 0; id_in_batch < windows_to_process; id_in_batch += 1) {

       int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
       int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin) - 1;
       int32_t thread_id = omp_get_thread_num();

       if (parameters.do_erc == false && thread_id == 0) {
         LOG_MEDHIGH("\r(thread_id = %d) Processing window: %ld bp to %ld bp (%.2f%%)", thread_id, window_start, window_end, 100.0 * ((float) window_start / (float) contig->get_data_length()));
       }

       // Cut a window out of all aligned sequences. This will be fed to an MSA algorithm.
       std::vector<std::string> windows_for_msa;
       std::vector<std::string> quals_for_msa;
       std::vector<uint32_t> starts_for_msa;
       std::vector<uint32_t> ends_for_msa;
       std::vector<uint32_t> starts_on_read;
       std::vector<uint32_t> ends_on_read;
       std::vector<const SingleSequence *> refs_for_msa;

       // When realigning reads, we cannot use the QV filtering because chunks of reads would not get realigned.
       // Realignment mode has been removed for now, so the upper comment can be ignored. Left here for future reference, so it doesn't get forgotten.
       ExtractWindowFromAlns(contig, ctg_alns, aln_lens_on_ref, aln_interval_tree,
                             window_start, window_end, parameters.qv_threshold, parameters.use_contig_qvs,
                             windows_for_msa, quals_for_msa, refs_for_msa,
                             starts_for_msa, ends_for_msa, starts_on_read, ends_on_read, NULL);
       //       ExtractWindowFromAlns(contig, ctg_alns, aln_lens_on_ref, aln_interval_tree, window_start, window_end, qv_threshold, windows_for_msa, quals_for_msa, refs_for_msa, starts_for_msa, ends_for_msa, NULL);
       if (parameters.do_erc == false && thread_id == 0) { LOG_MEDHIGH_NOHEADER(", coverage: %ldx", windows_for_msa.size()) }
       if (windows_for_msa.size() == 0) {
         consensus_windows[id_in_batch] = "";
         ERROR_REPORT(ERR_UNEXPECTED_VALUE, "windows_for_msa.size() == 0!");
       } else if (quals_for_msa.size() != windows_for_msa.size()) {
         ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Quality values not specified for input reads! Please use FASTQ or another format which provides qualities.");
       }

       auto graph = construct_partial_order_graph(windows_for_msa, quals_for_msa, starts_for_msa, ends_for_msa,
                                                  SPOA::AlignmentParams(parameters.match, parameters.mismatch,
                                                  parameters.gap_open, parameters.gap_ext, (SPOA::AlignmentType) parameters.aln_type));

       if (windows_for_msa.size() <= 2) {  // In case the coverage is too low, just pick the first sequence in the window.
           consensus_windows[id_in_batch] = windows_for_msa[0];
       } else {
         std::vector<uint32_t> coverages;
           std::string cons_window = graph->generate_consensus(coverages);

           int32_t start_offset = 0, end_offset = cons_window.size() - 1;
           for (;start_offset<cons_window.size(); start_offset++) {
             if (coverages[start_offset] >= ((windows_for_msa.size() - 1) / 2)) { break; }
           }
           for (; end_offset >= 0; end_offset--) {
             if (coverages[start_offset] >= ((windows_for_msa.size() - 1) / 2)) { break; }
           }

           consensus_windows[id_in_batch] = cons_window.substr(start_offset, (end_offset - start_offset + 1));

//           std::vector<std::string> msa;
//           graph->generate_msa(msa, true);
//           for (int64_t i=0; i<msa.size(); i++) { printf ("%s\n", msa[i].c_str()); }
//           FilterOverhangsFromMsa(msa, consensus_windows[id_in_batch]);
//           MajorityVoteFromMSA(msa, consensus_windows[id_in_batch]);
       }

//       std::vector<std::string> debug_msa;
//       graph->generate_msa(debug_msa, true);
//       FILE *fp_debug = fopen("temp/debug.msa", "w");
//       for (int64_t i=0; i<debug_msa.size(); i++) {
//         fprintf (fp_debug, "%s", debug_msa[i].c_str());
//         if ((i + 1) < debug_msa.size()) {
//           fprintf (fp_debug, "\t%s", refs_for_msa[i]->get_header());
//         }
//         fprintf (fp_debug, "\n");
//       }
//       fclose(fp_debug);

    }

    if (parameters.do_erc == false) {
      LOG_MEDHIGH_NOHEADER("\n");
      LOG_MEDHIGH("Batch checkpoint: Performed consensus on all windows, joining the windows now.\n");
    }

    for (int64_t id_in_batch = 0; id_in_batch < parameters.batch_of_windows && id_in_batch < num_windows; id_in_batch += 1) {
      int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
      int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin);
      ss_cons << consensus_windows[id_in_batch];
      if (parameters.do_erc == false && fp_out_cons) {
        fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
        fflush(fp_out_cons);
      }
    }

    // Deprecated, used for window overlapping.
//     LOG_MEDHIGH_NOHEADER("\n");
//     LOG_MEDHIGH("Batch checkpoint: Performed consensus on all windows, joining the windows now.\n");
//     for (int64_t id_in_batch = 0; id_in_batch < parameters.batch_of_windows && id_in_batch < num_windows; id_in_batch += 1) {
//       int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
//       int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin);
//
//       if (id_in_batch == 0) {
//         ss_cons << consensus_windows[id_in_batch];
//         if (fp_out_cons) {
//           fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
//           fflush(fp_out_cons);
//         }
//
//       } else {
//         if (parameters.win_ovl_margin <= 0.0) {  // If overlapping windows is turned off.
//           ss_cons << consensus_windows[id_in_batch];
//           if (fp_out_cons) {
//             fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
//             fflush(fp_out_cons);
//           }
//
//         } else {     // Otherwise, do the overlap alignment.
//           fprintf (stderr, "Overlapping windows.\n");
//           fflush(stderr);
//
//           std::string trimmed_window = consensus_windows[id_in_batch-1].substr((1.0 - parameters.win_ovl_margin * 3) * consensus_windows[id_in_batch-1].size());
//
//           std::vector<std::string> windows_for_alignment = {trimmed_window, consensus_windows[id_in_batch]};
//           std::vector<std::string> msa;
//
//           SPOA::generate_msa(msa, windows_for_alignment, SPOA::AlignmentParams(1, -1, -1, -1, SPOA::AlignmentType::kOV), false);
//
//           std::stringstream ss_clipped_window;
//           int32_t clip_pos = 0;
//           int32_t trimmed_id = 0, curr_window_id = 1;
//           for (clip_pos=(msa[trimmed_id].size()-1); clip_pos>=0 && msa[trimmed_id][clip_pos] == '-'; clip_pos--);
//
//           for (clip_pos++; clip_pos<msa[curr_window_id].size(); clip_pos++) {
//             if (msa[curr_window_id][clip_pos] != '-') { ss_clipped_window << msa[curr_window_id][clip_pos]; }
//           }
//           std::string clipped_window = ss_clipped_window.str();
//
//           if (clipped_window.size() > 0) {
//             ss_cons << clipped_window;
//             if (fp_out_cons) {
//               fprintf (fp_out_cons, "%s", clipped_window.c_str());
//               fflush(fp_out_cons);
//             }
////               printf ("[good] window_start = %ld, window_end = %ld, clipped_window.size() = %ld\n", window_start, window_end, clipped_window.size());
////               fflush(stdout);
//           } else {
//             printf ("\n");
//             printf ("[bad] window_start = %ld, window_end = %ld, clipped_window.size() = %ld\n", window_start, window_end, clipped_window.size());
//             printf ("\n");
//             for (int32_t i2=0; i2<windows_for_alignment.size(); i2++) {
//               printf ("%s\n\n", windows_for_alignment[i2].c_str());
//             }
//             printf ("\n");
//             printf ("Alignment:\n\n");
//             fflush(stdout);
//             for (int32_t i2=0; i2<msa.size(); i2++) {
//               printf ("%s\n\n", msa[i2].c_str());
//             }
//             printf ("\n");
//
//             fflush(stdout);
//             exit(1);
//           }
//
////             exit(1);
//         }
//       }
//     }

//     LOG_MEDHIGH_NOHEADER("\n");
    if (parameters.do_erc == false) {
      LOG_MEDHIGH("Batch checkpoint: Processed %ld windows and exported the consensus.\n", parameters.batch_of_windows);
    }
  }

  ret_consensus = ss_cons.str();
  if (parameters.do_erc == false && fp_out_cons) {
    fprintf (fp_out_cons, "\n");
  }
}
