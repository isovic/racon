/*
 * consensus-sam.cc
 *
 *  Created on: Nov 27, 2016
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
#include "tictoc.h"

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

// Used when alignments are specified via a SAM file.
int ConsensusFromAln(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &alns) {
  LOG_MEDHIGH("Running consensus - directly from alignments.\n");

  int32_t num_read_threads = (parameters.do_erc) ? (parameters.num_threads) : 1;
  int32_t num_window_threads = (!parameters.do_erc) ? (parameters.num_threads) : 1;

  std::vector<std::string> ctg_names;
  std::map<std::string, std::vector<const SingleSequence *> > all_ctg_alns;

  // Separate alignments into groups for each contig.
  // Alignments which are called unmapped will be skipped in this step.
  // Also, alignments are filtered by the base quality if available.
  LOG_MEDHIGH("Separating alignments to individual contigs.\n");
  GroupAlignmentsToContigs(alns, -1.0, ctg_names, all_ctg_alns);

  // Hash the sequences by their name.
  std::map<std::string, const SingleSequence *> rname_to_seq;
  for (int32_t i=0; i<contigs.get_sequences().size(); i++) {
    rname_to_seq[contigs.get_sequences()[i]->get_header()] = contigs.get_sequences()[i];
    rname_to_seq[TrimToFirstSpace(contigs.get_sequences()[i]->get_header())] = contigs.get_sequences()[i];
  }

  // Verbose.
  // If we are doing error correction, parallelization is per-read and not per-window.
  // We need to disable some of the debug info.
  if (parameters.do_erc == false) {
    LOG_MEDHIGH("In total, there are %ld contigs for consensus, each containing:\n", ctg_names.size());
    for (int32_t i=0; i<ctg_names.size(); i++) {
      LOG_MEDHIGH("\t[%ld] %s %ld alignments, contig len: %ld\n", i, ctg_names[i].c_str(), all_ctg_alns.find(ctg_names[i])->second.size(), rname_to_seq[ctg_names[i]]->get_sequence_length());
    }
  } else {
    LOG_MEDHIGH("In total, there are %ld sequences for error correction.\n", ctg_names.size());
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
  #pragma omp parallel for num_threads(num_read_threads) schedule(dynamic, 1)
  for (int32_t i=0; i<ctg_names.size(); i++) {
    int32_t thread_id = omp_get_thread_num();

    const SingleSequence *contig = rname_to_seq[ctg_names[i]];
    auto it = all_ctg_alns.find(ctg_names[i]);
    if (it == all_ctg_alns.end()) {
      FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Something strange happened. Contig name, which was extracted from alignments, cannot be found in the std::map containing those same alignments.");
      // Exits.
    }
    // Get alignments for current contig.
    std::vector<SingleSequence *> &ctg_alns = (std::vector<SingleSequence *> &) it->second;

    // This sorts ascending by the pos field.
    std::sort(ctg_alns.begin(), ctg_alns.end(), seqaln_sort_key());

    // If we are doing error correction, parallelization is per-read and not per-window.
    // We need to disable some of the debug info.
    if (parameters.do_erc == false) {
      LOG_ALL("Starting consensus for contig %ld / %ld (%.2f%%): %s\n", (i + 1), ctg_names.size(), 100.0*((float) (i + 1)) / ((float) ctg_names.size()), contig->get_header());
    }

    FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "a");
    std::string consensus;
    if (parameters.do_pileup == false) {
      if (parameters.do_erc == false) {
        CreateConsensusAln(parameters, num_window_threads, contig, ctg_alns, aln_lens_on_ref, consensus, fp_out_cons);

      } else {
        if (thread_id == 0) {
          LOG_MEDHIGH("\r(thread_id = %d) Processing contig %ld / %ld (%.2f%%), len: %10ld", thread_id, (i + 1), ctg_names.size(), 100.0f*(((float) (i)) / ((float) ctg_names.size())), contig->get_sequence_length());
        }

        CreateConsensusAln(parameters, num_window_threads, contig, ctg_alns, aln_lens_on_ref, consensus, NULL);
        #pragma omp critical
        {
          fprintf (fp_out_cons, ">Consensus_%s\n%s\n", contig->get_header(), consensus.c_str());
//          fflush (fp_out_cons);
        }
      }

    } else {
      Pileup pileup(contig, ctg_alns);
//      pileup.Verbose(stdout);
      pileup.GenerateConsensus(5, consensus);
      #pragma omp critical
      fprintf (fp_out_cons, ">Consensus_%s\n%s\n", contig->get_header(), consensus.c_str());
      #pragma omp critical
      fflush (fp_out_cons);
    }
    fclose(fp_out_cons);

    ///////////////////////////////////////
//    LOG_MEDHIGH_NOHEADER("\n");
    if (parameters.do_erc == false) {
      LOG_ALL("Processed %ld bp of %ld bp (100.00%%)\n", contig->get_data_length(), contig->get_data_length());
      LOG_MEDHIGH_NOHEADER("\n");
    }
  }

  return 0;
}

void ExtractWindowFromAlns(const SingleSequence *contig, const std::vector<SingleSequence *> &alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens,
                           IntervalTreeSS &aln_interval_tree, int64_t window_start, int64_t window_end, double qv_threshold, bool use_contig_qvs,
                           std::vector<std::string> &window_seqs, std::vector<std::string> &window_qv, std::vector<const SingleSequence *> &window_refs,
                           std::vector<uint32_t> &window_starts, std::vector<uint32_t> &window_ends,
                           std::vector<uint32_t> &starts_on_read, std::vector<uint32_t> &ends_on_read, FILE *fp_window) {
  if (window_start > window_end) {
    return;
  }

  int64_t temp_window_end = std::min((int64_t) window_end, (int64_t) (contig->get_sequence_length()-1));
  window_refs.push_back(contig);
  window_seqs.push_back(GetSubstring((char *) (contig->get_data() + window_start), (temp_window_end - window_start + 1)));
  if (use_contig_qvs == false || contig->get_quality() == NULL || contig->get_quality_length() == 0) {
    std::string dummy_quals((temp_window_end - window_start + 1), '!');
    window_qv.push_back(dummy_quals);
  } else {
    window_qv.push_back(GetSubstring((char *) (contig->get_quality() + window_start), (temp_window_end - window_start + 1)));
  }
  window_starts.push_back(0);
  window_ends.push_back(temp_window_end - window_start);
  starts_on_read.push_back(window_start);
  ends_on_read.push_back(window_end - 1);

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
      start_cig_id = 0;

    } else if (start_seq < 0) {
      fprintf (stderr, "ERROR: start_seq is < 0 and != -1! start_seq = %ld\n", start_seq); exit(1);
    }

    if (aln.get_cigar()[start_cig_id].op == 'D' || aln.get_cigar()[start_cig_id].op == 'I' || aln.get_cigar()[start_cig_id].op == 'S') {

      for (; start_cig_id < aln.get_cigar().size(); start_cig_id++) {
        if (aln.get_cigar()[start_cig_id].op == 'M' || aln.get_cigar()[start_cig_id].op == '=' || aln.get_cigar()[start_cig_id].op == 'X') { break; }
      }
      start_seq = aln.get_cigar()[start_cig_id].pos_query;
      seq_start_in_window = (aln.get_cigar()[start_cig_id].pos_ref + aln.get_pos() - 1) - window_start;
      seq_start_in_window = std::max((uint32_t) 0, (uint32_t) ((int32_t) seq_start_in_window - 0));
    }

    if (end_seq == -2) {
      end_seq = seq->get_data_length() - 1 - aln.GetClippedBasesBack();
      seq_end_in_window = (aln.get_pos() - 1 + aln.GetReferenceLengthFromCigar()) - window_start;
      seq_end_in_window = std::min((uint32_t) (temp_window_end - window_start), (uint32_t) ((int32_t) seq_end_in_window + 0));
      end_cig_id = aln.get_cigar().size() - 1;
    } else if (end_seq < 0) {
      fprintf (stderr, "ERROR: end_seq is < 0 and != -2!\n"); exit(1);
    }

    if (aln.get_cigar()[end_cig_id].op == 'D' || aln.get_cigar()[end_cig_id].op == 'I' || aln.get_cigar()[end_cig_id].op == 'S') {

      for (; end_cig_id >= 0; end_cig_id--) {
        if (aln.get_cigar()[end_cig_id].op == 'M' || aln.get_cigar()[end_cig_id].op == '=' || aln.get_cigar()[end_cig_id].op == 'X') { break; }
      }
      end_seq = aln.get_cigar()[end_cig_id].pos_query + aln.get_cigar()[end_cig_id].count - 1;
      seq_end_in_window = (aln.get_cigar()[end_cig_id].pos_ref + aln.get_cigar()[end_cig_id].count - 1 + aln.get_pos() - 1) - window_start;
      seq_end_in_window = std::max((uint32_t) 0, (uint32_t) ((int32_t) seq_end_in_window - 0));
    }

    std::string seq_data = GetSubstring((char *) (seq->get_data() + start_seq), end_seq - start_seq + 1);
    std::string seq_qual = (seq->get_quality() != NULL) ? (GetSubstring((char *) (seq->get_quality() + start_seq), end_seq - start_seq + 1)) : (std::string((end_seq - start_seq + 1), '!' + 0));

    // Safety percaution.
    if (seq_data.size() < 2) { continue; }

    double avg_qual = 0.0f;
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

      starts_on_read.push_back(start_seq);
      ends_on_read.push_back(end_seq);
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

void CreateConsensusAln(const ProgramParameters &parameters, int32_t num_window_threads, const SingleSequence *contig, const std::vector<SingleSequence *> &ctg_alns, std::map<const SingleSequence *, int64_t> &aln_lens_on_ref, std::string &ret_consensus, FILE *fp_out_cons) {
  std::stringstream ss_cons;

  if (parameters.do_erc == false && fp_out_cons) {
    fprintf (fp_out_cons, ">Consensus_%s\n", contig->get_header());
    fflush (fp_out_cons);
  }

  int64_t num_windows = ceil((float) contig->get_sequence_length() / (float) parameters.window_len);
  if (parameters.do_erc == false) {
    LOG_DEBUG ("current_contig->get_sequence_length() = %ld, parameters.window_len = %ld, num_windows = %ld\n", contig->get_sequence_length(), parameters.window_len, num_windows);
  }

  TicToc clk1;
  clk1.start();
  // Build the interval tree for fast overlap calculation.
  std::vector<IntervalSS> aln_intervals;
  for (int64_t i=0; i<ctg_alns.size(); i++) {
    int64_t aln_start = ctg_alns[i]->get_aln().get_pos() - 1;
    int64_t aln_end = aln_start + ctg_alns[i]->get_aln().GetReferenceLengthFromCigar() - 1;
    aln_intervals.push_back(IntervalSS(aln_start, aln_end, ctg_alns[i]));
  }
  IntervalTreeSS aln_interval_tree(aln_intervals);
  clk1.stop();
  LOG_DEBUG("CPU time for building the IntervalTree: %f sec.\n", clk1.get_secs());

  double clk_total_spoa = 0.0;
  double clk_total_extract = 0.0;

  // Process the genome in windows, but also process windows in batches. Each batch is processed in multiple threads,
  // then the results are collected and output to file. After that, a new batch is loaded.
  for (int64_t window_batch_start = parameters.start_window, num_batches = 0; window_batch_start < num_windows && (parameters.num_batches < 0 || num_batches < parameters.num_batches); window_batch_start += parameters.batch_of_windows, num_batches++) {
    std::vector<std::string> consensus_windows;
    consensus_windows.resize(parameters.batch_of_windows);
    int64_t windows_to_process = std::min(parameters.batch_of_windows, num_windows - window_batch_start);

//    #pragma omp parallel for num_threads(parameters.num_threads) schedule(dynamic, 1)
    #pragma omp parallel for num_threads(num_window_threads) reduction(+:clk_total_extract, clk_total_spoa) schedule(dynamic, 1)
    for (int64_t id_in_batch = 0; id_in_batch < windows_to_process; id_in_batch += 1) {

      TicToc clk2;
      clk2.start();
       int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
       int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin) - 1;
       int32_t thread_id = omp_get_thread_num();

       if (parameters.do_erc == false && thread_id == 0) {
         LOG_MEDHIGH("\r(thread_id = %d) Processing window: %ld bp to %ld bp (%.2f%%)",
                     thread_id, window_start, window_end, 100.0 * ((float) window_start / (float) contig->get_data_length()));
         LOG_DEBUG(", clk_total_extract = %f, clk_total_spoa = %f", clk_total_extract, clk_total_spoa);
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
       clk2.stop();
//       LOG_NEWLINE;
//       LOG_ALL("CPU time for extracting a window: %f sec.\n", clk2.get_secs());
       clk_total_extract += clk2.get_secs();

       TicToc clk3;
       clk3.start();
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

       clk3.stop();
//       LOG_ALL("CPU time for consensus on a window: %f sec.\n", clk3.get_secs());
       clk_total_spoa += clk3.get_secs();

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
//      LOG_MEDHIGH_NOHEADER("\n");
      LOG_MEDHIGH("Batch checkpoint: Performed consensus on all windows, joining the windows now.\n");
    }

    // If window overlapping is not used, use this section instead.
//    for (int64_t id_in_batch = 0; id_in_batch < parameters.batch_of_windows && id_in_batch < num_windows; id_in_batch += 1) {
//      int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
//      int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin);
//      ss_cons << consensus_windows[id_in_batch];
//      if (parameters.do_erc == false && fp_out_cons) {
//        fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
//        fflush(fp_out_cons);
//      }
//    }

    // Deprecated, used for window overlapping.
//     LOG_MEDHIGH_NOHEADER("\n");
//     LOG_MEDHIGH("Batch checkpoint: Performed consensus on all windows, joining the windows now.\n");
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
           // Extract last 30% of sequence from the last consensus window.
           std::string trimmed_window = consensus_windows[id_in_batch-1].substr(std::max((1.0 - parameters.win_ovl_margin * 1.3), 0.0) * consensus_windows[id_in_batch-1].size());

           // Create a vector which will be fed to SPOA.
           std::vector<std::string> windows_for_alignment = {trimmed_window, consensus_windows[id_in_batch]};
           std::vector<std::string> msa;

           // Do the alignment.
           SPOA::generate_msa(msa, windows_for_alignment, SPOA::AlignmentParams(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_ext, SPOA::AlignmentType::kOV), false);
           // SPOA MSA will contain 2 lines and look something like this (the first line is the trimmed window, and the last line is the current window):
           // GCTGGAGTGAGTGGGGAAGAGCGCCACGGACAGTATGTCGCAG-GT---AAA--TGCAGCCACGC--A--TTTGAGCGGTGCAT------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           // --------------------------CGGACAGTATGTCGCAGTGTCACAAACTTGCAGCCACGCAGACCTTTGA-TGGT--ATTGCACAGAATATGGCGGGCGATGCTGACCGGCAGTGAGCAGAACTGGCGCAGCTTCACCCGTTTCGTGCTGTCATGATGACAGAAATTCTGCTTAAGCAGGCAATGGTGGGGATTGTCGGGAGTATCGGCAGCGCCATTGGCGGGCTGTTGGTGGCGGCATCGCGTCAGGCGGTACAGCCATTCAGGCCGCTGCGGCGAAATTTCCATTTTAAACCGGAGGATTTACGGGAACCGGCGGCAAATATGAGCTACAGCAGGGGATTGTTCACGTGGTGAGTTTGTCTTCCACGAAGGAGGCAACCAGCCGGATTGGCGTGGGAATTGTTACGGCTGATGCGCGGCTGCCACCGGCGGTTATGTCGGTCACACCGGGGCAGCATGGCAGACAGCCGGTTCGCAGGCGTCGGGACGTTTGAGCAGAATAACCATGTGGTGATTAACAACGACGGCACGAACGGCAGATAGTCCGGCTGCTGCTGGAAGGCGGCGTATGACAT

//           printf ("\n");
//           for (int32_t i2=0; i2<msa.size(); i2++) {
//             printf ("%s\n", msa[i2].c_str());
//           }
//           fflush(stdout);

           // SPOA returned a padded MSA. This part finds the alignment start position.
           std::stringstream ss_clipped_window;
           int32_t clip_pos = 0;
           int32_t trimmed_window_id = 0, curr_window_id = 1;
           // Find the last non '-' base on the first line of MSA.
           for (clip_pos=(msa[trimmed_window_id].size()-1); clip_pos>=0 && msa[trimmed_window_id][clip_pos] == '-'; clip_pos--);
           // Extract all non '-' bases from the second line. This part will naturally concatenate the consensus.
           // Although there should be no more '-' bases, you never know. Maybe the alignment changes at some point.
           for (clip_pos++; clip_pos<msa[curr_window_id].size(); clip_pos++) {
             if (msa[curr_window_id][clip_pos] != '-') { ss_clipped_window << msa[curr_window_id][clip_pos]; }
           }
           // This is our clipped current window. The part which overlaps the previous window has been removed.
           std::string clipped_window = ss_clipped_window.str();

           // Sanity check if everything is alright. The window should be larger than zero.
           if (clipped_window.size() > 0) { // All went fine.
             ss_cons << clipped_window;
             if (fp_out_cons) {
               fprintf (fp_out_cons, "%s", clipped_window.c_str());
               fflush(fp_out_cons);
             }

           } else { // Something strange happened. Maybe the region is a giant repeat?
             // Just output the window as if there were no overlapping windows.
             ss_cons << consensus_windows[id_in_batch];
             if (fp_out_cons) {
               fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
               fflush(fp_out_cons);
             }

             // Report the problem.
             LOG_DEBUG_NOHEADER("\n");
             LOG_DEBUG("[bad] window_start = %ld, window_end = %ld, clipped_window.size() = %ld\n", window_start, window_end, clipped_window.size());
             LOG_DEBUG_NOHEADER("\n");
             for (int32_t i2=0; i2<windows_for_alignment.size(); i2++) {
               LOG_DEBUG_NOHEADER("%s\n\n", windows_for_alignment[i2].c_str());
             }
             LOG_DEBUG_NOHEADER("\n");
             LOG_DEBUG("Alignment:\n\n");
             fflush(stdout);
             for (int32_t i2=0; i2<msa.size(); i2++) {
               LOG_DEBUG_NOHEADER("%s\n\n", msa[i2].c_str());
             }
             LOG_DEBUG_NOHEADER("\n");
           }
         }
       }
     }

//     LOG_MEDHIGH_NOHEADER("\n");
    if (parameters.do_erc == false) {
      LOG_MEDHIGH("Batch checkpoint: Processed %ld windows and exported the consensus.", parameters.batch_of_windows);
    }
  }

  LOG_DEBUG("CPU total time for creating an IntervalTree: %f\n", clk1.get_secs());
  LOG_DEBUG("CPU total time for extracting windows: %f\n", clk_total_extract);
  LOG_DEBUG("CPU total time for SPOA: %f\n", clk_total_spoa);

  ret_consensus = ss_cons.str();
  if (parameters.do_erc == false && fp_out_cons) {
    fprintf (fp_out_cons, "\n");
  }
}
