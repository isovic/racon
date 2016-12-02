/*
 * consensus_sampling.cc
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
#include "libs/edlibcigar.h"
#include "tictoc.h"

// Used for consensus when overlaps are input to the main program (MHAP, PAF) instead of alignments (SAM).
int ConsensusFromOverlapsWSampling(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &reads,
                                const std::map<std::string, int64_t> &qname_to_ids, const std::vector<OverlapLine> &sorted_overlaps) {
  LOG_MEDHIGH("Running consensus.\n");

  int32_t num_read_threads = (parameters.do_erc) ? (parameters.num_threads) : 1;
  int32_t num_window_threads = (!parameters.do_erc) ? (parameters.num_threads) : 1;

  // For a given contig name (qname), the value is a range of indexes in the sorted overlaps vector.
  std::map<int64_t, ContigOverlapLocation> map_ctg_to_overlaps;

  // Separate overlaps into groups for each contig.
  // Alignments which are called unmapped will be skipped in this step.
  // Also, alignments are filtered by the base quality if available.
  LOG_MEDHIGH("Separating overlaps to individual contigs.\n");
  GroupOverlapsToContigs(sorted_overlaps, map_ctg_to_overlaps);

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
    LOG_MEDHIGH("In total, there are %ld contigs for consensus, each containing:\n", contigs.get_sequences().size());
    for (int32_t i=0; i<contigs.get_sequences().size(); i++) {
      std::string contig_name = contigs.get_sequences()[i]->get_header();
      int64_t contig_id = qname_to_ids.find(contig_name)->second + 1;
      auto it = map_ctg_to_overlaps.find(contig_id);
      if (it == map_ctg_to_overlaps.end()) {
        LOG_MEDHIGH("\t[%ld] %s %ld alignments, contig len: %ld\n", i, contigs.get_sequences()[i]->get_header(), 0, contigs.get_sequences()[i]->get_sequence_length());
      } else {
        auto &ovl_range = map_ctg_to_overlaps[contig_id];
        LOG_MEDHIGH("\t[%ld] %s %ld alignments, contig len: %ld\n", i, contigs.get_sequences()[i]->get_header(), ovl_range.end - ovl_range.start + 1, contigs.get_sequences()[i]->get_sequence_length());
      }
    }
  } else {
    LOG_MEDHIGH("In total, there are %ld sequences for error correction.\n", contigs.get_sequences().size());
  }

  LOG_NEWLINE;

  // Clear the output file for consensus.
  FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "w");
  fclose(fp_out_cons);

  double clk_total_sampling = 0.0;
  double clk_total_cons = 0.0;
  double clk_total_interval_tree = 0.0;

  // For each contig (draft can contain multiple contigs), process alignments which only map to that particular contig.
  #pragma omp parallel for num_threads(num_read_threads) schedule(dynamic, 1)
  for (int32_t i=0; i<contigs.get_sequences().size(); i++) {
    const SingleSequence *contig = contigs.get_sequences()[i];
    std::string contig_name = contig->get_header();
    int64_t contig_id = qname_to_ids.find(contig_name)->second + 1;
    int32_t thread_id = omp_get_thread_num();

    // If we are doing error correction, parallelization is per-read and not per-window.
    // We need to disable some of the debug info.
    if (parameters.do_erc == false) {
      LOG_ALL("Started processing contig %ld / %ld (%.2f%%): %s\n", (i + 1), contigs.get_sequences().size(), 100.0*((float) (i + 1)) / ((float) contigs.get_sequences().size()), contig->get_header());
    }

    auto it = map_ctg_to_overlaps.find(contig_id);

    if (it == map_ctg_to_overlaps.end()) {
      if (parameters.do_erc == false || (parameters.do_erc == true && thread_id == 0)) {
        LOG_MEDHIGH("Contig %ld has 0 overlaps, contig len: %ld, name: '%s'\n", i, contig->get_sequence_length(), contig->get_header());
      }
      continue;
    }

    if (parameters.do_erc == false || (parameters.do_erc == true && thread_id == 0)) {
      LOG_ALL("(thread_id = %d) Sampling overlaps for contig %ld / %ld (%.2f%%): %s\n", thread_id, (i + 1), contigs.get_sequences().size(), 100.0*((float) (i + 1)) / ((float) contigs.get_sequences().size()), contig->get_header());
    }

    TicToc clock_sampling;
    clock_sampling.start();
    // Extract all overlaps for the current contig, and initialize the objects.
    std::vector<OverlapLine> extracted_overlaps(sorted_overlaps.begin()+it->second.start, sorted_overlaps.begin()+it->second.end);
    std::vector<std::shared_ptr<SampledAlignment>> sampled_overlaps;
    if (parameters.do_erc == false) {
      PrepareAndSampleOverlaps(contigs, reads, extracted_overlaps, it->second.start, it->second.end, parameters.num_threads, sampled_overlaps, parameters.window_len, true);
    } else {
      PrepareAndSampleOverlaps(contigs, reads, extracted_overlaps, it->second.start, it->second.end, 1, sampled_overlaps, parameters.window_len, thread_id == 0);
    }
    clock_sampling.stop();
    LOG_DEBUG("CPU time spent on sampling: %.2lf sec.\n", clock_sampling.get_secs());
    clk_total_sampling += clock_sampling.get_secs();

    TicToc clock_interval_tree;
    clock_interval_tree.start();
    // Create an interval tree.
    std::vector<IntervalSampled> intervals;
    for (int64_t i=0; i<sampled_overlaps.size(); i++) {
      if (sampled_overlaps[i]->qpos.size() > 0) {
        intervals.push_back(IntervalSampled(sampled_overlaps[i]->mhap.Bstart, sampled_overlaps[i]->mhap.Bend-1, sampled_overlaps[i]));
      }
    }
    IntervalTreeSampled interval_tree(intervals);
    clock_interval_tree.stop();
    LOG_DEBUG("CPU time spent on building interval tree: %.2lf sec.\n", clock_interval_tree.get_secs());
    clk_total_interval_tree += clock_interval_tree.get_secs();

//    const std::vector<std::shared_ptr<SampledAlignment>> &sampled_overlaps,

//    // This sorts ascending by the pos field.
//    alns.Sort();
//    std::sort(overlaps_final.begin(), overlaps_final.end(), [](const OverlapLine &a, const OverlapLine &b){ return (a.Bid < b.Bid) || (a.Bid == b.Bid && a.Bstart < b.Bstart); } );

    TicToc clock_cons;
    clock_cons.start();
    FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "a");
    std::string consensus;
    if (parameters.do_pileup == false) {
      if (parameters.do_erc == false) {
        CreateConsensusSampling(parameters, num_window_threads, contig, interval_tree, consensus, fp_out_cons);

      } else {
        assert(parameters.do_erc == false && "ERC not implemented in this mode yet!");
//        if (thread_id == 0) {
//          LOG_MEDHIGH("\r(thread_id = %d) Processing contig %ld / %ld (%.2f%%), len: %10ld", thread_id, (i + 1), contigs.get_sequences().size(), 100.0f*(((float) (i)) / ((float) contigs.get_sequences().size())), contig->get_sequence_length());
//        }
//
//        CreateConsensus(parameters, num_window_threads, contig, alns.get_sequences(), aln_lens_on_ref, consensus, NULL);
//        #pragma omp critical
//        {
//          fprintf (fp_out_cons, ">Consensus_%s\n%s\n", contig->get_header(), consensus.c_str());
////          fflush (fp_out_cons);
//        }
      }

    } else {
      assert(parameters.do_pileup == false && "Pileup cannot be implemented in this mode!");

//      Pileup pileup(contig, alns.get_sequences());
////      pileup.Verbose(stdout);
//      pileup.GenerateConsensus(5, consensus);
//      #pragma omp critical
//      fprintf (fp_out_cons, ">Consensus_%s\n%s\n", contig->get_header(), consensus.c_str());
//      #pragma omp critical
//      fflush (fp_out_cons);
    }
    fclose(fp_out_cons);
    clock_cons.stop();
    clk_total_cons += clock_cons.get_secs();

    LOG_DEBUG("CPU time spent on consensus: %.2lf sec.\n", clock_cons.get_secs());
    LOG_DEBUG("Total CPU time spent on this contig: %.2lf sec.\n", (clock_sampling.get_secs() + clock_cons.get_secs()));

    ///////////////////////////////////////
//    LOG_MEDHIGH_NOHEADER("\n");
    if (parameters.do_erc == false) {
      LOG_ALL("Processed %ld bp of %ld bp (100.00%%)\n", contig->get_data_length(), contig->get_data_length());
    }
    LOG_MEDHIGH_NOHEADER("\n");
  }

  return 0;
}

void PerformSampling2(std::shared_ptr<SampledAlignment> &new_sampled_ovl, const SingleSequence* ref, int64_t window_len) {
  auto& mhap = new_sampled_ovl->mhap;

  const char* target = (const char *) (ref->get_data() + mhap.Bstart);
  int32_t target_len = mhap.Bend - mhap.Bstart;

  const char *read = (const char *) new_sampled_ovl->seq->get_data();
  int32_t read_len = read_len;
//    // Reverse the read, because we'll need this anyway.
//    char *r_read = Reverse(read, read_len);

  // Get only parts of the read as the query (only parts covered by the overlap).
  // Forward direction.
  const char *query = (const char *) (read + mhap.Astart);
  int32_t query_len = mhap.Aend - mhap.Astart;

  /////////////////////////
  /// Forward direction ///
  /////////////////////////
  EdlibAlignConfig config_fwd = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE);
  // Offset (fwd) is the distance from the start on ref to the first following window position.
  int32_t fwd_offset = ((int32_t) std::ceil(((double) new_sampled_ovl->mhap.Bstart) / ((double) window_len))) * window_len - new_sampled_ovl->mhap.Bstart - 1;
  int32_t num_rev_to_skip = 0;
  if (fwd_offset < 0) {
    fwd_offset += window_len;
    num_rev_to_skip = 1;
  }
  config_fwd.subscoresOffset = fwd_offset;
  config_fwd.subscoresDistance = window_len;
//  LOG_DEBUG("result_fwd\n");
  EdlibAlignResult result_fwd = edlibAlign((const char *) query, query_len, (const char *) target, target_len, config_fwd);

  /////////////////////////
  /// Reverse direction ///
  /////////////////////////
  EdlibAlignConfig config_rev = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE);
  int32_t rev_offset = new_sampled_ovl->mhap.Bend - ((int32_t) std::floor(((double) new_sampled_ovl->mhap.Bend) / ((double) window_len))) * window_len - 2 + 1;
  int32_t num_fwd_to_skip = 0;
  if (rev_offset < 0) {
    // If we got here, that means that the last base of the target (Bend) is actually on the window border.
    // This means that this is the end of the NW alignment. We can then skip the last fwd score (because it's only one line above this one)
    // and report the last base.
    rev_offset += window_len;
    num_fwd_to_skip = 1;
  }
  config_rev.subscoresOffset = rev_offset;
  config_rev.subscoresDistance = window_len;

  // Reverse direction.
//    const char *r_query = (const char *) (r_read + (mhap.Alen-mhap.Aend));
  char *r_query = Reverse((char *) query, query_len);
  // Reverse the target.
  char *r_target = Reverse((char *) target, target_len);

//  LOG_DEBUG("result_rev\n");
  EdlibAlignResult result_rev = edlibAlign((const char *) r_query, query_len, (const char *) r_target, target_len, config_rev);

  if (r_query) {
    delete[] r_query; r_query = NULL;
  }
  if (r_target) {
    delete[] r_target; r_target = NULL;
  }

//  if (std::string(new_sampled_ovl->seq->get_header()) == "31802056-9f77-4a08-afbc-83dc4f9feda5_Basecall_2D_000_2d LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch408_file87_strand_twodirections:MAP006-1_downloads/pass/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch408_file87_strand.fast5") {
//    printf ("\n");
//    printf ("%s\n", new_sampled_ovl->seq->get_header());
//
//    // Checking the result.
//    printf ("\nconfig_fwd.subscoresOffset = %d, config_fwd.subscoresDistance = %d\n", config_fwd.subscoresOffset, config_fwd.subscoresDistance);
//    printf ("sampling_ovl->mhap.Bstart = %d, sampling_ovl->mhap.Bend = %d, query_length = %d, target_length = %d\n", new_sampled_ovl->mhap.Bstart, new_sampled_ovl->mhap.Bend, query_len, target_len);
//    for (int32_t i=0; i<result_fwd.numSubscores; i++) {
//      printf ("[col = %d] startIdx = %d, length = %d, columnIdx = %d, abs_columnIdx = %d\n",
//              i, result_fwd.subscores[i].startIdx, result_fwd.subscores[i].length, result_fwd.subscores[i].columnIdx, result_fwd.subscores[i].columnIdx + new_sampled_ovl->mhap.Bstart);
//    }
//    printf ("\nconfig_rev.subscoresOffset = %d, config_rev.subscoresDistance = %d\n", config_rev.subscoresOffset, config_rev.subscoresDistance);
//    for (int32_t i=0; i<result_rev.numSubscores; i++) {
//      printf ("[col = %d] startIdx = %d, length = %d, columnIdx = %d, columnIdx_fwd = %d, abs_columnIdx_fwd = %d\n", i, result_rev.subscores[i].startIdx,
//              result_rev.subscores[i].length, result_rev.subscores[i].columnIdx, target_len - result_rev.subscores[i].columnIdx - 1,
//              target_len - result_rev.subscores[i].columnIdx - 1 + new_sampled_ovl->mhap.Bstart);
//    }
//    printf ("\n");
//  }

  ////////////////////////
  /// The funky magic. ///
  ////////////////////////
  if ((result_fwd.numSubscores + num_rev_to_skip) != (result_rev.numSubscores + num_fwd_to_skip)) { //TODO: Ovdje jos ide + num_rev_to_skip
    LOG_NEWLINE;
    LOG_DEBUG("Something went wrong here!");

    printf ("\n");
    printf ("%s\n", new_sampled_ovl->seq->get_header());

    // Checking the result.
    printf ("\nconfig_fwd.subscoresOffset = %d, config_fwd.subscoresDistance = %d\n", config_fwd.subscoresOffset, config_fwd.subscoresDistance);
    printf ("sampling_ovl->mhap.Bstart = %d, sampling_ovl->mhap.Bend = %d, query_length = %d, target_length = %d\n", new_sampled_ovl->mhap.Bstart, new_sampled_ovl->mhap.Bend, query_len, target_len);
    for (int32_t i=0; i<result_fwd.numSubscores; i++) {
      printf ("[col = %d] startIdx = %d, length = %d, columnIdx = %d, abs_columnIdx = %d\n",
              i, result_fwd.subscores[i].startIdx, result_fwd.subscores[i].length, result_fwd.subscores[i].columnIdx, result_fwd.subscores[i].columnIdx + new_sampled_ovl->mhap.Bstart);
    }
    printf ("\nconfig_rev.subscoresOffset = %d, config_rev.subscoresDistance = %d\n", config_rev.subscoresOffset, config_rev.subscoresDistance);
    for (int32_t i=0; i<result_rev.numSubscores; i++) {
      printf ("[col = %d] startIdx = %d, length = %d, columnIdx = %d, columnIdx_fwd = %d, abs_columnIdx_fwd = %d\n", i, result_rev.subscores[i].startIdx,
              result_rev.subscores[i].length, result_rev.subscores[i].columnIdx, target_len - result_rev.subscores[i].columnIdx - 1,
              target_len - result_rev.subscores[i].columnIdx - 1 + new_sampled_ovl->mhap.Bstart);
    }
    printf ("\n");
    printf ("Query:\n%s\n\n", GetSubstring((char *) query, query_len).c_str());
    printf ("Target:\n%s\n\n", GetSubstring((char *) target, target_len).c_str());
    printf ("Query length: %d\n", query_len);
    printf ("Target length: %d\n", target_len);
    printf ("Edit distance: %d\n", result_fwd.editDistance);
    printf ("Alignment length: %d\n", result_fwd.alignmentLength);
    printf ("Num subscores: %d\n", result_fwd.numSubscores);
    fflush(stdout);
//    exit(1);

    // TODO: 100% CPU-a se trosi na prazan window!? I ne ide dalje s njega.

    //  } else if (result_fwd.numSubscores == (result_rev.numSubscores + 1)) {
    //    assert(result_fwd.numSubscores == (result_rev.numSubscores));
    //  } else {
    //    assert(result_fwd.numSubscores == (result_rev.numSubscores));
        // The special snowflake. This is the leftover boundary coordinate.

  } else {

//  if (result_fwd.numSubscores == (result_rev.numSubscores + num_fwd_to_skip)) {
    for (int32_t i=0; i<(result_fwd.numSubscores - num_fwd_to_skip); i++) {

      int32_t sub_fwd_id = i;
      int32_t sub_rev_id = result_rev.numSubscores-i-1 - num_rev_to_skip;

      const int32_t large_int = (std::numeric_limits<int32_t>::max()) / 1000; // Compensate for overflow (allow summation of at least several such values).
      std::vector<int32_t> u(query_len, large_int);
      for (int32_t j=0; j<result_fwd.subscores[sub_fwd_id].length; j++) {
        u[result_fwd.subscores[sub_fwd_id].startIdx + j] = result_fwd.subscores[sub_fwd_id].scores[j];
      }

      std::vector<int32_t> d(query_len, large_int);
      for (int32_t j=0; j<result_rev.subscores[sub_rev_id].length; j++) {
        d[(query_len - result_rev.subscores[sub_rev_id].startIdx) - j - 1] = result_rev.subscores[sub_rev_id].scores[j];
      }

      int32_t min_u_id = 0, max_d_id = 0, min_sum = u[0] + d[0];
      for (int32_t j=0; j<(query_len-1); j++) {
        int32_t sum_vert = (u[j] + d[j]);
        if (sum_vert < min_sum) { // Edit distance alignment - find the minimum.
          min_sum = sum_vert; min_u_id = max_d_id = j;
        }

        int32_t sum_diag = (u[j] + d[j+1]);
        if (sum_diag < min_sum) { // Edit distance alignment - find the minimum.
          min_sum = sum_diag; min_u_id = j; max_d_id = j + 1;
        }
      }
      int32_t sum_vert = (u[query_len-1] + d[query_len-1]);
      if (sum_vert < min_sum) { // Edit distance alignment - find the minimum.
        min_sum = sum_vert; min_u_id = max_d_id = query_len-1;
      }

      int32_t u_target_id = result_fwd.subscores[i].columnIdx + new_sampled_ovl->mhap.Bstart;
      int32_t d_target_id = u_target_id + 1;

      new_sampled_ovl->qpos[u_target_id] = min_u_id + new_sampled_ovl->mhap.Astart;
      new_sampled_ovl->qpos[d_target_id] = max_d_id + new_sampled_ovl->mhap.Astart;

    }

    if (num_fwd_to_skip == 1) {
      int32_t u_target_id = result_fwd.subscores[result_fwd.numSubscores-1].columnIdx + new_sampled_ovl->mhap.Bstart;
      new_sampled_ovl->qpos[u_target_id] = query_len-1 + new_sampled_ovl->mhap.Astart;
    }

    if (num_rev_to_skip == 1) {
      int32_t d_target_id = target_len - result_rev.subscores[result_rev.numSubscores-1].columnIdx - 1 + new_sampled_ovl->mhap.Bstart;
      new_sampled_ovl->qpos[d_target_id] = 0 + new_sampled_ovl->mhap.Astart;
    }
  }



  edlibFreeAlignResult(result_fwd);
  edlibFreeAlignResult(result_rev);
}

void PrepareAndSampleOverlaps(const SequenceFile &refs, const SequenceFile &reads,
                             const std::vector<OverlapLine> &sorted_overlaps, int64_t start_overlap, int64_t end_overlap, int32_t num_threads,
                             std::vector<std::shared_ptr<SampledAlignment>> &extracted_overlaps, int32_t window_len, bool verbose_debug) {

  extracted_overlaps.clear();
  extracted_overlaps.reserve(end_overlap - start_overlap);

  int64_t num_skipped_overlaps = 0;

  #pragma omp parallel for num_threads(num_threads) shared(extracted_overlaps) schedule(dynamic, 1)
  for (int64_t i=0; i<sorted_overlaps.size(); i++) {
    int32_t thread_id = omp_get_thread_num();

    // Just verbose.
    if (verbose_debug == true && thread_id == 0) {
      LOG_ALL("\rSampling overlap: %ld / %ld (%.2f\%), skipped %ld / %ld (%.2f\%)", i, sorted_overlaps.size(), 100.0f*((float) (i + 1)) / ((float) sorted_overlaps.size()), num_skipped_overlaps, sorted_overlaps.size(), 100.0f*((float) (num_skipped_overlaps)) / ((float) sorted_overlaps.size()));
      fflush(stderr);
    }

    // Fetch the relevant data.
    auto mhap = sorted_overlaps[i];

    // Get the ref sequence.
    const SingleSequence* ref = refs.get_sequences()[mhap.Bid - 1];
    // Get the read sequence.
    const SingleSequence* ss_read = reads.get_sequences()[mhap.Aid - 1];


    auto new_sampled_ovl = std::make_shared<SampledAlignment>(ss_read, mhap);

    PerformSampling2(new_sampled_ovl, ref, window_len);
//    PerformDummySampling(new_sampled_ovl, ref, window_len);






//    std::shared_ptr<SingleSequence> seq(new SingleSequence);
//    if (mhap.Brev) {
//    }
//
//    if (!mhap.Brev) {
//      seq->InitAllFromAscii((char *) ss_read->get_header(), ss_read->get_header_length(),
//                            (int8_t *) read, (int8_t *) ss_read->get_quality(), read->get_data_length(),
//                            extracted_overlaps.size(), extracted_overlaps.size());
//    }
//
//
//    std::string ref_name(ref->get_header(), ref->get_header_length());
//
//    char *r_read = Reverse((char *) read->get_data(), read->get_data_length());
//    char *r_ref = Reverse((char *) (ref->get_data()), ref->get_data_length());
//
//    // Forward direction.
//    const char *query = (const char *) (sampling_ovl->seq->get_data() + sampling_ovl->mhap.Astart);
//    int32_t query_length = sampling_ovl->mhap.Aend - sampling_ovl->mhap.Astart;
//    const char *target = (const char *) (ref->get_data() + sampling_ovl->mhap.Bstart);
//    int32_t target_length = sampling_ovl->mhap.Bend - sampling_ovl->mhap.Bstart;
//    // Make a config for the alignment.
//    EdlibAlignConfig config_fwd = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE);
//    // Offset (fwd) is the distance from the start on ref to the first following window position.
//    config_fwd.subscoresOffset = ((int32_t) std::ceil(((double) sampling_ovl->mhap.Bstart) / ((double) window_len))) * window_len - sampling_ovl->mhap.Bstart - 1;
//    config_fwd.subscoresDistance = window_len;
//    EdlibAlignResult result_fwd = edlibAlign((const char *) query, query_length,
//                                         (const char *) target, target_length,
//                                         config_fwd);



//    // Copy the sequence. It might have to be rev-complemented.
//    std::shared_ptr<SingleSequence> seq(new SingleSequence);
//    seq->InitAllFromAscii((char *) read->get_header(), read->get_header_length(),
//                          (int8_t *) read->get_data(), (int8_t *) read->get_quality(), read->get_data_length(),
//                          extracted_overlaps.size(), extracted_overlaps.size());
//
//    // Initialize the sampled alignment object with the new sequence.
//    std::shared_ptr<SampledAlignment> new_ext_ovl(new SampledAlignment(seq, mhap));
//
//    // Check if reversal is needed and actually perform it.
//    if (new_ext_ovl->mhap.Brev) {
//      seq->ReverseComplement();
//      new_ext_ovl->mhap.Astart = mhap.Alen - mhap.Aend;
//      new_ext_ovl->mhap.Aend = mhap.Alen - mhap.Astart - 1;
//    }
//
//    // Sample each overlap.
//    PerformSampling(new_ext_ovl, ref,  window_len);
////    PerformDummySampling(new_ext_ovl, ref, window_len);

    #pragma omp critical
    {
      extracted_overlaps.push_back(new_sampled_ovl);
    }
  }

  LOG_NEWLINE;
}

char *Reverse(const char *s, int32_t len) {
  if (len <= 0) { return NULL; }
  char *rs = new char[(len+1)];
  for (int32_t i=0; i<len; i++) { rs[i] = s[len-i-1]; }
  rs[len] = '\0';
  return rs;
}

void PerformSampling(std::shared_ptr<SampledAlignment> sampling_ovl, const SingleSequence* ref, int64_t window_len) {

  // Forward direction.
  const char *query = (const char *) (sampling_ovl->seq->get_data() + sampling_ovl->mhap.Astart);
  int32_t query_length = sampling_ovl->mhap.Aend - sampling_ovl->mhap.Astart;

  const char *target = (const char *) (ref->get_data() + sampling_ovl->mhap.Bstart);
  int32_t target_length = sampling_ovl->mhap.Bend - sampling_ovl->mhap.Bstart;

  //  config_fwd.subscoresOffset = sampling_ovl->mhap.Bstart;

  EdlibAlignConfig config_fwd = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE);
  // Offset (fwd) is the distance from the start on ref to the first following window position.
  config_fwd.subscoresOffset = ((int32_t) std::ceil(((double) sampling_ovl->mhap.Bstart) / ((double) window_len))) * window_len - sampling_ovl->mhap.Bstart - 1;
  config_fwd.subscoresDistance = window_len;
  EdlibAlignResult result_fwd = edlibAlign((const char *) query, query_length,
                                       (const char *) target, target_length,
                                       config_fwd);

  // Reverse direction.
  char *r_query = Reverse((char *) query, query_length);
  char *r_target = Reverse((char *) (target), target_length);

  EdlibAlignConfig config_rev = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE);
//  int rev_strt = ((ref->get_sequence_length() - strt) % ffst) - 2 > 0 ? ((targetLength - strt) % ffst) - 2 : ((targetLength - strt) % ffst) - 2 + ffst;
//  int rev_offset = ((target_length - sampling_ovl->mhap.Bstart) % window_len) - 2 > 0 ?
//                  ((target_length - sampling_ovl->mhap.Bstart) % window_len) - 2 : ((target_length - sampling_ovl->mhap.Bstart) % window_len) - 2 + window_len;
  // Offset (rev) is the distance from the end on ref (Bend) to the first window before Bend.
  int32_t rev_offset = sampling_ovl->mhap.Bend - ((int32_t) std::floor(((double) sampling_ovl->mhap.Bend) / ((double) window_len))) * window_len - 2 + 1;
  int32_t num_fwd_to_skip = 0;
  if (rev_offset < 0) {
    // If we got here, that means that the last base of the target (Bend) is actually on the window border.
    // This means that this is the end of the NW alignment. We can then skip the last fwd score (because it's only one line above this one)
    // and report the last base.
    rev_offset += window_len;
//    sampling_ovl->qpos[sampling_ovl->mhap.Bend] = sampling_ovl->mhap.Aend;
    num_fwd_to_skip = 1;
  }
  config_rev.subscoresOffset = rev_offset;
  config_rev.subscoresDistance = window_len;

  EdlibAlignResult result_rev = edlibAlign((const char *) r_query, query_length,
                                       (const char *) r_target, target_length,
                                       config_rev);

  if (r_query) { delete[] r_query; r_query = NULL; }
  if (r_target) { delete[] r_target; r_target = NULL; }

//  // Checking the result.
//  printf ("\nconfig_fwd.subscoresOffset = %d, config_fwd.subscoresDistance = %d\n", config_fwd.subscoresOffset, config_fwd.subscoresDistance);
//  printf ("sampling_ovl->mhap.Bstart = %d, sampling_ovl->mhap.Bend = %d, query_length = %d, target_length = %d\n", sampling_ovl->mhap.Bstart, sampling_ovl->mhap.Bend, query_length, target_length);
//  for (int32_t i=0; i<result_fwd.numSubscores; i++) {
//    printf ("[col = %d] startIdx = %d, length = %d, columnIdx = %d, abs_columnIdx = %d\n",
//            i, result_fwd.subscores[i].startIdx, result_fwd.subscores[i].length, result_fwd.subscores[i].columnIdx, result_fwd.subscores[i].columnIdx + sampling_ovl->mhap.Bstart);
//  }
//  printf ("\nconfig_rev.subscoresOffset = %d, config_rev.subscoresDistance = %d\n", config_rev.subscoresOffset, config_rev.subscoresDistance);
//  for (int32_t i=0; i<result_rev.numSubscores; i++) {
//    printf ("[col = %d] startIdx = %d, length = %d, columnIdx = %d, columnIdx_fwd = %d, abs_columnIdx_fwd = %d\n", i, result_rev.subscores[i].startIdx,
//            result_rev.subscores[i].length, result_rev.subscores[i].columnIdx, target_length - result_rev.subscores[i].columnIdx - 1,
//            target_length - result_rev.subscores[i].columnIdx - 1 + sampling_ovl->mhap.Bstart);
//  }
//  printf ("\n");

  // Check if something funky happened.

//  if (num_fwd_to_skip == 1) {
//    // Pair-up the last cell to the one before it.
//    int32_t sub_fwd_id = result_fwd.numSubscores - num_fwd_to_skip - 1;
//    // This could be done faster, but for simplicity just fill everything and then check.
//    std::vector<int32_t> u(query_length, 0);
//    for (int32_t j=0; j<result_fwd.subscores[sub_fwd_id].length; j++) {
//      u[result_fwd.subscores[sub_fwd_id].startIdx + j] = result_fwd.subscores[sub_fwd_id].scores[j];
//    }
//    int32_t u_target_id = result_fwd.subscores[sub_fwd_id].columnIdx + sampling_ovl->mhap.Bstart;
//    if (u[u.size()-1] > u[u.size()-2]) {
//      // Deletion.
//      sampling_ovl->qpos[u_target_id] = u.size() - 1  + sampling_ovl->mhap.Astart;
//    } else {
//      // Match/mismatch.
//      sampling_ovl->qpos[u_target_id] = u.size() - 2  + sampling_ovl->mhap.Astart;
//    }
//  }

  if (result_fwd.numSubscores == (result_rev.numSubscores + num_fwd_to_skip)) {
    for (int32_t i=0; i<(result_fwd.numSubscores - num_fwd_to_skip); i++) {

      int32_t sub_fwd_id = i;
      int32_t sub_rev_id = result_rev.numSubscores-i-1;

      const int32_t large_int = (std::numeric_limits<int32_t>::max()) / 1000; // Compensate for overflow (allow summation of at least several such values).
      std::vector<int32_t> u(query_length, large_int);
      for (int32_t j=0; j<result_fwd.subscores[sub_fwd_id].length; j++) {
        u[result_fwd.subscores[sub_fwd_id].startIdx + j] = result_fwd.subscores[sub_fwd_id].scores[j];
      }

      std::vector<int32_t> d(query_length, large_int);
      for (int32_t j=0; j<result_rev.subscores[sub_rev_id].length; j++) {
        d[(query_length - result_rev.subscores[sub_rev_id].startIdx) - j - 1] = result_rev.subscores[sub_rev_id].scores[j];
      }

      int32_t min_u_id = 0, max_d_id = 0, min_sum = u[0] + d[0];
      for (int32_t j=0; j<(query_length-1); j++) {
        int32_t sum_vert = (u[j] + d[j]);
        if (sum_vert < min_sum) { // Edit distance alignment - find the minimum.
          min_sum = sum_vert; min_u_id = max_d_id = j;
        }

        int32_t sum_diag = (u[j] + d[j+1]);
        if (sum_diag < min_sum) { // Edit distance alignment - find the minimum.
          min_sum = sum_diag; min_u_id = j; max_d_id = j + 1;
        }
      }
      int32_t sum_vert = (u[query_length-1] + d[query_length-1]);
      if (sum_vert < min_sum) { // Edit distance alignment - find the minimum.
        min_sum = sum_vert; min_u_id = max_d_id = query_length-1;
      }

      int32_t u_target_id = result_fwd.subscores[i].columnIdx + sampling_ovl->mhap.Bstart;
      int32_t d_target_id = u_target_id + 1;


      sampling_ovl->qpos[u_target_id] = min_u_id + sampling_ovl->mhap.Astart;
      sampling_ovl->qpos[d_target_id] = max_d_id + sampling_ovl->mhap.Astart;

    }

//  } else if (result_fwd.numSubscores == (result_rev.numSubscores + 1)) {
//    assert(result_fwd.numSubscores == (result_rev.numSubscores));

  } else {
//    assert(result_fwd.numSubscores == (result_rev.numSubscores));

  }

  edlibFreeAlignResult(result_fwd);
  edlibFreeAlignResult(result_rev);

//  for (int32_t i=0; i<result_fwd.numSubscores; i++) {
//
//    for (int32_t j=0; j<target_length; j++) {
//    }
//  }

//    // moras rucno postaviti parametre u config
//    // strt je start, ffst je offset
//
//    std::string query_rev(query, queryLength);
//    std::reverse(query_rev.begin(), query_rev.end());
//    std::string target_rev(target, targetLength);
//    std::reverse(target_rev.begin(), target_rev.end());
//    rev_start = ((target - start) % offset ) - 2
//
//        ovako sam racunao reverzni start
//        int rev_strt = ((targetLength - strt) % ffst) - 2 > 0 ? ((targetLength - strt) % ffst) - 2 : ((targetLength - strt) % ffst) - 2 + ffst;
//
//        strt je obicni start, ffst je offset

}

void PerformDummySampling(std::shared_ptr<SampledAlignment> sampling_ovl, const SingleSequence* ref, int64_t window_len) {
//  EdlibAlignConfig config = edlibNewAlignConfig(k, modeCode, alignTask);
//  config.subscoresOffset = strt;
//  config.subscoresDistance = ffst;
//
//  // moras rucno postaviti parametre u config
//  // strt je start, ffst je offset
//
//  std::string query_rev(query, queryLength);
//  std::reverse(query_rev.begin(), query_rev.end());
//  std::string target_rev(target, targetLength);
//  std::reverse(target_rev.begin(), target_rev.end());
//  rev_start = ((target - start) % offset ) - 2

  std::string ref_name = ref->get_header();

  int64_t aln_start = 0, aln_end = 0, editdist = 0;
  std::vector<unsigned char> alignment;
  int rcaln = EdlibNWWrapper(sampling_ovl->seq->get_data() + sampling_ovl->mhap.Astart, (sampling_ovl->mhap.Aend - sampling_ovl->mhap.Astart),
                          ref->get_data() + sampling_ovl->mhap.Bstart, (sampling_ovl->mhap.Bend - sampling_ovl->mhap.Bstart),
                          &aln_start, &aln_end, &editdist, alignment);

  if (!rcaln) {
    char *cigar_cstring = NULL;
    int rccig = edlibAlignmentToCigar(&alignment[0], alignment.size(), EDLIB_CIGAR_EXTENDED, &cigar_cstring);

    std::string cigar_string(cigar_cstring);

    SequenceAlignment aln;
    aln.SetCigarFromString(cigar_string);
    aln.set_pos(sampling_ovl->mhap.Bstart + aln_start + 1);
    aln.set_flag((int32_t) (16 * sampling_ovl->mhap.Brev));     // 0 if fwd and 16 if rev.
    aln.set_mapq(40);
    aln.set_rname(ref_name);

    // Remove insertions at front.
    for (int32_t j=0; j<aln.cigar().size(); j++) {
      if (aln.cigar()[j].count == 0) { continue; }
      if (aln.cigar()[j].op != 'D') { break; }
      aln.cigar()[j].count = 0;
    }
    for (int32_t j=0; j<aln.cigar().size(); j++) {
      if (aln.cigar()[j].count == 0) { continue; }
      if (aln.cigar()[j].op != 'I') { break; }
      aln.cigar()[j].op = 'S';
    }
    // Remove insertions at back.
    for (int32_t j=(aln.cigar().size()-1); j>=0; j--) {
      if (aln.cigar()[j].count == 0) { continue; }
      if (aln.cigar()[j].op != 'D') { break; }
      aln.cigar()[j].count = 0;
    }
    for (int32_t j=(aln.cigar().size()-1); j>=0; j--) {
      if (aln.cigar()[j].count == 0) { continue; }
      if (aln.cigar()[j].op != 'I') { break; }
      aln.cigar()[j].op = 'S';
    }

    // If the overlap does not cover the entire read (and most likely it does not).
    if (sampling_ovl->mhap.Astart > 0) {
      if (aln.cigar().size() > 0 && aln.cigar().front().op == 'S') { aln.cigar().front().count += sampling_ovl->mhap.Astart; }
      else { CigarOp new_op; new_op.op = 'S'; new_op.count = sampling_ovl->mhap.Astart; aln.cigar().insert(aln.cigar().begin(), new_op); }
    }
    if ((sampling_ovl->mhap.Aend) < (sampling_ovl->mhap.Alen)) {
      if (aln.cigar().size() > 0 && aln.cigar().back().op == 'S') { aln.cigar().back().count += (sampling_ovl->mhap.Alen - sampling_ovl->mhap.Aend); }
      else { CigarOp new_op; new_op.op = 'S'; new_op.count = (sampling_ovl->mhap.Alen - sampling_ovl->mhap.Aend); aln.cigar().insert(aln.cigar().end(), new_op); }
    }

    aln.RecalcCigarPositions();

    int64_t m=0, x=0, eq=0, ins=0, del=0;
    for (int32_t j=0; j<aln.cigar().size(); j++) {
      if (aln.cigar()[j].op == 'M') { m += aln.cigar()[j].count; }
      if (aln.cigar()[j].op == '=') { eq += aln.cigar()[j].count; }
      if (aln.cigar()[j].op == 'X') { x += aln.cigar()[j].count; }
      if (aln.cigar()[j].op == 'I') { ins += aln.cigar()[j].count; }
      if (aln.cigar()[j].op == 'D') { del += aln.cigar()[j].count; }
    }
    aln.optional().push_back(FormatString("X1:Z:equal=%ld_x=%ld_ins=%ld_del=%ld", eq, x, ins, del));

//    sampling_ovl->seq->InitAlignment(aln);
    // This should normally *NEVER* be performed, but since this is a dummy function meant for debugging.
    ((SingleSequence *) sampling_ovl->seq)->InitAlignment(aln);


    std::vector<int32_t> sampling_rpos;
    int32_t start_sp = ((int32_t) std::floor(((double) sampling_ovl->mhap.Bstart) / ((double) window_len))) * 500;
    int32_t end_sp = ((int32_t) std::ceil(((double) sampling_ovl->mhap.Bend) / ((double) window_len))) * 500;
    for (int32_t rp=start_sp; rp<=end_sp; rp+=window_len) {
      if (rp >= sampling_ovl->mhap.Bstart && rp <= sampling_ovl->mhap.Bend) {
//        sampling_rpos.push_back(rp);
        int64_t start_cig_id = 0;
        int32_t qp = aln.FindBasePositionOnRead(rp, &start_cig_id);
        sampling_ovl->qpos[rp] = qp;

        if ((rp-1 >= 0)) {
          int32_t qp1 = aln.FindBasePositionOnRead(rp-1, &start_cig_id);
          sampling_ovl->qpos[rp-1] = qp1;
        }
//        printf ("(%ld, %ld), ", rp, qp);
      }
    }
//    printf ("\n[start_sp=%d, end_sp=%d, Astart=%d, Aend=%d, Bstart=%d, Bend=%d] ", start_sp, end_sp, sampling_ovl->mhap.Astart, sampling_ovl->mhap.Aend, sampling_ovl->mhap.Bstart, sampling_ovl->mhap.Bend);
//    for (auto itm = sampling_ovl->qpos.begin(); itm != sampling_ovl->qpos.end(); itm++) {
//      printf ("(%d, %d), ", itm->first, itm->second);
//    }
//    printf (";\n");
//    fflush(stdout);
  }
}

void CreateConsensusSampling(const ProgramParameters &parameters, int32_t num_window_threads, const SingleSequence *contig,
                             IntervalTreeSampled &interval_tree, std::string &ret_consensus, FILE *fp_out_cons) {
  std::stringstream ss_cons;

  if (parameters.do_erc == false && fp_out_cons) {
    fprintf (fp_out_cons, ">Consensus_%s\n", contig->get_header());
    fflush (fp_out_cons);
  }

  int64_t num_windows = ceil((float) contig->get_sequence_length() / (float) parameters.window_len);
  if (parameters.do_erc == false) {
    LOG_DEBUG ("current_contig->get_sequence_length() = %ld, parameters.window_len = %ld, num_windows = %ld\n", contig->get_sequence_length(), parameters.window_len, num_windows);
  }

  // Process the genome in windows, but also process windows in batches. Each batch is processed in multiple threads,
  // then the results are collected and output to file. After that, a new batch is loaded.
  for (int64_t window_batch_start = parameters.start_window, num_batches = 0; window_batch_start < num_windows && (parameters.num_batches < 0 || num_batches < parameters.num_batches); window_batch_start += parameters.batch_of_windows, num_batches++) {
    std::vector<std::string> consensus_windows;
    consensus_windows.resize(parameters.batch_of_windows);
    int64_t windows_to_process = std::min(parameters.batch_of_windows, num_windows - window_batch_start);

    #pragma omp parallel for num_threads(num_window_threads) schedule(dynamic, 1)
    for (int64_t id_in_batch = 0; id_in_batch < windows_to_process; id_in_batch += 1) {

       int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
       int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin) - 1;
       int32_t thread_id = omp_get_thread_num();

       if (parameters.do_erc == false && thread_id == 0) {
         LOG_ALL("\r(thread_id = %d) Processing window: %ld bp to %ld bp (%.2f%%)", thread_id, window_start, window_end, 100.0 * ((float) window_start / (float) contig->get_data_length()));
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
       ExtractWindowFromSampledSeqs(contig, interval_tree,
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
//
//       if (id_in_batch >= 30) {
//         break;
//       }
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

    if (parameters.do_erc == false) {
      LOG_MEDHIGH("Batch checkpoint: Processed %ld windows and exported the consensus.\n", parameters.batch_of_windows);
    }
  }

  ret_consensus = ss_cons.str();
  if (parameters.do_erc == false && fp_out_cons) {
    fprintf (fp_out_cons, "\n");
  }
}

void ExtractWindowFromSampledSeqs(const SingleSequence *contig, IntervalTreeSampled &interval_tree,
                                  int64_t window_start, int64_t window_end, double qv_threshold, bool use_contig_qvs,
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
  std::vector<IntervalSampled> intervals;
  interval_tree.findOverlapping(window_start, temp_window_end, intervals);

//  printf ("\nIntervaltree returned %ld hits:\n", intervals.size());
//  fflush(stdout);

  // For each seq, extract its segment which falls into the window.
  for (int64_t i=0; i<intervals.size(); i++) {
    std::shared_ptr<SampledAlignment> sal = intervals[i].value;
//    std::shared_ptr<SingleSequence> seq = sal->seq;
    const SingleSequence *seq = sal->seq;

    int32_t start_seq = sal->find(window_start);
    int32_t end_seq = sal->find(temp_window_end);

    auto aln = sal->seq->get_aln();

    uint32_t seq_start_in_window = 0;
    uint32_t seq_end_in_window = temp_window_end - window_start;

    // If the window starts within the interval (after the window_start location).
    if (start_seq == -1 && intervals[i].start > window_start) {
      start_seq = intervals[i].value->mhap.Astart;
      seq_start_in_window = intervals[i].start - window_start;
      seq_start_in_window = std::max((uint32_t) 0, (uint32_t) ((int32_t) seq_start_in_window - 0));
    } else if (start_seq < 0) {

      LOG_DEBUG("start_seq < 0 and it's not a window start! Skipping the read.\n");

      printf ("Problematic read:\n%s\n", intervals[i].value->seq->get_header());
      printf ("All reads in the window:\n");
      for (int32_t j=0; j<intervals.size(); j++) {
        printf ("%s\n", intervals[j].value->seq->get_header());
      }
      fflush(stdout);

      printf ("=============> start_seq = %d, end_seq = %d, Astart = %d, Aend = %d, Alen = %d, Bstart = %d, Bend = %d\n",
              start_seq, end_seq, intervals[i].value->mhap.Astart, intervals[i].value->mhap.Aend, intervals[i].value->mhap.Alen,
              intervals[i].value->mhap.Bstart, intervals[i].value->mhap.Bend);
      for (auto itm = sal->qpos.begin(); itm != sal->qpos.end(); itm++) {
        printf ("(%d, %d), ", itm->first, itm->second);
      }
      printf (";\n");
      fflush(stdout);
      fflush(stdout);

//      assert(start_seq >= 0 && "ERROR: start_seq is < 0");
//      exit(1);
      continue;
    }

    if (end_seq == -1 && intervals[i].stop < window_end) {
      end_seq = intervals[i].value->mhap.Aend;
      seq_end_in_window = intervals[i].stop - window_start;
      seq_end_in_window = std::min((uint32_t) (temp_window_end - window_start), (uint32_t) ((int32_t) seq_end_in_window + 0));
    } else if (end_seq < 0) {
//      assert(end_seq >= 0 && ("ERROR: start_seq is < 0"));
      LOG_DEBUG("end_seq < 0 and it's not a window end! Skipping the read.\n");

      printf ("Problematic read:\n%s\n", intervals[i].value->seq->get_header());
      printf ("All reads in the window:\n");
      for (int32_t j=0; j<intervals.size(); j++) {
        printf ("%s\n", intervals[j].value->seq->get_header());
      }
      fflush(stdout);

      printf ("=============> start_seq = %d, end_seq = %d, Astart = %d, Aend = %d, Alen = %d, Bstart = %d, Bend = %d\n",
              start_seq, end_seq, intervals[i].value->mhap.Astart, intervals[i].value->mhap.Aend, intervals[i].value->mhap.Alen,
              intervals[i].value->mhap.Bstart, intervals[i].value->mhap.Bend);
      for (auto itm = sal->qpos.begin(); itm != sal->qpos.end(); itm++) {
        printf ("(%d, %d), ", itm->first, itm->second);
      }
      printf (";\n");
      fflush(stdout);
      fflush(stdout);

      continue;
    }

////    printf ("%s\tstart_seq = %d, end_seq = %d, Astart = %d, Aend = %d, Alen = %d\n", seq->get_header(), start_seq, end_seq, intervals[i].value->mhap.Astart, intervals[i].value->mhap.Aend, intervals[i].value->mhap.Alen);
//    printf ("=============> start_seq = %d, end_seq = %d, Astart = %d, Aend = %d, Alen = %d, Bstart = %d, Bend = %d\n",
//            start_seq, end_seq, intervals[i].value->mhap.Astart, intervals[i].value->mhap.Aend, intervals[i].value->mhap.Alen,
//            intervals[i].value->mhap.Bstart, intervals[i].value->mhap.Bend);
//    for (auto itm = sal->qpos.begin(); itm != sal->qpos.end(); itm++) {
//      printf ("(%d, %d), ", itm->first, itm->second);
//    }
//    printf (";\n");
//    fflush(stdout);
//    fflush(stdout);
//
//    if ((end_seq - start_seq) > 700) {
//      printf ("\nNesto cudno!!\n---\n");
//      fflush(stdout);
//    }

//    if ((start_seq < 0 && intervals[i].start <= window_start) || (end_seq < 0 && intervals[i].stop >= window_end)) {
//      continue;
//    }

    std::string seq_data = GetSubstring((char *) (seq->get_data() + start_seq), end_seq - start_seq + 1);
    std::string seq_qual = (seq->get_quality() != NULL) ? (GetSubstring((char *) (seq->get_quality() + start_seq), end_seq - start_seq + 1)) : (std::string((end_seq - start_seq + 1), '!' + 0));

    // Safety percaution.
    if (seq_data.size() < 2) { continue; }

    double avg_qual;
    for (int64_t j=0; j<seq_qual.size(); j++) {
      avg_qual += (double) (seq_qual[j] - '!');
    }
    avg_qual /= std::max((double) seq_qual.size(), 1.0);

    if (avg_qual >= qv_threshold) {
      window_refs.push_back(&(*seq));
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

//  printf ("\n");
//  fflush(stdout);
}
