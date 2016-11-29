/*
 * consensus-overlap.cc
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

int GroupOverlapsToContigs(const std::vector<OverlapLine> &sorted_overlaps, std::map<int64_t, ContigOverlapLocation> &map_ctg_to_overlaps) {
  map_ctg_to_overlaps.clear();

  if (sorted_overlaps.size() == 0) { return 1; }

  ContigOverlapLocation ctg_loc;
  ctg_loc.ctg_id = sorted_overlaps[0].Bid;
  for (int64_t i=1; i<sorted_overlaps.size(); i++) {
    if (sorted_overlaps[i].Bid == ctg_loc.ctg_id) {
      ctg_loc.end = i;
    } else {
      map_ctg_to_overlaps[sorted_overlaps[i-1].Bid] = ctg_loc;
      ctg_loc.start = ctg_loc.end = i;
      ctg_loc.ctg_id = sorted_overlaps[i].Bid;
    }
  }

  // Last update and push the last streak to the map.
  map_ctg_to_overlaps[sorted_overlaps.back().Bid] = ctg_loc;

  return 0;
}

// Used for consensus when overlaps are input to the main program (MHAP, PAF) instead of alignments (SAM).
int ConsensusFromOverlaps(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &reads,
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
      LOG_ALL("(thread_id = %d) Aligning overlaps for contig %ld / %ld (%.2f%%): %s\n", thread_id, (i + 1), contigs.get_sequences().size(), 100.0*((float) (i + 1)) / ((float) contigs.get_sequences().size()), contig->get_header());
    }

    TicToc clock_aln;
    clock_aln.start();
    SequenceFile alns;
    std::vector<OverlapLine> extracted_overlaps(sorted_overlaps.begin()+it->second.start, sorted_overlaps.begin()+it->second.end);
    if (parameters.do_erc == false) {
      AlignOverlaps(contigs, reads, extracted_overlaps, parameters.num_threads, alns, true);
    } else {
      AlignOverlaps(contigs, reads, extracted_overlaps, 1, alns, thread_id == 0);
    }
    // Hash all the alignment lengths (which will be used a lot).
    std::map<const SingleSequence *, int64_t> aln_lens_on_ref;
    for (int64_t i=0; i<alns.get_sequences().size(); i++) {
      aln_lens_on_ref[alns.get_sequences()[i]] = alns.get_sequences()[i]->get_aln().GetReferenceLengthFromCigar();
    }
    // This sorts ascending by the pos field.
    alns.Sort();
    clock_aln.stop();
    LOG_ALL("CPU time spent on alignment: %.2lf sec.\n", clock_aln.get_secs());

    TicToc clock_cons;
    clock_cons.start();
    FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "a");
    std::string consensus;
    if (parameters.do_pileup == false) {
      if (parameters.do_erc == false) {
        CreateConsensusAln(parameters, num_window_threads, contig, alns.get_sequences(), aln_lens_on_ref, consensus, fp_out_cons);

      } else {
        if (thread_id == 0) {
          LOG_MEDHIGH("\r(thread_id = %d) Processing contig %ld / %ld (%.2f%%), len: %10ld", thread_id, (i + 1), contigs.get_sequences().size(), 100.0f*(((float) (i)) / ((float) contigs.get_sequences().size())), contig->get_sequence_length());
        }

        CreateConsensusAln(parameters, num_window_threads, contig, alns.get_sequences(), aln_lens_on_ref, consensus, NULL);
        #pragma omp critical
        {
          fprintf (fp_out_cons, ">Consensus_%s\n%s\n", contig->get_header(), consensus.c_str());
//          fflush (fp_out_cons);
        }
      }

    } else {
      Pileup pileup(contig, alns.get_sequences());
//      pileup.Verbose(stdout);
      pileup.GenerateConsensus(5, consensus);
      #pragma omp critical
      fprintf (fp_out_cons, ">Consensus_%s\n%s\n", contig->get_header(), consensus.c_str());
      #pragma omp critical
      fflush (fp_out_cons);
    }
    fclose(fp_out_cons);
    clock_cons.stop();

    LOG_ALL("CPU time spent on consensus: %.2lf sec.\n", clock_cons.get_secs());
    LOG_ALL("Total CPU time spent on this contig: %.2lf sec.\n", (clock_aln.get_secs() + clock_cons.get_secs()));

    ///////////////////////////////////////
//    LOG_MEDHIGH_NOHEADER("\n");
    if (parameters.do_erc == false) {
      LOG_ALL("Processed %ld bp of %ld bp (100.00%%)\n", contig->get_data_length(), contig->get_data_length());
    }
    LOG_MEDHIGH_NOHEADER("\n");
  }

  return 0;
}
