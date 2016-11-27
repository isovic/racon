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
        CreateConsensus(parameters, num_window_threads, contig, ctg_alns, aln_lens_on_ref, consensus, fp_out_cons);

      } else {
        if (thread_id == 0) {
          LOG_MEDHIGH("\r(thread_id = %d) Processing contig %ld / %ld (%.2f%%), len: %10ld", thread_id, (i + 1), ctg_names.size(), 100.0f*(((float) (i)) / ((float) ctg_names.size())), contig->get_sequence_length());
        }

        CreateConsensus(parameters, num_window_threads, contig, ctg_alns, aln_lens_on_ref, consensus, NULL);
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
