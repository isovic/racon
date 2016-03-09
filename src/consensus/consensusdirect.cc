/*
 * consensusdirect.cc
 *
 *  Created on: Feb 29, 2016
 *      Author: isovic
 */

#include "consensus/consensus.h"
#include "../../codebase/spoa/src/poa.hpp"
#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include <stdint.h>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <omp.h>

std::vector<size_t> soort(const std::vector<std::string>& windows) {
    std::vector<size_t> indices(windows.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return windows[a].size() > windows[b].size(); }
    );
    return indices;
}

void ExtractWindowFromAlns(const std::vector<const SingleSequence *> &alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens, int64_t window_start, int64_t window_end, std::vector<std::string> &window_seqs, std::vector<std::string> &window_qv, FILE *fp_window) {
  std::vector<SingleSequence *> candidates;
  for (int64_t i=0; i<alns.size(); i++) {
    auto aln = alns[i]->get_aln();
    int64_t aln_ref_len = aln_ref_lens.find(alns[i])->second;
    if ((aln.pos - 1) > window_end) {
      break;
    } else if ((aln.pos - 1) <= window_start && (aln.pos - 1 + aln_ref_len) >= window_end) {

      int64_t start_cig_id = 0, end_cig_id = 0;
      int64_t start_seq = aln.FindBasePositionOnRead(aln.cigar, window_start, &start_cig_id);
      int64_t end_seq = aln.FindBasePositionOnRead(aln.cigar, window_end, &end_cig_id);

      if (start_seq < 0) { start_seq = 0; }
      if (end_seq < 0) { end_seq = alns[i]->get_data_length() - 1; }

      window_seqs.push_back(GetSubstring((char *) (alns[i]->get_data() + start_seq), end_seq - start_seq + 1));
      if (alns[i]->get_quality() != NULL) {
        window_qv.push_back(GetSubstring((char *) (alns[i]->get_quality() + start_seq), end_seq - start_seq + 1));
      }

      if (fp_window) {
        fprintf (fp_window, ">%s Window_%d_to_%d\n%s\n", alns[i]->get_header(), window_start, window_end, window_seqs.back().c_str());
//        fprintf (fp_window, "@%s Window_%d_to_%d\n%s\n", alns[i]->get_header(), window_start, window_end, window_seqs.back().c_str());
//        fprintf (fp_window, "+\n%s\n", window_qv.back().c_str());
      }

    }
  }
}

int ConsensusDirectFromAln(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &alns) {
  LOG_ALL("Running consensus - directly from alignments.\n");

  std::vector<std::string> ctg_names;
  std::map<std::string, std::vector<const SingleSequence *> > all_ctg_alns;

  // Separate alignments into groups for each contig.
  // Alignments which are called unmapped will be skipped in this step.
  LOG_ALL("Separating alignments to individual contigs.\n");
  AlignmentsToContigs(alns, ctg_names, all_ctg_alns);

  // Verbose.
  LOG_ALL("In total, there are %ld contigs, each containing:\n", ctg_names.size());
  for (int32_t i=0; i<ctg_names.size(); i++) {
    LOG_ALL("\t[%ld] %s %ld alignments\n", i, ctg_names[i].c_str(), all_ctg_alns.find(ctg_names[i])->second.size());
  }

  // Hash the sequences by their name.
  std::map<std::string, const SingleSequence *> rname_to_seq;
  for (int32_t i=0; i<contigs.get_sequences().size(); i++) {
    rname_to_seq[contigs.get_sequences()[i]->get_header()] = contigs.get_sequences()[i];
    rname_to_seq[TrimToFirstSpace(contigs.get_sequences()[i]->get_header())] = contigs.get_sequences()[i];
  }

  // Hash all the alignment lengths (which will be used a lot).
  std::map<const SingleSequence *, int64_t> aln_ref_lens;
  for (int64_t i=0; i<alns.get_sequences().size(); i++) {
    aln_ref_lens[alns.get_sequences()[i]] = alns.get_sequences()[i]->get_aln().GetReferenceLengthFromCigar();
  }

  // Debug output of alternate contigs, aligned to the raw contig (input sequence), in SAM format.
  FILE *fp_alt_contig_path = NULL;
  if (parameters.temp_alt_contig_path != "") {
    fp_alt_contig_path = fopen(parameters.temp_alt_contig_path.c_str(), "w");
    fprintf (fp_alt_contig_path, "@HD\tVN:1.0\tSO:unknown\n");
    for (int32_t i=0; i<ctg_names.size(); i++) {
      fprintf (fp_alt_contig_path, "@SQ\tSN:%s\tLN:%ld\n", ctg_names[i].c_str(), rname_to_seq[ctg_names[i]]->get_data_length());
    }
    fprintf (fp_alt_contig_path, "@PG\tID:consise PN:consise\n");
  }

  // Clear the output file for consensus.
  FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "w");
  fclose(fp_out_cons);

  // For each contig (draft can contain multiple contigs), process alignments to obtain alternate contigs.
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
    fprintf (fp_out_cons, ">Consensus_%d %s\n", i, ctg_names[i].c_str());

    int64_t num_windows = ceil((float) contig->get_sequence_length() / (float) parameters.window_len);
    LOG_DEBUG ("current_contig->get_sequence_length() = %ld, parameters.window_len = %ld, num_windows = %ld\n", contig->get_sequence_length(), parameters.window_len, num_windows);

    // Process the genome in windows, but also process windows in batches. Each batch is processed in multiple threads,
    // then the results are collected and output to file. After that, a new batch is loaded.
    for (int64_t window_batch_start = 0; window_batch_start < num_windows; window_batch_start += parameters.batch_of_windows) {
      std::vector<std::string> consensus_windows;
      consensus_windows.resize(parameters.batch_of_windows);
      int64_t windows_to_process = std::min(parameters.batch_of_windows, num_windows - window_batch_start);

      #pragma omp parallel for num_threads(parameters.num_threads)
       for (int64_t id_in_batch = 0; id_in_batch < windows_to_process; id_in_batch += 1) {
         int64_t window_start = std::max((int64_t) 0, (int64_t) ((window_batch_start + id_in_batch) * parameters.window_len - (parameters.window_len * parameters.win_ovl_margin)));
         int64_t window_end = window_start + parameters.window_len + (parameters.window_len * parameters.win_ovl_margin);
         int32_t thread_id = omp_get_thread_num();

         if (thread_id == 0) { LOG_ALL("\r(thread_id = %d) Processing window: %ld bp to %ld bp (%.2f%%)", thread_id, window_start, window_end, 100.0 * ((float) window_start / (float) contig->get_data_length())); }

         // Cut a window out of all aligned sequences. This will be fed to an MSA algorithm.
         std::vector<std::string> windows_for_msa;
         std::vector<std::string> quals_for_msa;
//         ExtractWindow(alt_contig_seqs, window_start, window_end, windows_for_msa, fp_window);

         // Chosing the MSA algorithm, and running the consensus on the window.
         if (parameters.msa == "poa") {
           ExtractWindowFromAlns(ctg_alns, aln_ref_lens, window_start, window_end, windows_for_msa, quals_for_msa, NULL);
           if (windows_for_msa.size() == 0) {
               consensus_windows[id_in_batch] = "";
           } else {
               auto indices = soort(windows_for_msa);

               GraphSharedPtr graph = createGraph(windows_for_msa[indices[0]], quals_for_msa[indices[0]]);
               //GraphSharedPtr graph = createGraph(windows_for_msa[0], quals_for_msa[0]);
               graph->topological_sort();
               for (uint32_t w = 1; w < windows_for_msa.size(); ++w) {
                   //auto alignment = createAlignment(windows_for_msa[w], graph,
                   auto alignment = createAlignment(windows_for_msa[indices[w]], graph,
                       AlignmentParams(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_ext, (AlignmentType) parameters.aln_type));
                   alignment->align_sequence_to_graph();
                   alignment->backtrack();
                   graph->add_alignment(alignment->alignment_node_ids(),
                       alignment->alignment_seq_ids(), windows_for_msa[indices[w]], quals_for_msa[indices[w]]);
                       //alignment->alignment_seq_ids(), windows_for_msa[w], quals_for_msa[w]);
               }
               consensus_windows[id_in_batch] = graph->generate_consensus();
           }
         } else {
           FILE *fp_window = NULL;
           std::string window_path = FormatString("%s.%ld", parameters.temp_window_path.c_str(), thread_id);
           if (parameters.temp_window_path != "") {
             fp_window = fopen(window_path.c_str(), "w");
           }
           if (fp_window == NULL) {
             ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Window file not opened!\n");
           }
           ExtractWindowFromAlns(ctg_alns, aln_ref_lens, window_start, window_end, windows_for_msa, quals_for_msa, fp_window);
           fclose(fp_window);
           RunMSAFromSystemLocal(parameters, window_path, consensus_windows[id_in_batch]);
         }
       }

       // This worked for non-overlapping windows.
//       for (int64_t start_window_id = window_batch_start; start_window_id < (window_batch_start + parameters.batch_of_windows) && start_window_id < num_windows; start_window_id += 1) {
//         int64_t id_in_batch = start_window_id - window_batch_start;
//         fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
//         fflush(fp_out_cons);
//         printf ("\nTu sam 1!\n");
//       }
       printf ("\n");

       for (int64_t id_in_batch = 0; id_in_batch < parameters.batch_of_windows && id_in_batch < num_windows; id_in_batch += 1) {
         if (id_in_batch == 0) {
           fprintf (fp_out_cons, "%s", consensus_windows[id_in_batch].c_str());
           fflush(fp_out_cons);
         } else {
//           fprintf (fp_out_cons, "\n>Window2\n%s", consensus_windows[id_in_batch].c_str());
//           fflush(fp_out_cons);
           std::string trimmed_window = consensus_windows[id_in_batch-1].substr((1.0 - parameters.win_ovl_margin * 2) * consensus_windows[id_in_batch-1].size());
           GraphSharedPtr graph = createGraph(trimmed_window);
           graph->topological_sort();
           printf ("%s\n\n", trimmed_window.c_str());

           auto alignment = createAlignment(trimmed_window, graph,
               AlignmentParams(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_ext, AlignmentType::kOV));
           alignment->align_sequence_to_graph();
           alignment->backtrack();
           graph->add_alignment(alignment->alignment_node_ids(), alignment->alignment_seq_ids(), consensus_windows[id_in_batch]);
//           consensus_windows[id_in_batch] = graph->generate_consensus();
           std::vector<std::string> msa;
           graph->generate_msa(msa);
           printf ("msa.size() = %ld\n", msa.size());
           printf ("Sequence %d:\n%s\n", 0, msa[0].c_str());
           for (int32_t i1=0; i1<(msa[0].size()-1); i1++) {
             if (msa[0][i1] == '-' || msa[1][i1] == '-') { printf (" "); }
             else if (msa[0][i1] == msa[1][i1]) { printf ("|"); }
             else { printf ("X"); }
           }
           printf ("\n");
           printf ("Sequence %d:\n%s\n", 1, msa[1].c_str());
           printf ("\n");
           exit(1);

         }
       }
       fflush(stdout);

       LOG_NOHEADER("\n");
       LOG_ALL("Batch checkpoint: Processed %ld windows and exported the consensus.\n", parameters.batch_of_windows);
//       LOG_ALL("Debug exiting.\n");
//       exit(1);
    }

    fprintf (fp_out_cons, "\n");
    fclose(fp_out_cons);

    ///////////////////////////////////////
    LOG_NOHEADER("\n");
    LOG_ALL("Processed %ld bp of %ld bp (100.00%%)\n", contig->get_data_length(), contig->get_data_length());
    LOG_NOHEADER("\n");
  }

  if (fp_alt_contig_path) { fclose(fp_alt_contig_path); }

  return 0;
}

int MajorityVoteFromMSALocal(std::string pir_path, std::string *cons) {
  SequenceFile pir(SEQ_FORMAT_FASTQ, pir_path);

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
//  std::vector<int64_t> offset_start;
//  std::vector<int64_t> offset_end;
//  offset_start.resize(seqs.size());
//  offset_end.resize(seqs.size());
//  for (int64_t i=0; i<seqs.size(); i++) {
//    offset_start[i] = 0;
//    offset_end[i] = seq_len;
//
////    for (int64_t j=0; j<seq_len; j++) {
////      if (seqs[i]->get_data()[j] != '-' && seqs[i]->get_data()[j] != '.') break;
////      offset_start[i] += 1;
////    }
////    for (int64_t j=(seq_len-1); j>=0; j--) {
////      if (seqs[i]->get_data()[j] != '-' && seqs[i]->get_data()[j] != '.') break;
////      offset_end[i] -= 1;
////    }
//  }

//  printf ("seq_len = %ld\n", seq_len);
//  for (int64_t i=0; i<seqs.size(); i++) {
//    printf ("offset_start[%ld] = %ld, offset_end[%ld] = %ld\n", i, offset_start[i], i, offset_end[i]);
//  }



  for (int64_t i=0; i<seq_len; i++) {
    // Count occurances for the column.
    int32_t base_counts[256] = {0};
    for (int32_t j=0; j<seqs.size(); j++) {
//      if (i < offset_start[j] || i >= offset_end[j]) { continue; }
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

    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --retree 1 --maxiterate 0 --nofft --op 0 --ep 1 --quiet %s > %s", // AlignedBases           48306(99.60%)       47482(100.00%)  AvgIdentity                    96.87                96.87
                        parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), window_path.c_str(), msa_path.c_str()).c_str());
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
