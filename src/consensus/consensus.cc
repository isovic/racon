/*
b * consensus.cc
 *
 *  Created on: Feb 15, 2016
 *      Author: isovic
 */

#include "consensus/consensus.h"
#include "../../codebase/poa/src/poa.hpp"
#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include <stdint.h>
#include <algorithm>
#include <sstream>
#include <stdlib.h>

int AlignmentsToContigs(const SequenceFile &alns, std::vector<std::string> &ctg_names, std::map<std::string, std::vector<const SingleSequence *> > &ctg_alns) {
  ctg_names.clear();
  ctg_alns.clear();

  for (int64_t i=0; i<alns.get_sequences().size(); i++) {
    if (alns.get_sequences()[i]->get_aln().IsMapped() == false) continue;

    auto it = ctg_alns.find(alns.get_sequences()[i]->get_aln().rname);
    if (it != ctg_alns.end()) {
      it->second.push_back((const SingleSequence *) (alns.get_sequences()[i]));
    } else {
      ctg_alns[alns.get_sequences()[i]->get_aln().rname] = std::vector<const SingleSequence *> {(const SingleSequence *) alns.get_sequences()[i]};
      ctg_names.push_back(alns.get_sequences()[i]->get_aln().rname);
    }
  }

  return 0;
}

int ExtractAltContigs(std::vector<const SingleSequence *> &ctg_alns, int64_t raw_ctg_len, double coverage_threshold, double percent_overlap, double qv_threshold, std::vector<std::vector<const SingleSequence *> *> &ret_alt_contigs, std::vector<const SingleSequence *> &rejected_alns) {
  ret_alt_contigs.clear();

  // This sorts ascending by the pos field.
  std::sort(ctg_alns.begin(), ctg_alns.end(), seqaln_sort_key());

  // Initialize the current set of alignments to process the full set of mapped reads.
  // Skip the unmapped reads. Also, skip reads with base qualities below some threshold.
  std::vector<const SingleSequence *> alns_to_process;
  alns_to_process.reserve(ctg_alns.size());
  for (int64_t i=0; i<ctg_alns.size(); i++) {
    double average_bq = ctg_alns[i]->CalcAverageBQ();
    if (average_bq >= 0 && average_bq < qv_threshold) { continue; }
//    printf ("ctg_alns[i]->CalcAverageBQ() = %f\n", ctg_alns[i]->CalcAverageBQ());

    if (ctg_alns[i]->get_aln().IsMapped()) { alns_to_process.push_back(ctg_alns[i]); }
  }

  // Hash all the alignment lengths (which will be used a lot).
  std::map<const SingleSequence *, int64_t> aln_ref_lens;
  for (int64_t i=0; i<alns_to_process.size(); i++) {
    aln_ref_lens[alns_to_process[i]] = alns_to_process[i]->get_aln().GetReferenceLengthFromCigar();
  }

  int32_t coverage = 0;

  while ((MAX_COV <= 0 || (MAX_COV > 0 && coverage < MAX_COV)) and alns_to_process.size() > 0) {
    coverage += 1;
    LOG_ALL("Coverage: %ld, alns_to_process.size() = %ld\n", coverage, alns_to_process.size());

    // Reserve space for the alternate contig.
    std::vector<const SingleSequence *> *new_alt_ctg_alns = new std::vector<const SingleSequence *>;
    new_alt_ctg_alns->reserve(ctg_alns.size());

    // Alignments which weren't used in this iteration.
    std::vector<const SingleSequence *> unused_alns;
    unused_alns.reserve(ctg_alns.size());

    // Initialize some counters.
    int64_t aln_ref_len = aln_ref_lens[alns_to_process[0]];
    int64_t new_ctg_bases = aln_ref_len;                // Number of bases on the reference (original raw contig) covered by new alignments.
    int64_t old_ctg_bases = alns_to_process[0]->get_aln().pos - 1;         // Number of bases in between alignments (or, those comming from the original raw contig).
    int64_t previous_start = alns_to_process[0]->get_aln().pos - 1;
    int64_t previous_end = previous_start + aln_ref_len;
    int64_t prev_ref_len = aln_ref_len;
    int64_t prev_candidate_id = 0;
    new_alt_ctg_alns->push_back(alns_to_process[0]);

    for (int64_t i=0; i<alns_to_process.size(); i++) {
      const SingleSequence *candidate = alns_to_process[i];
      const SequenceAlignment &candidate_aln = alns_to_process[i]->get_aln();
      aln_ref_len = aln_ref_lens[alns_to_process[i]];

      if ((candidate_aln.pos - 1) >= (previous_end - prev_ref_len * percent_overlap) and ((candidate_aln.pos - 1) + aln_ref_len) > previous_end) {
        new_ctg_bases += ((candidate_aln.pos - 1) >= previous_end) ? (aln_ref_len) : (aln_ref_len - (previous_end - candidate_aln.pos + 1));
        old_ctg_bases += ((candidate_aln.pos - 1) >= previous_end) ? (candidate_aln.pos - 1 - previous_end) : 0;
        unused_alns.insert(unused_alns.end(), alns_to_process.begin() + (prev_candidate_id + 1), alns_to_process.begin() + i);
        prev_candidate_id = i;
        previous_start = candidate_aln.pos - 1;
        previous_end = previous_start + aln_ref_len;
        prev_ref_len = aln_ref_len;
        new_alt_ctg_alns->push_back(candidate);
      }
    }
    unused_alns.insert(unused_alns.end(), alns_to_process.begin() + (prev_candidate_id + 1), alns_to_process.end());
    old_ctg_bases += raw_ctg_len - previous_end;

    // Check the threshold for the covered bases. If above threshold, accept the contig.
    if ((((float) new_ctg_bases) / ((float) new_ctg_bases + old_ctg_bases)) < coverage_threshold) {
#ifdef EXTREME_DEBUG_
      printf ("[REJECTED]:\n");
      printf ("new_ctg_bases = %ld\n", new_ctg_bases);
      printf ("new_ctg_bases + old_ctg_bases = %ld\n", new_ctg_bases + old_ctg_bases);
      printf ("contig_coverage = %f\n", (((float) new_ctg_bases) / ((float) new_ctg_bases + old_ctg_bases)));
      printf ("\n");
#endif

      unused_alns.clear();
      unused_alns.insert(unused_alns.end(), alns_to_process.begin() + 1, alns_to_process.end());
      rejected_alns.push_back(alns_to_process[0]);
      new_alt_ctg_alns->clear();
      if (new_alt_ctg_alns) delete new_alt_ctg_alns;
      new_alt_ctg_alns = NULL;

    } else {
#ifdef EXTREME_DEBUG_
      printf ("\t[ACCEPTED]:\n");
      printf ("\tcontig_coverage = %f\n", (((float) new_ctg_bases) / ((float) new_ctg_bases + old_ctg_bases)));
#endif
      ret_alt_contigs.push_back(new_alt_ctg_alns);
    }

    alns_to_process = unused_alns;
    unused_alns.clear();

  }

  return 0;
}

void VerboseExtractedAlignments(std::vector<std::vector<const SingleSequence *> *> &alt_contigs, std::vector<const SingleSequence *> rejected_alns, std::string out_alt_ctg_path, std::string out_rejected_path) {

  FILE *fp_alt_contigs = fopen(out_alt_ctg_path.c_str(), "w");
  for (int64_t j=0; j<alt_contigs.size(); j++) {
    fprintf (fp_alt_contigs, "Contig %ld\n", j);
    for (int64_t k=0; k<alt_contigs[j]->size(); k++) {
      const SingleSequence *seq = (*alt_contigs[j])[k];
//        fprintf (fp_alt_contigs, "%s\n", seq->MakeSAMLine().c_str());
      fprintf (fp_alt_contigs, "%s\t%ld\t%ld\n", seq->get_header(), (seq->get_aln().pos - 1), (seq->get_aln().pos - 1 + seq->get_aln().GetReferenceLengthFromCigar()));
    }
    fprintf (fp_alt_contigs, "\n");
  }
  fclose(fp_alt_contigs);

  FILE *fp_rejected_alns = fopen(out_rejected_path.c_str(), "w");
  for (int64_t j=0; j<rejected_alns.size(); j++) {
    const SingleSequence *seq = (rejected_alns[j]);
//        fprintf (fp_rejected_alns, "%s\n", seq->MakeSAMLine().c_str());
    fprintf (fp_rejected_alns, "%s\t%ld\t%ld\n", seq->get_header(), (seq->get_aln().pos - 1), (seq->get_aln().pos - 1 + seq->get_aln().GetReferenceLengthFromCigar()));
  }
  fclose(fp_rejected_alns);
}

int GetCigarBetweenPositions(const SingleSequence *seq_aln, int64_t start_seq_pos, int64_t start_cig_id, int64_t end_seq_pos, int64_t end_cig_id, std::vector<CigarOp> *ret_cigar) {
  // Check if the beginning CIGAR operation is fully captured, or we need to split it and take only the right part.
  if (start_cig_id == end_cig_id) {
    int64_t dist_start = start_seq_pos - seq_aln->get_aln().cigar[start_cig_id].pos_query;
    int64_t dist_end = end_seq_pos - seq_aln->get_aln().cigar[end_cig_id].pos_query;
    CigarOp partial_op = seq_aln->get_aln().cigar[start_cig_id];
    partial_op.count = (end_seq_pos - start_seq_pos);
    if (is_cigar_read(partial_op.op)) { partial_op.pos_query += dist_start; }
    if (is_cigar_ref(partial_op.op)) { partial_op.pos_ref += dist_start; }
    ret_cigar->insert(ret_cigar->end(), &seq_aln->get_aln().cigar[0] + start_cig_id + 1, &seq_aln->get_aln().cigar[0] + end_cig_id - 1);
    return 0;
  }

  int64_t dist_start = start_seq_pos - seq_aln->get_aln().cigar[start_cig_id].pos_query;
  CigarOp first_op = seq_aln->get_aln().cigar[start_cig_id];
  first_op.count -= dist_start;
  if (is_cigar_read(first_op.op)) { first_op.pos_query += dist_start; }
  if (is_cigar_ref(first_op.op)) { first_op.pos_ref += dist_start; }
  if (first_op.count > 0) {
    ret_cigar->push_back(first_op);
  }

  if (end_cig_id > (start_cig_id + 1)) {
    ret_cigar->insert(ret_cigar->end(), &seq_aln->get_aln().cigar[0] + start_cig_id + 1, &seq_aln->get_aln().cigar[0] + end_cig_id);
  }

  int64_t dist_end = end_seq_pos - seq_aln->get_aln().cigar[end_cig_id].pos_query + 1;
  CigarOp last_op = seq_aln->get_aln().cigar[end_cig_id];
  last_op.count = dist_end;
  if (last_op.count > 0) {
    ret_cigar->push_back(last_op);
  }

  return 0;
}

int ConstructContigFromAlns(const SingleSequence &orig_contig, const std::vector<const SingleSequence *> *seq_alns, const std::map<const SingleSequence *, int64_t> &aln_ref_lens, SingleSequence &new_contig) {
  //, std::vector<int8_t> *new_contig, std::vector<CigarOp> *new_contig_cigar) {
  LOG_ALL("Constructing an alternate contig sequence from a set of alignments.\n");

  std::vector<int8_t> new_contig_seq;
  new_contig_seq.reserve(orig_contig.get_data_length());

  SequenceAlignment &new_contig_aln = (SequenceAlignment &) new_contig.get_aln();

  new_contig.Clear();
  new_contig.InitHeader(std::string("AlternateContig"));
  new_contig_aln.rname = TrimToFirstSpace(orig_contig.get_header()).c_str();
  new_contig_aln.pos = 1;
  new_contig_aln.flag = 0;
  new_contig_aln.cigar.clear();
  new_contig_aln.as = orig_contig.get_sequence_length();
  new_contig_aln.mapq = 40;
  new_contig_aln.evalue = 0.0;

  new_contig_aln.cigar.reserve(orig_contig.get_data_length()/2 + 1);

  int64_t prev_end_ref = 0;
  for (int64_t i=0; i<seq_alns->size(); i++) {
    const SingleSequence *seq_aln = seq_alns->at(i);
    const SequenceAlignment &aln = seq_aln->get_aln();
    int64_t pos = aln.pos;

    int64_t start_ref = ((pos - 1) < prev_end_ref) ? (prev_end_ref) : (pos - 1);
    int64_t end_ref = (pos - 1) + aln_ref_lens.find(seq_alns->at(i))->second - 1;

    int64_t start_cig_id = 0, end_cig_id = 0;
    int64_t start_seq = seq_aln->get_aln().FindBasePositionOnRead(aln.cigar, start_ref, &start_cig_id);
    int64_t end_seq = seq_aln->get_aln().FindBasePositionOnRead(aln.cigar, end_ref, &end_cig_id);

#ifdef EXTREME_DEBUG_
    LOG_DEBUG_NOHEADER("  [%ld] new_contig->size() = %ld\n", i, new_contig_seq.size());
#endif

    // If there was a gap between the current and the previous read, fill it with the corresponding part of the given contig_raw.
    // Also, if the first sam_line did not start at position 0, the leading chunk will be filled in.
    if (start_ref > prev_end_ref) {
#ifdef EXTREME_DEBUG_
      LOG_DEBUG_NOHEADER("  [%ld] Filling gap: start_ref = %ld, end_ref = %ld, prev_end_ref = %ld, start_seq = %ld, end_seq = %ld\n", i, start_ref, end_ref, prev_end_ref, start_seq, end_seq);
      LOG_DEBUG_NOHEADER("      Inserting sequence from: %ld to: %ld\n", prev_end_ref, start_ref);
#endif

      CigarOp match_op;
      match_op.op = '=';  match_op.count = start_ref - prev_end_ref - 1;  match_op.pos_query = new_contig_seq.size();  match_op.pos_ref = prev_end_ref;
      if (match_op.count > 0) { new_contig_aln.cigar.push_back(match_op); }

      new_contig_seq.insert(new_contig_seq.end(), orig_contig.get_data() + prev_end_ref + 1 , orig_contig.get_data() + start_ref);
#ifdef EXTREME_DEBUG_
      LOG_DEBUG_NOHEADER("  [%ld] new_contig->size() = %ld\n", i, new_contig_seq.size());
#endif
    }

    // Copy the new sequence into the alternate contig.
#ifdef EXTREME_DEBUG_
    LOG_DEBUG_NOHEADER("  [%ld] Adding SAM: start_ref = %ld, end_ref = %ld, prev_end_ref = %ld, start_seq = %ld, end_seq = %ld, new_contig_len (before) = %ld, new_contig_len (after) = %ld, len(ref) = %ld, len(seq) = %ld\n", i, start_ref, end_ref, prev_end_ref, start_seq, end_seq, new_contig_seq.size(), (new_contig_seq.size() + (end_seq - start_seq)), orig_contig.get_data_length(), seq_aln->get_data_length());
#endif
    GetCigarBetweenPositions(seq_aln, start_seq, start_cig_id, end_seq, end_cig_id, &new_contig_aln.cigar);
    new_contig_seq.insert(new_contig_seq.end(), seq_aln->get_data() + start_seq, seq_aln->get_data() + end_seq + 1);
#ifdef EXTREME_DEBUG_
    LOG_DEBUG_NOHEADER("  [%ld] new_contig->size() = %ld\n", i, new_contig_seq.size());
#endif

    // Check if this was the last sam_line. If it was, and it didn't end on the end of the reference, fill in the missing part on the back.
    if ((i + 1) == seq_alns->size()) {
      if (end_ref < orig_contig.get_data_length()) {
#ifdef EXTREME_DEBUG_
        LOG_DEBUG_NOHEADER("  [%ld] Fixing the end: (I) with contig_raw: start_ref = %ld, end_ref = %ld, prev_end_ref = %ld\n", i, start_ref, end_ref, prev_end_ref);
#endif

        CigarOp match_op;
        match_op.op = '=';  match_op.count = orig_contig.get_data_length() - end_ref - 1;  match_op.pos_query = new_contig_seq.size();  match_op.pos_ref = end_ref;
        if (match_op.count > 0) { new_contig_aln.cigar.push_back(match_op); }
        new_contig_seq.insert(new_contig_seq.end(), orig_contig.get_data() + end_ref + 1, orig_contig.get_data() + orig_contig.get_data_length());

#ifdef EXTREME_DEBUG_
        LOG_DEBUG_NOHEADER("  [%ld] match_op.count = %ld, orig_contig.get_data_length() = %ld, end_ref = %ld\n", i, match_op.count, orig_contig.get_data_length(), end_ref);
        LOG_DEBUG_NOHEADER("  [%ld] new_contig->size() = %ld\n", i, new_contig_seq.size());
#endif
      } else {
#ifdef EXTREME_DEBUG_
        LOG_DEBUG_NOHEADER("  [%ld] Fixing the end: (II) with contig_raw: start_ref = %ld, end_ref = %ld, prev_end_ref = %ld", i, start_ref, end_ref, prev_end_ref);
#endif
        GetCigarBetweenPositions(seq_aln, end_seq, end_cig_id, seq_aln->get_data_length(), seq_aln->get_aln().cigar.size() - 1, &new_contig_aln.cigar);
        new_contig_seq.insert(new_contig_seq.end(), seq_aln->get_data() + end_seq + 1, seq_aln->get_data() + seq_aln->get_data_length());
#ifdef EXTREME_DEBUG_
        LOG_DEBUG_NOHEADER("  [%ld] new_contig->size() = %ld\n", i, new_contig_seq.size());
#endif
      }
    }
    prev_end_ref = end_ref;
#ifdef EXTREME_DEBUG_
    LOG_DEBUG_NOHEADER("  [%ld] new_contig->size() = %ld\n", i, new_contig_seq.size());
    LOG_DEBUG_NOHEADER("\n");
#endif
  }

  new_contig.InitDatafromAscii(&new_contig_seq[0], (uint64_t) new_contig_seq.size());

  int64_t pos_ref = new_contig_aln.pos - 1, pos_query = 0;
  for (int64_t i=0; i<new_contig_aln.cigar.size(); i++) {
    new_contig_aln.cigar[i].pos_query = pos_query;
    new_contig_aln.cigar[i].pos_ref = pos_ref;
    if (is_cigar_ref(new_contig_aln.cigar[i].op)) pos_ref += new_contig_aln.cigar[i].count;
    if (is_cigar_read(new_contig_aln.cigar[i].op)) pos_query += new_contig_aln.cigar[i].count;
  }

#ifdef EXTREME_DEBUG_
  LOG_DEBUG_NOHEADER("new_contig->size() = %ld\n", (new_contig_seq.size()));
#endif

  LOG_DEBUG("Finished constructing the alternate contig sequence.\n");

  return 0;
}


int ConsensusFromMSA(std::string pir_path, std::string *cons) {
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

  for (int64_t i=0; i<seq_len; i++) {
    // Count occurances for the column.
    int32_t base_counts[256] = {0};
    for (int32_t j=0; j<seqs.size(); j++) {
      base_counts[toupper(seqs[j]->get_data()[i])] += 1;
    }
    // Find the majority vote.
    int8_t max_base = (int8_t) 'A';
    for (int32_t j=0; j<256; j++) {
      if (base_counts[j] > base_counts[max_base]) max_base = (char) j;
    }

    if (max_base != '.' && max_base != '-') {
      ss << (char) max_base;
    }
  }

  *cons = ss.str();

  return 0;
}

int RunMSAFromSystem(const ProgramParameters &parameters, std::string &cons) {
  std::string msa_path = FormatString("%s.msa", parameters.temp_window_path.c_str());
  if (parameters.msa == "mafft") {
//    system(FormatString("export MAFFT_BINARIES=/home/isovic/work/eclipse-workspace/git/consise/tools/mafft-7.273-with-extensions/binaries/; /home/isovic/work/eclipse-workspace/git/consise/tools/mafft-7.273-with-extensions/scripts/mafft --op 0 --ep 1 --quiet %s > %s",
//                        parameters.temp_window_path.c_str(), msa_path.c_str()).c_str());
    int32_t rc = system(FormatString("export MAFFT_BINARIES=$PWD/%s/%s/binaries/; %s/%s/scripts/mafft --op 0 --ep 1 --quiet %s > %s",
                        parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), parameters.program_folder.c_str(), parameters.mafft_folder.c_str(), parameters.temp_window_path.c_str(), msa_path.c_str()).c_str());
  } else if (parameters.msa == "poav2") {
//    system(FormatString("/home/isovic/work/eclipse-workspace/git/consise/tools/poaV2/poa -do_global -do_progressive -read_fasta %s -pir %s /home/isovic/work/eclipse-workspace/git/consise/test-data/all1.mat",
//                        parameters.temp_window_path.c_str(), msa_path.c_str()).c_str());
    int32_t rc = system(FormatString("%s/%s/poa -do_global -do_progressive -read_fasta %s -pir %s %s/../settings/all1-poav2.mat",
                        parameters.program_folder.c_str(), parameters.poav2_folder.c_str(), parameters.temp_window_path.c_str(), msa_path.c_str(), parameters.program_folder.c_str()).c_str());
  } else {
    ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Unrecognized MSA option '%s'! Exiting.",parameters. msa.c_str());
    return 1;
  }

  ConsensusFromMSA(msa_path, &cons);

  return 0;
}

int Consensus(const ProgramParameters &parameters, const SequenceFile &contigs, const SequenceFile &alns) {
  LOG_ALL("Running consensus.\n");

  std::vector<std::string> ctg_names;
  std::map<std::string, std::vector<const SingleSequence *> > all_ctg_alns;

  // Separate alignments into groups for each contig.
  LOG_ALL("Separating alignments to individual contigs.\n");
  AlignmentsToContigs(alns, ctg_names, all_ctg_alns);
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

  // Clear the output file.
  FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "w");
  fclose(fp_out_cons);

  // For each contig, process alignments to obtain alternate contigs.
  for (int32_t i=0; i<ctg_names.size(); i++) {
    auto it = all_ctg_alns.find(ctg_names[i]);
    if (it == all_ctg_alns.end()) {
      FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Something strange happened. Contig name, which was extracted from alignments, cannot be found in the std::map containing those same alignments.");
      // Exits.
    }

    const SingleSequence *current_contig = rname_to_seq[ctg_names[i]];

    // Get alternate contigs in form of vectors of alignments. Each vector presents one alternate contig.
    std::vector<const SingleSequence *> &ctg_alns = it->second;
    std::vector<std::vector<const SingleSequence *> *> alt_contigs;
    std::vector<const SingleSequence *> rejected_alns;

    ExtractAltContigs(ctg_alns, rname_to_seq.find(it->first)->second->get_sequence_length(), 0.80, 0.01, parameters.qv_threshold, alt_contigs, rejected_alns);

    // Generate the alternate contig sequences from the sets of alignments (alt_contigs).
    SingleSequence alt_contig_seqs[alt_contigs.size()];
    for (int32_t j=0; j<alt_contigs.size(); j++) {
      LOG_DEBUG("Constructing contig %d / %d...\n", (j + 1), alt_contigs.size());

      ConstructContigFromAlns(*rname_to_seq[ctg_names[i]], (const std::vector<const SingleSequence *> *) alt_contigs[j], aln_ref_lens, alt_contig_seqs[j]);
      alt_contig_seqs[j].InitHeader(FormatString("AlternateContig_%d", j));

      LOG_DEBUG("Done with construction of alternate contig %d / %d.\n\n", (j + 1), alt_contigs.size());

      // Debug output. Verbose the alternate contigs with their alignments to an output temp FASTA file.
      if (fp_alt_contig_path) {
        LOG_DEBUG("Writing the alternate contigs to file: '%s'.\n", parameters.temp_alt_contig_path.c_str());
        fprintf (fp_alt_contig_path, "%s\n", alt_contig_seqs[j].MakeSAMLine().c_str());
      }
    }

    LOG_DEBUG("Alternate contigs generated. Starting to process windows.\n");
    LOG_NOHEADER("\n");

    ///////////////////////////////////////
    /// At this point we have obtained all alternate contig sequences.
    /// Now we need to process them with a sliding (non-overlapping) window and POA.
    ///////////////////////////////////////
    FILE *fp_out_cons = fopen(parameters.consensus_path.c_str(), "a");
    fprintf (fp_out_cons, ">Consensus_%d %s\n", i, ctg_names[i].c_str());
    for (int64_t window_start = 0; window_start < current_contig->get_sequence_length(); window_start += parameters.window_len) {
      int64_t window_end = window_start + parameters.window_len;
      LOG_ALL("\rProcessing window: %ld bp to %ld bp (%.2f%%)", window_start, window_end, 100.0 * ((float) window_start / (float) current_contig->get_data_length()));

      FILE *fp_window = NULL;

      if (parameters.temp_window_path != "") {
        fp_window = fopen(parameters.temp_window_path.c_str(), "w");
      }

      std::vector<const char *> sequences_for_poa;
      std::vector<int64_t> sequences_for_poa_lengths;
      for (int64_t j=0; j<alt_contigs.size(); j++) {
        const SequenceAlignment &aln = alt_contig_seqs[j].get_aln();
        int64_t start_seq = aln.FindBasePositionOnRead(aln.cigar, window_start);
        int64_t end_seq = aln.FindBasePositionOnRead(aln.cigar, window_end);

        // If the last requested position was out of bounds of the sequence, position it at the last base.
        if (end_seq == -1) {
          end_seq = (aln.pos - 1) + aln.GetReferenceLengthFromCigar() - 1;  // The extra -1 is to point to the last base (and not one after it).
        }

        if ((end_seq - start_seq) <= 0) { continue; };

        sequences_for_poa.push_back((const char*) (alt_contig_seqs[j].get_data() + start_seq));
        sequences_for_poa_lengths.push_back(end_seq - start_seq + 1);

        if (fp_window) {
          fprintf (fp_window, ">%s__Window_%d_to_%d\n", alt_contig_seqs[j].get_header(), window_start, window_end);
          fprintf (fp_window, "%s\n", GetSubstring((char *) sequences_for_poa.back(), sequences_for_poa_lengths.back()).c_str());
        }
      }

      if (parameters.msa == "poa") {
        std::vector<std::string> poa_sequences;
        for (auto s = 0; s < sequences_for_poa.size(); ++s) {
            poa_sequences.emplace_back(sequences_for_poa[s], sequences_for_poa_lengths[s]);
        }
        auto consensus_sequence = POA::poa_consensus(poa_sequences);
        fprintf (fp_out_cons, "%s", consensus_sequence.c_str());

      } else if (fp_window) {
        fclose(fp_window);
        std::string cons;
        std::string a, b;
        RunMSAFromSystem(parameters, cons);
        fprintf (fp_out_cons, "%s", cons.c_str());
      } else {
        ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Window file not opened!\n");
      }

    }
    fprintf (fp_out_cons, "\n");
    fclose(fp_out_cons);
//    exit(1);
    ///////////////////////////////////////
    LOG_NOHEADER("\n");
    LOG_ALL("Processed %ld bp of %ld bp (100.00%%)\n", current_contig->get_data_length(), current_contig->get_data_length());
    LOG_NOHEADER("\n");

    VerboseExtractedAlignments(alt_contigs, rejected_alns, "temp/alt_contigs.csv", "temp/rejected_alns.csv");
//    VerboseExtractedAlignments(alt_contigs1, rejected_alns1, "temp/alt_contigs1.csv", "temp/rejected_alns1.csv");

    for (int64_t j=0; j<alt_contigs.size(); j++) {
      if (alt_contigs[j]) delete alt_contigs[j];
      alt_contigs[j] = NULL;
    }
    alt_contigs.clear();
  }

  if (fp_alt_contig_path) { fclose(fp_alt_contig_path); }

  return 0;
}
