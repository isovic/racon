/*
 * pileup.cc
 *
 *  Created on: May 16, 2016
 *      Author: isovic
 */

#include "pileup.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include <assert.h>
#include "spoa.hpp"
#include "graph.hpp"

Pileup::Pileup(const SingleSequence* ref, std::vector<const SingleSequence*>& ref_alns) : ref_(ref) {
  bases_.clear();
  bases_.resize(ref->get_sequence_length());

  for (int64_t i=0; i<ref->get_sequence_length(); i++) {
    bases_[i].ref_base = ref->get_data()[i];
  }

  for (int64_t i=0; i<ref_alns.size(); i++) {
    if ((i % 1000) == 0) { fprintf (stderr, "Adding alignment %d to the pileup...\n", i); }
    AddAlignment(ref, ref_alns[i], i);
  }
}

void Pileup::AddAlignment(const SingleSequence* ref, const SingleSequence *seq, int64_t qid) {
  auto& aln = seq->get_aln();
  auto& cigar = aln.get_cigar();
  for (int64_t j=0; j<cigar.size(); j++) {
    if (cigar[j].op == 'S' || cigar[j].op == 'H') { continue; }
    // Just for easier handling.
    int64_t rpos = cigar[j].pos_ref + aln.get_pos() - 1;     // pos_ value is 1-based, just like in the SAM file.
    int64_t qpos = cigar[j].pos_query;

    Event new_event;
    char op = cigar[j].op;
    new_event.t = (op == 'M' || op == '=' || op == 'X') ? (kEventBase) :
                  (op == 'I') ? (kEventInsertion) :
                  (op == 'D') ? (kEventDeletion) :
                  (kEventUnknown);

    // If there is an unknown alignment event...well, we are currently lazy and don't feel like handling it.
    assert(new_event.t != kEventUnknown);

    // ID of the read which generated the event.
    new_event.qid = qid;

    // Initialize the bases and qualities.
    if (new_event.t == kEventInsertion) {
      new_event.e = seq->GetSequenceAsString(qpos, qpos + cigar[j].count);
      if (seq->get_quality() != NULL) { new_event.q = seq->GetQualityAsString(qpos, qpos + cigar[j].count); }
      // Done with generating the event, push it int.
      bases_[rpos].events.push_back(new_event);
    } else if (new_event.t == kEventBase) {
      for (int32_t k=0; k<cigar[j].count; k++) {
        new_event.e = seq->GetSequenceAsString(qpos+k, qpos+k+1);
        if (seq->get_quality() != NULL) { new_event.q = seq->GetQualityAsString(qpos+k, qpos+k+1); }
            // Done with generating the event, push it in.
            bases_[rpos+k].events.push_back(new_event);
            bases_[rpos+k].coverage += 1;
            bases_[rpos+k].base_coverage += 1;
      }
    } else if (new_event.t == kEventDeletion) {
      new_event.e = ref->GetSequenceAsString(rpos, rpos + cigar[j].count);
      new_event.q = std::string(new_event.e.size(), ('!'-1));
      // Done with generating the event, push it int.
      bases_[rpos].events.push_back(new_event);

      new_event.t = kEventBase;
      for (int32_t k=1; k<=cigar[j].count; k++) {
        new_event.e = '-';
        if (seq->get_quality() != NULL) { new_event.q = ('!'-1); }
          // Done with generating the event, push it int.
          bases_[rpos+k].events.push_back(new_event);
          bases_[rpos+k].coverage += 1;
          bases_[rpos+k].del_coverage += 1;
      }
    }
  }
}

Pileup::~Pileup() {
  bases_.clear();
}

void Pileup::GenerateConsensus(int32_t cov_threshold, std::string& cons) {
  std::stringstream ss;

  int32_t temp = 0;

  for (int64_t i=0; i<bases_.size(); i++) {

    if (bases_[i].coverage < cov_threshold) {
//      ss << bases_[i].ref_base;
//      printf ("  bases_[i].coverage < cov_threshold\n");
//      fflush(stdout);
      continue;
    }

    int32_t bc[256] = {0};
    std::map<std::string, int32_t> ic;
    std::map<std::string, int32_t> dc;
    std::vector<std::string> insertions;
    int32_t bcov = 0, icov = 0, dcov = 0;

    for (int64_t j=0; j<bases_[i].events.size(); j++) {
      auto &event = bases_[i].events[j];
      if (event.t == kEventBase) {
        bc[(int32_t) event.e.c_str()[0]] += 1;
        bcov += 1;

      } else if (event.t == kEventInsertion) {
        insertions.push_back(event.e);
        auto it = ic.find(event.e);
        if (it == ic.end()) { ic[event.e] = 1; }
        else { it->second += 1; }
        icov += 1;

      } else if (event.t == kEventDeletion) {
        auto it = dc.find(event.e);
        if (it == dc.end()) { dc[event.e] = 1; }
        else { it->second += 1; }
        dcov += 1;
      }
    }

    bool base_is_del = false;

    // This part handles the current base calling. It considers deletions
    // as normal bases, and picks the majority vote.
    std::vector<int32_t> bc_actgd = {bc['A'] + bc['a'], bc['C'] + bc['c'], bc['T'] + bc['t'], bc['G'] + bc['g'], bc['-']};  // d stands for deletion
    char actgd[] = {'A', 'C', 'T', 'G', '-'};
    std::vector<size_t> bc_actgd_indices;
    ordered_sort_array(&bc_actgd[0], bc_actgd.size(), bc_actgd_indices);
    int32_t sum_all_bc = 0;
    for (int32_t j=0; j<bc_actgd.size(); j++) { sum_all_bc += bc_actgd[j]; }
    int32_t num_bases = sum_all_bc - bc['-'];
//    printf ("  sum_all_bc = %d, sum_bases = %d, bc['-'] = %d\n", sum_all_bc, sum_bases, bc['-']);
    int32_t num_dels = bc['-'];

//    if (bc['-'] > num_bases) {
    if (actgd[bc_actgd_indices.back()] == '-') {
//      printf ("  base_is_del = true\n");
//      fflush(stdout);
      base_is_del = true;
      // Pass. The base is a deletion. Do not handle this case.
    } else {
      // Call the base.
      bool is_base_called = false;
      for (int32_t j=(bc_actgd_indices.size()-1); j>=0; j--) {
        if ((actgd[bc_actgd_indices[j]]) != '-' && bc_actgd[bc_actgd_indices[j]] > 0) {
          ss << actgd[bc_actgd_indices[j]];
          is_base_called = true;
          temp += 1;
//          printf ("  + is_base_called = true, '%c', count = %d\n", (char) bc_indices[j], bc[bc_indices[j]]);
//          fflush(stdout);
          break;
        }
      }
      if (is_base_called == false) {
        ss << bases_[i].ref_base;
//        printf ("Tu sam 1!\n");
//        exit(1);
//        printf ("  - bases_[i].ref_base = '%c'\n", bases_[i].ref_base);
//        fflush(stdout);
      }
    }

//    if (i == 7173) {
//      exit(1);
//    }

//    printf ("[i = %ld] coverage = %d, base_coverage = %d, del_coverage = %d\n", i, bases_[i].coverage, bases_[i].base_coverage, bases_[i].del_coverage);
//    fflush(stdout);
//    printf ("  bcov = %d, bc['-'] = %d, sum_bases = %d, icov = %d, {A: %d, C: %d, T: %d, G: %d}\t", bcov, bc['-'], sum_bases, icov, bc['A'], bc['C'], bc['T'], bc['G']);
//    for (int64_t j=0; j<insertions.size(); j++) {
//      if (j > 0) { printf (", "); }
//      else { printf ("  "); }
//      printf ("%s", insertions[j].c_str());
//    }
//    printf ("\n");
//    fflush(stdout);

    if (base_is_del == false && icov > bases_[i].coverage) {
      auto graph = construct_partial_order_graph(insertions, SPOA::AlignmentParams(1, -1, -1, -1, (SPOA::AlignmentType) 0));
      graph->generate_consensus();
      std::vector<std::string> msa;
      graph->generate_msa(msa, false);
      std::string ins_cons;
      MajorityVoteFromMSA(msa, ins_cons);
      ss << ins_cons;
//      printf ("  insertion! '%s'\n", ins_cons.c_str());
//      fflush(stdout);
    }
  }

  cons = ss.str();

//  printf ("temp = %d\n", temp);
//  fflush(stdout);
}

void Pileup::Verbose(FILE* fp) {
  for (int64_t i=0; i<bases_.size(); i++) {
    fprintf (fp, "%s\t%ld\t%s\t%ld\t", ref_->get_header(), i, bases_[i].ref_base.c_str(), bases_[i].events.size());
    int32_t base_cov = 0;
    for (int64_t j=0; j<bases_[i].events.size(); j++) {
//      fprintf (fp, "\t(%d, %s)", bases_[i].events[j].t, bases_[i].events[j].e.c_str());
      if (bases_[i].events[j].t == kEventBase) {
        if (bases_[i].events[j].e == bases_[i].ref_base) {
          fprintf (fp, ".");
          base_cov += 1;
        } else if (bases_[i].events[j].e == "-") {
          // Pass the deletions for now.
        } else {
          fprintf (fp, "%s", bases_[i].events[j].e.c_str());
          base_cov += 1;
        }
      }
    }

    fprintf (fp, " (%d)\t", base_cov);

    int32_t del_cov = 0;
    for (int64_t j=0; j<bases_[i].events.size(); j++) {
      if (bases_[i].events[j].t == kEventBase && bases_[i].events[j].e == "-") {
        fprintf (fp, "*");
        del_cov += 1;
      }
    }

    fprintf (fp, " (%d)\t", del_cov);

    int32_t ins_cov = 0;
    for (int64_t j=0; j<bases_[i].events.size(); j++) {
      if (bases_[i].events[j].t == kEventInsertion) {
        fprintf (fp, "+%d%s", bases_[i].events[j].e.size(), bases_[i].events[j].e.c_str());
        ins_cov += 1;
      }
    }

    fprintf (fp, " (%d)\n", ins_cov);
  }
}

int Pileup::MajorityVoteFromMSA(std::vector<std::string> &msa, std::string &consensus) {
  if (msa.size() == 0) { return 1; }

  int64_t seq_len = msa[0].size();
  std::stringstream ss;

  for (int64_t i=0; i<seq_len; i++) {
    // Count occurrences for the column.
    int32_t base_counts[256] = {0};
    for (int32_t j=0; j<msa.size(); j++) {
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
