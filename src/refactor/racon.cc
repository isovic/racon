/*
 * racon.cc
 *
 *  Created on: January 18, 2017
 *      Author: Ivan Sovic
 */

#include "racon.h"
#include <assert.h>
#include <iostream>
#include "utility/utility_general.h"
#include "log_system/log_system.h"
#include "overlaps.h"

namespace is {

std::shared_ptr<Racon> createRacon(const std::shared_ptr<Parameters> param) {
  return std::shared_ptr<Racon>(new Racon(param));
}

Racon::~Racon() {
}

Racon::Racon(const std::shared_ptr<Parameters> param) : param_(param) {
}

Racon::Racon(const SequenceFile& reads, const SequenceFile& targets) {
}

void Racon::CreateConsensus() {
  if (param_->overlap_format().isPaf() || param_->overlap_format().isMhap()) {
    RunFromOverlaps_();
  }
}

//void Racon::PopulateJobsConsensus_(const SequenceFile &refs, int64_t win_len, JobQueue &jobs) const {
//  assert(win_len > 0);
//  for (int64_t i=0; i<((int64_t) refs.get_sequences().size()); i++) {
//  	auto s = refs.get_sequences()[i];
//  	int64_t s_len = (int64_t) s->get_sequence_length();
//    for (int64_t j=0; j<s_len; j+=win_len) {
//      int64_t start = j;
//      int64_t end = std::min(start + win_len, s_len);
//   	  jobs.push_back(createJob(i, start, end, win_len));
//    }
//  }
//}

//void Racon::PopulateJobsErc_(const SequenceFile &refs, int64_t win_len, JobQueue &jobs) const {
//  assert(win_len > 0);
//  for (int64_t i=0; i<((int64_t) refs.get_sequences().size()); i++) {
//    auto s = refs.get_sequences()[i];
//    int64_t s_len = (int64_t) s->get_sequence_length();
//    jobs.push_back(createJob(i, 0, s_len, win_len));
//  }
//}

void Racon::RunFromOverlaps_() {
  // Parse the backbone.
  SequenceFile targets(SEQ_FORMAT_AUTO, param_->contigs_path());

  // Parse the reads.
  SequenceFile queries(SEQ_FORMAT_AUTO, param_->reads_path());
  // Sanity check to see if the reads have quality values.
  if (queries.HasQV() == false) {
    fprintf (stderr, "ERROR: Reads are not specified in a format which contains quality information. Exiting.\n");
    exit(1);
  }

  // Hash the sequences by their name.
  MapId query_id, target_id;
  HashNames_(queries, query_id);
  HashNames_(targets, target_id);

  // Load the overlaps.
  LOG_ALL("Using %s for input alignments. (%s)\n",
          (param_->overlap_format().isPaf()) ? "PAF" : "MHAP", param_->aln_path().c_str())
  LOG_ALL("Started parsing the overlaps file.\n");

  Overlaps overlaps(param_->aln_path(), param_->overlap_format(), query_id, target_id, param_->error_rate(), param_->do_erc());
  overlaps.SortByTargetId();

  MapOverlapRange contig_overlaps;
  FindContigOverlaps_(overlaps, contig_overlaps);

  // Process each contig individually.
  auto& tseqs = targets.get_sequences();
  for (int64_t i=0; i<tseqs.size(); i++) {
    auto t = tseqs[i];
    std::string tname = TrimToFirstSpace(std::string(t->get_header()));

    // Retrieve all overlaps for the current target.
    auto it = contig_overlaps.find(i);

    if (it == contig_overlaps.end()) {  // This target has no overlaps. Put it in a special place.
      fprintf (stderr, "TODO: targets without overlaps not handled yet!\n");
      fflush(stderr);
      exit(1);
    }

//    AlignOverlaps_();
//    CreateIntervalTree_();
//    PopulateJobs_();
//    wait for threads to finish
  }

//  if (parameters.do_sparse == false || parameters.do_erc) {
//    LOG_ALL("Overlaps will be fully aligned.\n");
//    ConsensusFromOverlaps(parameters, seqs_gfa, seqs_reads, qname_to_ids, rname_to_ids, overlaps_final);
//  }
}

void Racon::RunFromAlignments_() {
}

void Racon::HashNames_(const SequenceFile &seqs, MapId &id) const {
  for (size_t i=0; i<seqs.get_sequences().size(); i++) {
    const auto& s = seqs.get_sequences()[i];
    std::string header = std::string(s->get_header());
    id[header] = i;
    id[TrimToFirstSpace(header)] = i;
    std::size_t found = header.find(":");
    id[header.substr(0, found)] = i;
  }
}

int Racon::FindContigOverlaps_(const Overlaps &sorted_overlaps, MapOverlapRange &contig_overlaps) const {
  // Brevity.
  auto& overlaps = sorted_overlaps.overlaps();

  // Sanity check.
  if (overlaps.size() == 0) { return 1; }

  // Make sure it's clear.
  contig_overlaps.clear();

  Range ctg_loc;
  auto ctg_id = overlaps.front().Bid();
  for (int64_t i=1; i<overlaps.size(); i++) {
    if (overlaps[i].Bid() == ctg_id) {   // Target is the same, just move the end.
      ctg_loc.end = i;
    } else {                                  // Different target, add to map and restart.
      contig_overlaps[ctg_id] = ctg_loc;
      ctg_loc.start = ctg_loc.end = i;
      ctg_id = overlaps[i].Bid();
    }
  }

  // Last update and push the last streak to the map.
  contig_overlaps[overlaps.back().Bid()] = ctg_loc;

  // All went fine.
  return 0;
}

} /* namespace is */
