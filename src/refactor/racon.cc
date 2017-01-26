/*
 * racon.cc
 *
 *  Created on: January 18, 2017
 *      Author: Ivan Sovic
 */

#include "racon.h"
#include <assert.h>
#include <iostream>
#include <vector>
#include <future>
#include "utility/utility_general.h"
#include "log_system/log_system.h"
#include "overlaps.h"

namespace is {

std::shared_ptr<Racon> createRacon(const std::shared_ptr<Parameters> param) {
  return std::shared_ptr<Racon>(new Racon(param));
}

Racon::~Racon() {

}

Racon::Racon(const std::shared_ptr<Parameters> param) : param_(param), thread_pool_(thread_pool::createThreadPool(param->num_threads())) {

}

void Racon::CreateConsensus() {
  if (param_->overlap_format().isPaf() || param_->overlap_format().isMhap()) {
    RunFromOverlaps_();
  }
}

void Racon::RunFromOverlaps_() {
  // Parse the backbone.
  SequenceFile targets(SEQ_FORMAT_AUTO, param_->contigs_path());

  // Parse the reads.
  SequenceFile queries(SEQ_FORMAT_AUTO, param_->reads_path());
  // Sanity check to see if the reads have quality values.
  if (queries.HasQV() == false) {
    FATAL_REPORT(ERR_WRONG_FILE_TYPE, "ERROR: Reads are not specified in a format which contains quality information. Exiting.\n");
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

  AlignAndGenerateWindows_(queries, targets, overlaps);

//  MapOverlapRange contig_overlaps;
//  FindContigOverlaps_(overlaps, contig_overlaps);

//  // Process each contig individually.
//  auto& tseqs = targets.get_sequences();
//  for (int64_t i=0; i<tseqs.size(); i++) {
//    auto t = tseqs[i];
//    std::string tname = TrimToFirstSpace(std::string(t->get_header()));
//
//    // Retrieve all overlaps for the current target.
//    auto it = contig_overlaps.find(i);
//
//    if (it == contig_overlaps.end()) {  // This target has no overlaps. Put it in a special place.
//      fprintf (stderr, "TODO: targets without overlaps not handled yet!\n");
//      fflush(stderr);
//      exit(1);
//    }
//
////    AlignOverlaps_();
////    CreateIntervalTree_();
////    PopulateJobs_();
////    wait for threads to finish
//  }
}

void Racon::RunFromAlignments_() {
}

int Racon::AlignAndGenerateWindows_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps) {
  // Create storage for return values.
  std::vector<std::future<int>> futures;

  for (int64_t i=0; i<queries.get_sequences().size(); i++) {
//    thread_futures.emplace_back(thread_pool->submit_task(function1, std::ref(data), index, ...)); // be sure to use std::ref() when passing references!
  }

  // Wait for threads to finish
  for (auto& it: futures) {
      it.wait();
  }

  return 0;
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
    if (overlaps[i].Bid() == ctg_id) {        // Target is the same, just move the end.
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
