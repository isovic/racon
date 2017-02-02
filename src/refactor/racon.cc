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
#include "alignment.h"
#include "spoa.hpp"
#include "graph.hpp"

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

  Overlaps overlaps(param_->aln_path(), param_->overlap_format(), query_id, target_id, param_->error_rate(), param_->do_erc() == false);
  overlaps.SortByTargetId();


  // Storage for sampled overlaps. This is actually all we need from the alignments.
  // Every position on the reference which is a multiple of param_->window_len() will
  // be stored, as well as 1 base before and 1 base after. If window_ext > 0, then
  // positions around k*window_len +- window_ext will also be stored.
  std::vector<std::shared_ptr<SampledOverlap>> sampled;

  AlignAndSampleOverlaps_(queries, targets, overlaps, sampled);

  ConstructWindows_(targets, overlaps, sampled, windows_);

  RunAllJobs_(queries, targets, overlaps, windows_);

  LOG_ALL("Done!\n");

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

int Racon::AlignAndSampleOverlaps_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps, std::vector<std::shared_ptr<SampledOverlap>> &sampled) {
  LOG_ALL("Aligning overlaps.\n");

  int64_t n_overlaps = overlaps.overlaps().size();                        // Number of overlaps to process.
  int64_t window_ext = param_->window_len() * param_->win_ovl_margin();   // Overlapping windows - extension (in bp) on both ends.

  LOG_ALL("n_overlaps = %ld\n", n_overlaps);

  // Create storage for return values.
  std::vector<std::future<int>> futures;
  futures.reserve(n_overlaps);

  sampled.clear();
  sampled.resize(n_overlaps, nullptr);

//  n_overlaps = 1;
  for (int64_t oid=0; oid<n_overlaps; oid++) {
    auto& overlap = overlaps.overlaps()[oid];

    // Get at the correct query. Aid is 1-based.
    int64_t qid = overlap.Aid() - 1;
    auto& query = *(queries.get_sequences()[qid]);

    // Get at the correct target. Bid is 1-based.
    int64_t tid = overlap.Bid() - 1;
    auto& target = *(targets.get_sequences()[tid]);

    // Initialize the structure for the sampled overlap.
    sampled[oid] = createSampledOverlap();

    // Enqueue a parallel job.
    futures.emplace_back(thread_pool_->submit_task(Alignment::AlignOverlap, std::cref(query), std::cref(target), std::cref(overlap), oid, param_->window_len(), window_ext, sampled[oid]));
  }

  // Wait for threads to finish.
  for (auto& it: futures) {
      it.wait();
  }

  return 0;
}

void Racon::ConstructWindows_(const SequenceFile &targets, const Overlaps &overlaps, const std::vector<std::shared_ptr<SampledOverlap>> &sampled_overlaps, std::vector<std::vector<Window>> &windows) const {
  LOG_ALL("Constructing windows.\n");

  int64_t num_targets = targets.get_sequences().size();
  int64_t window_len = param_->window_len();
  int64_t window_ext = param_->window_len() * param_->win_ovl_margin();

  // Allocate space for random access.
  windows.clear();
  windows.resize(num_targets);
  for (int64_t i=0; i<num_targets; i++) {
    auto target = targets.get_sequences()[i];
    int64_t tlen = target->get_sequence_length();

    // int64_t num_windows = (tlen + window_len + window_ext - 1) / (window_len + window_ext);
    int64_t num_windows = (tlen + window_len- 1) / (window_len);
    windows[i].resize(num_windows, Window(i));

    for (int64_t j=0; j<windows[i].size(); j++) {
      windows[i][j].start(window_start);
      windows[i][j].end(window_end);
    }
  }

  for (int64_t i=0; i<sampled_overlaps.size(); i++) {
    AddOverlapToWindows_(targets, overlaps, sampled_overlaps[i], window_len, window_ext, windows);
  }

}

void Racon::AddOverlapToWindows_(const SequenceFile &targets, const Overlaps &overlaps, std::shared_ptr<SampledOverlap> sampled_overlap, int64_t window_len, int64_t window_ext, std::vector<std::vector<Window>> &windows) const {
  int64_t overlap_id = sampled_overlap->overlap_id();
  auto& overlap = overlaps.overlaps()[overlap_id];

  int64_t Brev = overlap.Brev();
  int64_t tid = overlap.Bid() - 1;

  auto& tw = windows[tid];		// Target windows (tw) (vector).
  auto target = targets.get_sequences()[tid];
  int64_t tlen = target->get_sequence_length();

  // Process only a part of the target sequence right around overlap positions.
  int64_t tstart = (overlap.Bstart() / window_len) * window_len;
  int64_t tend = std::min(((overlap.Bend() + window_len - 1) / window_len) * window_len, tlen);
//  printf ("Processing overlap: %ld, Astart = %ld, Aend = %ld, Alen = %ld, Bstart = %ld, Bend = %ld, Bend = %ld, Brev = %ld\n", overlap_id, overlap.Astart(), overlap.Aend(), overlap.Alen(), overlap.Bstart(), overlap.Bend(), overlap.Blen(), overlap.Brev());
  for (int64_t i=tstart; i<tend; i+=window_len) {
    // TODO: Play with the window_end coordinate - maybe it would be better to remove the -1
    // and instead just trim the sequence to 1 base less?
  	int64_t window_start = std::max(i - window_ext, (int64_t) 0);
  	int64_t window_end = std::min(i + window_len + window_ext, tlen) - 1;
  	int64_t window_id = i / window_len;

  	int64_t qstart = sampled_overlap->find(window_start);	// Found query start, or -1 if invalid.
  	int64_t qend = sampled_overlap->find(window_end);		// Found query end, or -1 if invalid.

//    printf ("  Lookup: qstart = %ld, qend = %ld, window_start = %ld, window_end = %ld\n", qstart, qend, window_start, window_end);

  	if (qstart < 0 && i == tstart) { // The beginning of a query may fall in the middle of a window.
//  	  printf ("  qstart < 0 (qstart = %ld\n", qstart);
      qstart = (overlap.Brev() == 0) ? (overlap.Astart()) : (overlap.Alen() - overlap.Aend());
      window_start = overlap.Bstart();
  	}

  	if (qend < 0 && ((i + window_len) >= tend)) {   // The end of a query may also fall in the middle of a window.
//      printf ("  qend < 0 (qend = %ld\n", qend);
      qend = (overlap.Brev() == 0) ? (overlap.Aend()) : (overlap.Alen() - overlap.Astart() - 1);
      window_end = overlap.Bend();
  	}

  	assert (qstart >= 0 && qend >= 0);

//  	printf ("  Adding to window %ld: qstart = %ld, qend = %ld, window_start = %ld, window_end = %ld\n", window_id, qstart, qend, window_start, window_end);
//  	fflush(stdout);
    tw[window_id].add(overlap_id, qstart, qend, window_start, window_end);
  }
//  printf ("\n");
//  fflush(stdout);
}

void Racon::RunAllJobs_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps, const std::vector<std::vector<Window>> &windows) const {
  LOG_ALL("Running consensus on all windows.\n");

  // Create storage for return values.
  std::vector<std::vector<std::future<int>>> futures;     // For threading.
  std::vector<std::vector<std::string>> cons_seq;         // Consensus sequences, for each window.
  std::vector<std::vector<std::string>> cons_qual;        // Quality values, for each window.
  futures.resize(windows.size());       // Reserve space for all targets
  cons_seq.resize(windows.size());      // Ditto.
  cons_qual.resize(windows.size());     // Ditto.
  for (int64_t i=0; i<windows.size(); i++) {
    futures.reserve(windows[i].size());        // Will be emplaced back.
    cons_seq[i].resize(windows[i].size());     // Needs to be known during emplacement.
    cons_qual[i].resize(windows[i].size());     // Needs to be known during emplacement.
  }

  // Emplace all jobs.
  for (int64_t i=0; i<windows.size(); i++) {
    for (int64_t j=0; j<windows[i].size(); j++) {
      futures[i].emplace_back(thread_pool_->submit_task(Racon::WindowConsensus_, std::cref(queries), std::cref(targets), std::cref(overlaps), param_, std::cref(windows[i][j]), std::ref(cons_seq[i][j]), std::ref(cons_qual[i][j])));
    }
  }

  FILE *fp_out = fopen(param_->consensus_path().c_str(), "w");
  assert(fp_out);

  for (int64_t i=0; i<futures.size(); i++) {
    for (int64_t j=0; j<futures[i].size(); j++) {
      futures[i][j].wait();
    }

    auto target = targets.get_sequences()[i];
    fprintf (fp_out, ">Consensus_%s\n", target->get_header());
    for (int64_t j=0; j<futures[i].size(); j++) {
      fprintf (fp_out, "%s", cons_seq[i][j].c_str());
    }
    fflush(fp_out);
  }

  fclose(fp_out);
}

int Racon::WindowConsensus_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps, const std::shared_ptr<Parameters> param, const Window& window, std::string& cons_seq, std::string& cons_qual) {
  std::vector<std::string> seqs, quals;
  std::vector<uint32_t> starts, ends;

  // Get the actual sequences which will be fed to SPOA.
  Racon::ExtractSequencesForSPOA_(queries, targets, overlaps, param, window, seqs, quals, starts, ends);

  // In case the coverage is too low, just pick the first sequence in the window.
  if (seqs.size() <= 2) {
    cons_seq = seqs[0];
    return 1;
  }

  auto graph = SPOA::construct_partial_order_graph(seqs, quals, starts, ends,
                                             SPOA::AlignmentParams(param->match(), param->mismatch(),
                                             param->gap_open(), param->gap_ext(), (SPOA::AlignmentType) param->aln_type()));

  std::vector<uint32_t> coverages;
  cons_seq = graph->generate_consensus(coverages);

  // Unfortunately, POA is bad when there are errors, such as long insertions, at
  // the end of a sequence. The consensus walk will also have those overhang
  // nodes, which then need to be trimmed heuristically.
  int32_t start_offset = 0, end_offset = cons_seq.size() - 1;
  for (;start_offset<cons_seq.size(); start_offset++) {
    if (coverages[start_offset] >= ((seqs.size() - 1) / 2)) { break; }
  }
  for (; end_offset >= 0; end_offset--) {
    if (coverages[start_offset] >= ((seqs.size() - 1) / 2)) { break; }
  }

  cons_seq = cons_seq.substr(start_offset, (end_offset - start_offset + 1));

  return 0;
}

void Racon::ReverseInPlace_(std::string &seq) {
  int64_t len = seq.size();
  for (int64_t i=0; i<len/2; i++) {
    std::swap(seq[i], seq[len-i-1]);
  }
}

void Racon::ReverseComplementInPlace_(std::string &seq) {
  ReverseInPlace_(seq);

  int64_t len = seq.size();
  for (int64_t i=0; i<len; i++) {
    seq[i] = kBaseComplement[seq[i]];
  }
}

void Racon::ExtractSequencesForSPOA_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps, const std::shared_ptr<Parameters> param, const Window& window, std::vector<std::string>& seqs, std::vector<std::string>& quals, std::vector<uint32_t> &starts, std::vector<uint32_t> &ends) {
  seqs.clear();
  quals.clear();
  starts.clear();
  ends.clear();

  int64_t window_start = 0;
  int64_t window_end = 0;
  window.start
  int64_t tid = window.target_id();
  seqs.emplace_back(targets.get_sequences()[tid]->GetSequenceAsString(window_start, window_end));
//  quals.push_back(targets.get_sequences()[tid]->GetQualityAsString(window_start, window_end));
  std::string dummy_quals((window_end - window_start + 1), '!');
  quals.emplace_back(dummy_quals);

  auto entries = window.entries();
  for (int64_t i=0; i<entries.size(); i++) {
    auto& entry = entries[i];
    auto& overlap = overlaps.overlaps()[entry.overlap_id()];
    auto query = queries.get_sequences()[overlap.Aid() - 1];
    auto target = targets.get_sequences()[overlap.Bid() - 1];

    if (overlap.Brev() == 0) {  // The query is forward.
      seqs.emplace_back(query->GetSequenceAsString(entry.query().start, entry.query().end));

    } else {                    // The query is reverse-complemented.
      int64_t start = overlap.Alen() - entry.query().end;
      int64_t end = overlap.Alen() - entry.query().start - 1;
      std::string seq = query->GetSequenceAsString(start, end);
      ReverseComplementInPlace_(seq);
      seqs.emplace_back(seq);

      std::string qual = query->GetQualityAsString(start, end);
      ReverseInPlace_(qual);
      quals.emplace_back(qual);
    }
  }
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
