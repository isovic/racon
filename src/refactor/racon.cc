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

std::shared_ptr<Racon> createRacon(const Parameters& param) {
  return std::shared_ptr<Racon>(new Racon(param));
}

Racon::~Racon() {
}

Racon::Racon(const Parameters& param) {
  if (param.overlap_format().isPaf() || param.overlap_format().isMhap()) {
    RunFromOverlaps_(param);
  }
}

Racon::Racon(const SequenceFile& reads, const SequenceFile& targets) {
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

void Racon::RunFromOverlaps_(const Parameters& param) {
  // Parse the backbone.
  SequenceFile targets(SEQ_FORMAT_AUTO, param.contigs_path());

  // Parse the reads.
  SequenceFile queries(SEQ_FORMAT_AUTO, param.reads_path());
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
          (param.overlap_format().isPaf()) ? "PAF" : "MHAP", param.aln_path().c_str())
  LOG_ALL("Started parsing the overlaps file.\n");

  Overlaps overlaps(param.aln_path(), param.overlap_format(), query_id, target_id, param.error_rate(), param.do_erc());
  overlaps.SortByTargetId();

//  if (parameters.do_sparse == false || parameters.do_erc) {
//    LOG_ALL("Overlaps will be fully aligned.\n");
//    ConsensusFromOverlaps(parameters, seqs_gfa, seqs_reads, qname_to_ids, rname_to_ids, overlaps_final);
//  }
}

void Racon::RunFromAlignments_(const Parameters& param) {
}

//void Racon::VerboseJobs_() {
//  int64_t count;
////  for (auto& job: jobs_) {
//  for (auto it = jobs_.begin(); it != jobs_.end(); it++) {
//    auto& job = *it;
//    std::cout << "[" << count << "] seq_id = " << job->seq_id() << ", start = " << job->start() << ", end = " << job->end() << ", win_len = " << job->win_len() << std::endl;
//    count++;
//  }
//}

}
