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

namespace is {

std::shared_ptr<Racon> createRacon(const Parameters& param) {
  return std::shared_ptr<Racon>(new Racon(param));
}

Racon::~Racon() {
}

Racon::Racon(const Parameters& param) {
}

Racon::Racon(const SequenceFile& reads, const SequenceFile& targets) {
}

void Racon::PopulateJobsConsensus_(const SequenceFile &refs, int64_t win_len, JobQueue &jobs) const {
  assert(win_len > 0);
  for (int64_t i=0; i<((int64_t) refs.get_sequences().size()); i++) {
  	auto s = refs.get_sequences()[i];
  	int64_t s_len = (int64_t) s->get_sequence_length();
    for (int64_t j=0; j<s_len; j+=win_len) {
      int64_t start = j;
      int64_t end = std::min(start + win_len, s_len);
   	  jobs.push_back(createJob(i, start, end, win_len));
    }
  }
}

void Racon::PopulateJobsErc_(const SequenceFile &refs, int64_t win_len, JobQueue &jobs) const {
  assert(win_len > 0);
  for (int64_t i=0; i<((int64_t) refs.get_sequences().size()); i++) {
    auto s = refs.get_sequences()[i];
    int64_t s_len = (int64_t) s->get_sequence_length();
    jobs.push_back(createJob(i, 0, s_len, win_len));
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

void Racon::RunFromOverlaps_(const Parameters& param) {
  std::string overlaps_file = parameters.aln_path;
  if (parameters.is_paf == true) { LOG_ALL("Using PAF for input alignments. (%s)\n", overlaps_file.c_str());}
  else { LOG_ALL("Using MHAP for input alignments. (%s)\n", overlaps_file.c_str()); }

  std::vector<OverlapLine> overlaps, overlaps_filtered, overlaps_final;

  LOG_ALL("Loading reads.\n");
  SequenceFile seqs_reads(SEQ_FORMAT_AUTO, parameters.reads_path);

  // Sanity check to see if the reads have quality values.
  if (seqs_reads.HasQV() == false) {
    fprintf (stderr, "ERROR: Reads are not specified in a format which contains quality information. Exiting.\n");
    exit(1);
  }

  // Hash the read sequences by their name.
  LOG_ALL("Hashing qnames.\n");
  std::map<std::string, int64_t> qname_to_ids;
  std::map<std::string, int64_t> rname_to_ids;
  HashQnames(seqs_gfa, rname_to_ids);
  HashQnames(seqs_reads, qname_to_ids);

  LOG_ALL("Parsing the overlaps file.\n");
  if (overlaps_file == "-") { LOG_ALL("Stdin will be used to load the overlap lines.\n"); }
  OverlapFormat overlap_format = (parameters.is_paf) ? kOverlapFormatPAF : kOverlapFormatMHAP;
  if (parameters.do_erc == false) {
    LOG_ALL("Unique overlaps will be filtered on the fly.\n");
    ParseUniqueAndFilterErrors(overlaps_file, overlap_format, qname_to_ids, rname_to_ids, parameters.error_rate, overlaps_final);
  } else {
    ParseAndFilterErrors(overlaps_file, overlap_format, qname_to_ids, rname_to_ids, parameters.error_rate, overlaps_final);
  }

  std::sort(overlaps_final.begin(), overlaps_final.end(), [](const OverlapLine &a, const OverlapLine &b){ return (a.Bid < b.Bid); } );
  if (parameters.do_sparse == false || parameters.do_erc) {
    LOG_ALL("Overlaps will be fully aligned.\n");
    ConsensusFromOverlaps(parameters, seqs_gfa, seqs_reads, qname_to_ids, rname_to_ids, overlaps_final);
  }
}

void Racon::RunFromAlignments_() {
}

void Racon::VerboseJobs_() {
  int64_t count;
//  for (auto& job: jobs_) {
  for (auto it = jobs_.begin(); it != jobs_.end(); it++) {
    auto& job = *it;
    std::cout << "[" << count << "] seq_id = " << job->seq_id() << ", start = " << job->start() << ", end = " << job->end() << ", win_len = " << job->win_len() << std::endl;
    count++;
  }
}

}
