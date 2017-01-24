/*
 * racon.cc
 *
 *  Created on: January 18, 2017
 *      Author: Ivan Sovic
 */

#include "racon.h"
#include "utility/utility_general.h"

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

int Racon::PopulateJobsConsensus_(const SequenceFile &refs, int64_t window_len, std::deque<Job> &jobs) const {
  if (window_len <= 0) {
  	return 1;
  }

  for (int64_t i=0; i<((int64_t) refs.get_sequences().size()); i++) {
  	auto s = refs.get_sequences()[i];
  	int64_t s_len = (int64_t) s->get_sequence_length();
    for (size_t j=0; j<s_len; j+=window_len) {
      int64_t start = j;
      int64_t end = std::min(start + window_len, s_len);
   	  jobs.push_back(Job(i, start, end));
    }
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

}
