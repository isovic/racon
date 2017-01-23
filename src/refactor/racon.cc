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
