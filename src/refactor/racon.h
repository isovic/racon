/*
 * racon.h
 *
 *  Created on: January 18, 2017
 *      Author: Ivan Sovic
 */

#ifndef RACON_H_
#define RACON_H_

#include <stdint.h>
#include <map>
#include <memory>
#include <string>
#include "parameters.h"
#include "sequences/sequence_file.h"

namespace is {

class Racon;

std::shared_ptr<Racon> createRacon(const Parameters& param);

class Racon {
 public:
  ~Racon();

  friend std::shared_ptr<Racon> createRacon(const Parameters& param);

 private:
  typedef std::map<std::string, int64_t> MapId;

  Racon(const Racon&) = delete;
  const Racon& operator=(const Racon&) = delete;
  Racon(const Parameters& param);
  Racon(const SequenceFile& reads, const SequenceFile& targets);

  void HashNames_(const SequenceFile &seqs, MapId &id) const;
  MapId read_id_;
  MapId target_id_;
};

}

#endif
