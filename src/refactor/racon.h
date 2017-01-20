#ifndef RACON_H_
#define RACON_H_

#include <stdint.h>
#include <map>
#include <memory>
#include <string>

namespace is {

class Racon {
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
