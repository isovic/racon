/*
 * pileup.h
 *
 *  Created on: May 16, 2016
 *      Author: isovic
 */

#ifndef SRC_CONSENSUS_PILEUP_H_
#define SRC_CONSENSUS_PILEUP_H_

#include <string>
#include <vector>
#include <map>
#include "sequences/single_sequence.h"
#include "sequences/sequence_file.h"

enum EventType {
  kEventBase = 0,
  kEventInsertion = 1,
  kEventDeletion = 2,
  kEventUnknown = 3
};

typedef struct {
  EventType t = kEventBase;     // Event type: either a base, an insertion or a deletion.
  std::string e = "";           // Event bases.
  std::string q = "";           // Qualities for each base in the event.
  int32_t qid = -1;             // Read ID that the base originated from.
} Event;

typedef struct {
  std::vector<Event> events;
  int32_t coverage = 0;
  int32_t base_coverage = 0;
  int32_t del_coverage = 0;
  char ref_base = '?';
} RefBase;

class Pileup {
 public:
  Pileup(const SingleSequence *ref, std::vector<const SingleSequence *> &ref_alns);
  ~Pileup();

  void AddAlignment(const SingleSequence* ref, const SingleSequence *seq, int64_t qid);
  void GenerateConsensus(int32_t cov_threshold, std::string &cons);
  int MajorityVoteFromMSA(std::vector<std::string> &msa, std::string &consensus);

  void Verbose(FILE *fp);

 private:
  std::vector<RefBase> bases_;      // The actual pileup for every base of the reference.
  const SingleSequence* ref_;
};



#endif /* SRC_CONSENSUS_PILEUP_H_ */
