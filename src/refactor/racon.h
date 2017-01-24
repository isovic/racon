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
#include <deque>
#include "parameters.h"
#include "sequences/sequence_file.h"

namespace is {

class Job {
 public:
  explicit Job(int64_t ref_id, int64_t start, int64_t end) : ref_id_(ref_id), start_(start), end_(end) { }
  ~Job() { };

 private:
  int64_t ref_id_;
  int64_t start_;
  int64_t end_;
};

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

  /** A helper function to fill a map in which the key is a sequence
  	* name, and the value is the ordinal number of the sequence in
  	* the SequenceFile.
  */
  void HashNames_(const SequenceFile &seqs, MapId &id) const;

  /** Populates the job queue for consensus-type parallelization.
  	* This means that each window of the contig is presented as a
  	* single job.
  	* @refs Reference sequences from which to create jobs.
  	* @jobs The deque to which to append the jobs.
  */
  int PopulateJobsConsensus_(const SequenceFile &refs, int64_t window_len, std::deque<Job> &jobs) const;

  MapId read_id_;
  MapId target_id_;
  std::deque<Job> jobs_;
};

}

#endif
