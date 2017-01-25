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

#include "types.h"
// #include "job.h"
#include "parameters.h"
#include "sequences/sequence_file.h"
#include "overlaps.h"

namespace is {

class Racon;

std::shared_ptr<Racon> createRacon(const Parameters& param);

class Racon {
 public:
  ~Racon();

  friend std::shared_ptr<Racon> createRacon(const Parameters& param);

 private:
//  typedef std::deque<std::shared_ptr<Job>> JobQueue;

  Racon(const Racon&) = delete;
  const Racon& operator=(const Racon&) = delete;
  Racon(const Parameters& param);
  Racon(const SequenceFile& reads, const SequenceFile& targets);

  void RunFromOverlaps_(const Parameters& param);
  void RunFromAlignments_(const Parameters& param);

  /** A helper function to fill a map in which the key is a sequence
  	* name, and the value is the ordinal number of the sequence in
  	* the SequenceFile.
  */
  void HashNames_(const SequenceFile &seqs, MapId &id) const;

  /** Takes a sorted vector of all overlaps loaded from a file. The overlaps
   * are sorted ascending by their Bid (target ID). This means that all
   * overlaps for a particular contig will be contiguous on the list.
   * This function finds the range in the vector for each separate contig
   * and stores it into a map.
   */
  int FindContigOverlaps_(const Overlaps &sorted_overlaps, MapOverlapRange &contig_overlaps) const;

  /** Populates the job queue for consensus-type parallelization.
  	* This means that each window of the contig is presented as a
  	* single job.
  	* @refs Reference sequences from which to create jobs.
  	* @window_len Length of the window for processing. Must be >= 0.
  	* @jobs The deque to which to append the jobs.
  */
//  void PopulateJobsConsensus_(const SequenceFile &refs, int64_t win_len, JobQueue &jobs) const;

  /** Populates the job queue for error-correction type parallelization.
   * This means that each sequence will be a single threadpool job.
   * @refs Reference sequences from which to create jobs.
   * @window_len Length of the window for processing. Must be >= 0.
   * @jobs The deque to which to append the jobs.
   */
//  void PopulateJobsErc_(const SequenceFile &refs, int64_t win_len, JobQueue &jobs) const;

//  MapId query_id_;
//  MapId target_id_;
//  JobQueue jobs_;
};

} /* namespace is */

#endif
