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
#include <unordered_map>

#include "types.h"
// #include "job.h"
#include "parameters.h"
#include "sequences/sequence_file.h"
#include "overlaps.h"
#include "thread_pool.hpp"
#include "sampled_overlap.h"
#include "window.h"

namespace is {

class Racon;

std::shared_ptr<Racon> createRacon(const std::shared_ptr<Parameters> param);

class Racon {
 public:
  ~Racon();

  friend std::shared_ptr<Racon> createRacon(const std::shared_ptr<Parameters> param);
  void CreateConsensus();

 private:
  const std::shared_ptr<Parameters> param_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  std::vector<std::vector<Window>> windows_;

  Racon(const Racon&) = delete;
  const Racon& operator=(const Racon&) = delete;
  Racon(std::shared_ptr<Parameters> param);

  void RunFromOverlaps_();
  void RunFromAlignments_();

  int AlignAndSampleOverlaps_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps, std::vector<std::shared_ptr<SampledOverlap>> &sampled);

  void ConstructWindows_(const SequenceFile &targets, const Overlaps &overlaps, const std::vector<std::shared_ptr<SampledOverlap>> &sampled_overlaps, std::vector<std::vector<Window>> &windows) const;
  void AddOVerlapToWindows_(const SequenceFile &targets, const Overlaps &overlaps, std::shared_ptr<SampledOverlap> sampled_overlap, std::vector<std::vector<Window>> &windows) const;

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

};

} /* namespace is */

#endif
