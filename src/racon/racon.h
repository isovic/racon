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

std::unique_ptr<Racon> createRacon(const std::shared_ptr<Parameters> param);

class Racon {
 public:
  ~Racon();

  friend std::unique_ptr<Racon> createRacon(const std::shared_ptr<Parameters> param);
  void CreateConsensus();

 private:
  Racon(const Racon&) = delete;
  const Racon& operator=(const Racon&) = delete;
  Racon(std::shared_ptr<Parameters> param);

  void RunFromOverlaps_();
  void RunFromAlignments_();

  int AlignAndSampleOverlaps_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps,
                              std::vector<std::shared_ptr<SampledOverlap>> &sampled);

  void ConstructWindows_(const SequenceFile &targets, const Overlaps &overlaps, const std::vector<std::shared_ptr<SampledOverlap>> &sampled_overlaps,
                         std::vector<std::vector<Window>> &windows) const;

  void AddSampledOverlapToWindows_(const SequenceFile &targets, const Overlaps &overlaps, std::shared_ptr<SampledOverlap> sampled_overlap,
                            int64_t window_len, int64_t window_ext, std::vector<std::vector<Window>> &windows) const;

  void RunAllJobs_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps,
                   const std::vector<std::vector<Window>> &windows, ConsensusType cons_type) const;

  static int WindowConsensus_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps,
                              const std::shared_ptr<Parameters> param, const std::vector<Window>& windows, std::vector<std::string>& cons_seqs,
                              std::vector<std::string>& cons_quals, int64_t starting_window, int64_t ending_window);

  static void ExtractSequencesForSPOA_(const SequenceFile &queries, const SequenceFile &targets, const Overlaps &overlaps,
                                       const std::shared_ptr<Parameters> param, const Window& window, std::vector<std::string>& seqs,
                                       std::vector<std::string>& quals, std::vector<uint32_t> &starts, std::vector<uint32_t> &ends);

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

  static void ReverseInPlace_(std::string &seq);
  static void ReverseComplementInPlace_(std::string &seq);

  const std::shared_ptr<Parameters> param_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  std::vector<std::vector<Window>> windows_;
};

} /* namespace is */

#endif
