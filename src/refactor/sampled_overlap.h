/*
 * sampled_overlap.h
 *
 *  Created on: Jan 28, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_SAMPLED_OVERLAP_H_
#define SRC_REFACTOR_SAMPLED_OVERLAP_H_

#include <memory>
#include <unordered_map>

#include "sequences/sequence_alignment.h"
#include "libs/edlib.h"
#include "libs/edlibcigar.h"
#include "overlaps.h"

namespace is {

class SampledOverlap;

std::shared_ptr<SampledOverlap> createSampledOverlap();
std::shared_ptr<SampledOverlap> createSampledOverlap(const Overlap& overlap, int64_t overlap_id, const std::vector<uint8_t>& alignment, int64_t sample_step, int64_t extension);

class SampledOverlap {
 public:
  ~SampledOverlap();
  friend std::shared_ptr<SampledOverlap> createSampledOverlap();
  friend std::shared_ptr<SampledOverlap> createSampledOverlap(const Overlap& overlap, int64_t overlap_id, const std::vector<uint8_t>& alignment, int64_t sample_step, int64_t extension);

  /** Initializes the object in case it was built using a
   * plain constructor.
   */
  void set(const Overlap& overlap, int64_t overlap_id, const std::vector<uint8_t>& alignment, int64_t sample_step, int64_t extension);

  /** Find a query position for a given target position,
   * or return an invalid value (-1).
   */
  int64_t find(int64_t target_pos) const;

  const std::unordered_map<int64_t, int64_t>& pos() const {
    return pos_;
  }

  int64_t overlap_id() const {
    return overlap_id_;
  }

 private:
  SampledOverlap();
  SampledOverlap(const Overlap& overlap, int64_t overlap_id, const std::vector<uint8_t>& alignment, int64_t sample_step, int64_t extension);
  SampledOverlap(SampledOverlap&) = delete;
  SampledOverlap& operator=(SampledOverlap&) = delete;

  void PopulatePos_(const Overlap& overlap, int64_t overlap_id, const std::vector<uint8_t>& alignment, int64_t sample_step, int64_t extension);

  std::unordered_map<int64_t, int64_t> pos_;     // Key: position on target, value: position on query.
  int64_t overlap_id_;
};

} /* namespace is */

#endif /* SRC_REFACTOR_SAMPLED_OVERLAP_H_ */
