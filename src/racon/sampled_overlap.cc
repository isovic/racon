/*
 * sampled_overlap.cc
 *
 *  Created on: Jan 28, 2017
 *      Author: isovic
 */

#include "sampled_overlap.h"
#include <deque>

namespace is {

std::shared_ptr<SampledOverlap> createSampledOverlap(
    const Overlap& overlap, int64_t overlap_id, const std::vector<uint8_t>& alignment,
    int64_t sample_step, int64_t extension) {
  return std::shared_ptr<SampledOverlap>(
      new SampledOverlap(overlap, overlap_id, alignment, sample_step, extension));
}

std::shared_ptr<SampledOverlap> createSampledOverlap() {
  return std::shared_ptr<SampledOverlap>(
      new SampledOverlap());
}

SampledOverlap::SampledOverlap() : overlap_id_(-1) {
}

SampledOverlap::SampledOverlap(const Overlap& overlap, int64_t overlap_id,
                               const std::vector<uint8_t>& alignment,
                               int64_t sample_step, int64_t extension) {
  PopulatePos_(overlap, overlap_id, alignment, sample_step, extension);
}

SampledOverlap::~SampledOverlap() {

}

void SampledOverlap::set(const Overlap& overlap, int64_t overlap_id,
                         const std::vector<uint8_t>& alignment,
                         int64_t sample_step, int64_t extension) {
  PopulatePos_(overlap, overlap_id, alignment, sample_step, extension);
}

int64_t SampledOverlap::find(int64_t target_pos) const {
  auto it = pos_.find(target_pos);
  if (it == pos_.end()) {
    return -1;
  }
  return it->second;
}

void SampledOverlap::Verbose(std::ostream& os) {
  int64_t i = 0;
  os << "This overlap has " << pos_.size() << " sampled points:\n";
  for (auto& it: pos_) {
    os << "[" << i << "] " << it.first << " -> " << it.second << "\n";
    i += 1;
  }

}

void SampledOverlap::PopulatePos_(const Overlap& overlap, int64_t overlap_id,
                                  const std::vector<uint8_t>& alignment,
                                  int64_t sample_step, int64_t extension) {

  overlap_id_ = overlap_id;

  int64_t qpos = (overlap.Brev() == 0) ? (overlap.Astart()) : (overlap.Alen() - overlap.Aend());
  int64_t rpos = overlap.Bstart();
  int64_t rpos_end = overlap.Bend();
  // Find the first position of a window at >= rpos.
  int64_t window_pos = ((rpos + sample_step - 1) / sample_step) * sample_step;

  // Initialize a queue with positions of interest.
  // This includes the actual window boundaries, as well as the
  // locations of extended windows in the overlapping windows mode.
  std::deque<int64_t> pos_to_store;
  for (int64_t i = window_pos; i < rpos_end; i += sample_step) {
    if (extension > 0) {
      // Store the overlapping window extension coordinates.
      if ((i - extension - 1) >= 0) { pos_to_store.push_back(i - extension - 1); }
      if ((i - extension) >= 0) { pos_to_store.push_back(i - extension); }
      if ((i - extension + 1) < rpos_end) { pos_to_store.push_back(i - extension + 1); }
    }

    // Store the actual window boundary.
    if ((i - 1) >= 0) { pos_to_store.push_back(i - 1); }
    if (i < rpos_end) { pos_to_store.push_back(i); }
    if ((i + 1) < rpos_end) { pos_to_store.push_back(i + 1); }

    if (extension > 0) {
      // Store the overlapping window extension coordinates.
      if ((i + extension - 1) < rpos_end) { pos_to_store.push_back(i + extension - 1); }
      if ((i + extension) < rpos_end) { pos_to_store.push_back(i + extension); }
      if ((i + extension + 1) < rpos_end) { pos_to_store.push_back(i + extension + 1); }
    }
  }

  pos_.clear();

  for (auto& op : alignment) {
    if (pos_to_store.size() == 0) { break; }
    // Hash the positions.
    if (op == EDLIB_EDOP_MATCH || op == EDLIB_EDOP_MISMATCH
        || op == EDLIB_EDOP_DELETE) {
      if (rpos == pos_to_store.front()) {
        pos_[rpos] = qpos;
        pos_to_store.pop_front();
      }
    }

    // Move down the alignment.
    if (op == EDLIB_EDOP_MATCH || op == EDLIB_EDOP_MISMATCH) {
      qpos += 1;
      rpos += 1;
    } else if (op == EDLIB_EDOP_INSERT) {
      qpos += 1;
    } else if (op == EDLIB_EDOP_DELETE) {
      rpos += 1;
    }
  }

}

} /* namespace is */
