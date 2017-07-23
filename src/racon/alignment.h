/*
 * alignment.h
 *
 *  Created on: Jan 26, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_ALIGNMENT_H_
#define SRC_REFACTOR_ALIGNMENT_H_

#include <stdint.h>
#include <memory>
#include "libs/edlib.h"
#include "libs/edlibcigar.h"
#include "sequences/sequence_file.h"
#include "overlaps.h"
#include "sampled_overlap.h"

namespace is {

class Alignment {
 public:
  static int AlignOverlap(const SingleSequence& query, const SingleSequence& target, const Overlap& overlap, int64_t overlap_id, int64_t win_size, int64_t win_ext, std::shared_ptr<SampledOverlap> sampled);

  static int AlignNW(const int8_t *q, int64_t qlen, const int8_t *t,
                     int64_t tlen, int64_t* start, int64_t *end,
                     int64_t *eddist, std::vector<unsigned char> &alignment);
  static int32_t CalcTargetLen(const unsigned char *aln, int32_t len);

 private:
  Alignment() = delete;
  ~Alignment() = delete;
};

} /* namespace is */

#endif /* SRC_INTERVALTREE_ALIGNMENT_H_ */
