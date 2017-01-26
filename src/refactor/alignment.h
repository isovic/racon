/*
 * alignment.h
 *
 *  Created on: Jan 26, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_ALIGNMENT_H_
#define SRC_REFACTOR_ALIGNMENT_H_

#include <stdint.h>
#include "libs/edlib.h"
#include "libs/edlibcigar.h"
#include "sequences/sequence_file.h"

namespace is {

class Alignment {
 public:
//  int AlignOverlaps(const SequenceFile &refs, const SequenceFile &reads,
//                    const std::vector<OldOverlapLine> &overlaps,
//                    int32_t num_threads, SequenceFile &aligned,
//                    bool verbose_debug);

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
