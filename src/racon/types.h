/*
 * types.h
 *
 *  Created on: Jan 25, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_TYPES_H_
#define SRC_REFACTOR_TYPES_H_

#include <stdint.h>
#include <map>
#include <string>

namespace is {

class OverlapFormat;

// A generic range structure. Beats unnamed pairs/tuples.
typedef struct {
  int64_t start = 0, end = 0;
} Range;

typedef std::map<std::string, int64_t> MapId;
typedef std::map<int64_t, Range> MapOverlapRange;
// typedef std::map<std::string, const SingleSequence*> MapSeqs;

enum class ConsensusType {
  Contig,
  Read
};

class OverlapFormat {
 public:
  OverlapFormat(const OverlapFormat& op) { // Public copy constructor to allow setting
    x_ = op.x_;                                     // Only from static factory functions.
  }                                                 // Only assignment of int values should be hidden

  OverlapFormat& operator=(const OverlapFormat& op) { // Public operator= for assignments.
    if (&op == this) { return *this; };
    x_ = op.x_;
    return *this;
  }

  static OverlapFormat Paf() { return OverlapFormat(kPaf); }    // Factory functions.
  static OverlapFormat Mhap() { return OverlapFormat(kMhap); }
  static OverlapFormat Sam() { return OverlapFormat(kSam); }

  bool isPaf() const { return x_ == kPaf; }     // Testing values.
  bool isMhap() const { return x_ == kMhap; }
  bool isSam() const { return x_ == kSam; }

 private:
  OverlapFormat(int32_t x): x_(x) { }  // Don't allow direct construction. Easy to use wrong.
  int32_t x_;

  enum OvlFormatTypes {
    kPaf = 0,
    kMhap = 1,
    kSam = 2
  };
};

}




#endif /* SRC_REFACTOR_TYPES_H_ */
