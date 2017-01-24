/*
 * overlaps.h
 *
 *  Created on: Jan 24, 2017
 *      Author: isovic
 */

#ifndef SRC_REFACTOR_OVERLAPS_H_
#define SRC_REFACTOR_OVERLAPS_H_

#include <string>
#include <vector>

namespace is {

class OverlapLine {
 public:
  OverlapLine(const std::string &line);
  ~OverlapLine();

 private:
  int64_t Aid_, Bid_;
  std::string Aname_, Bname_;
  double perc_err_;
  int64_t shared_minmers_;
  int64_t Arev_, Astart_, Aend_, Alen_;     // start is zero-based, end points to a position right after the last inclusive base.
  int64_t Brev_, Bstart_, Bend_, Blen_;
};

class Overlaps {
 public:
  Overlaps();
  ~Overlaps();
};

} /* namespace is */

#endif /* SRC_REFACTOR_OVERLAPS_H_ */
