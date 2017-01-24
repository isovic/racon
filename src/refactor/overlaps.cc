/*
 * overlaps.cc
 *
 *  Created on: Jan 24, 2017
 *      Author: isovic
 */

#include "overlaps.h"
#include <sstream>

namespace is {

OverlapLine::OverlapLine(const std::string &line) {
  std::istringstream iss(line);
  if (!(iss >> Aid_ >> Bid_ >> perc_err_ >> shared_minmers_ >> Arev_ >> Astart_ >> Aend_ >> Alen_ >> Brev_ >> Bstart_ >> Bend_ >> Blen_)) {
//    ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Overlaps are not formatted in the MHAP format. Exiting.");
  }
//  Aname_ = FormatString("%ld", Aid_);
//  Bname_ = FormatString("%ld", Bid_);
}

Overlaps::Overlaps() {
}

Overlaps::~Overlaps() {
  // TODO Auto-generated destructor stub
}

} /* namespace is */
