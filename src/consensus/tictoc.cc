/*
 * tictoc.cc
 *
 *  Created on: Nov 29, 2016
 *      Author: isovic
 */

#include "tictoc.h"

namespace racon {
TicToc::TicToc() : start_(0), end_(0) {
}

TicToc::~TicToc() {
}

void TicToc::start() {
  start_ = clock();
}

void TicToc::stop() {
  end_ = clock();
}

double TicToc::get_secs() {
  double elapsed_secs = double(end_ - start_) / CLOCKS_PER_SEC;
  return elapsed_secs;
}

}
