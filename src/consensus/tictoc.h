/*
 * tictoc.h
 *
 *  Created on: Nov 29, 2016
 *      Author: isovic
 */

#ifndef SRC_CONSENSUS_TICTOC_H_
#define SRC_CONSENSUS_TICTOC_H_

#include <time.h>

namespace racon {
class TicToc {
 public:
  TicToc();
  ~TicToc();

  void start();
  void stop();
  double get_secs();

 private:
  clock_t start_;
  clock_t end_;
};

}

#endif /* SRC_CONSENSUS_TICTOC_H_ */
