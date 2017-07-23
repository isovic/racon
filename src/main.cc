/*
 * main.cc
 *
 *  Created on: Jan 24, 2017
 *      Author: isovic
 */

#ifndef RUN_ALL_TESTS_

#include "racon/parameters.h"
#include "racon/racon.h"

int main(int argc, char* argv[]) {
  std::shared_ptr<is::Parameters> p = std::move(is::createParameters(argc, argv));
  auto r = is::createRacon(p);
  r->CreateConsensus();

  return 0;
}

#endif
