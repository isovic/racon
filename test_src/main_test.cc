#ifdef RUN_ALL_TESTS_

#include <stdio.h>
#include <stdint.h>
#include "gtest/gtest.h"
#include "test_alignment.h"

GTEST_API_ int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#endif
