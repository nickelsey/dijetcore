#include "gtest/gtest.h"

#include "dijetcore/lib/containers.h"

#include <vector>
#include <limits>

using std::string;

TEST(FindFirst, Find) {
  std::vector<int> c{1, 6, 3, 6, 3, 8};
  size_t idx = dijetcore::FindFirst(c, 6);
  EXPECT_EQ(idx, 1);
}

TEST(FindFirst, NotPresent) {
  std::vector<int> c{1, 6, 3, 6, 3, 8};
  size_t idx = dijetcore::FindFirst(c, 10);
  EXPECT_EQ(idx, std::numeric_limits<size_t>::max());
}