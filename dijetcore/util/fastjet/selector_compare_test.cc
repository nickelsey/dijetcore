#include "gtest/gtest.h"

#include "dijetcore/util/fastjet/selector_compare.h"

#include "fastjet/Selector.hh"

TEST(SelectorCompare, SelectorCompare) {
  
  fastjet::Selector ptminsel1 = fastjet::SelectorPtMin(0.2);
  fastjet::Selector ptminsel2 = fastjet::SelectorPtMin(2.0);
  fastjet::Selector ptmaxsel1 = fastjet::SelectorPtMax(0.2);
  fastjet::Selector ptmaxsel2 = fastjet::SelectorPtMax(2.0);
  fastjet::Selector compound_sel = ptminsel1 && ptmaxsel2;
  
  EXPECT_TRUE(dijetcore::SelectorPtMinLessThan(ptminsel1, ptminsel2));
  EXPECT_FALSE(dijetcore::SelectorPtMinGreaterThan(ptminsel1, ptminsel2));
  EXPECT_FALSE(dijetcore::SelectorPtMinLessThan(ptminsel2, ptminsel1));
  EXPECT_TRUE(dijetcore::SelectorPtMinGreaterThan(ptminsel2, ptminsel1));
  
  EXPECT_TRUE(dijetcore::SelectorPtMaxLessThan(ptmaxsel1, ptmaxsel2));
  EXPECT_FALSE(dijetcore::SelectorPtMaxGreaterThan(ptmaxsel1, ptmaxsel2));
  EXPECT_FALSE(dijetcore::SelectorPtMaxLessThan(ptmaxsel2, ptmaxsel1));
  EXPECT_TRUE(dijetcore::SelectorPtMaxGreaterThan(ptmaxsel2, ptmaxsel1));
  
  EXPECT_FALSE(dijetcore::SelectorPtMinLessThan(ptminsel2, ptminsel1));
  EXPECT_TRUE(dijetcore::SelectorPtMinGreaterThan(ptminsel2, ptminsel1));
  EXPECT_TRUE(dijetcore::SelectorPtMinLessThan(ptminsel1, ptminsel2));
  EXPECT_FALSE(dijetcore::SelectorPtMinGreaterThan(ptminsel1, ptminsel2));
  
  EXPECT_FALSE(dijetcore::SelectorPtMaxLessThan(ptmaxsel2, ptmaxsel1));
  EXPECT_TRUE(dijetcore::SelectorPtMaxGreaterThan(ptmaxsel2, ptmaxsel1));
  EXPECT_TRUE(dijetcore::SelectorPtMaxLessThan(ptmaxsel1, ptmaxsel2));
  EXPECT_FALSE(dijetcore::SelectorPtMaxGreaterThan(ptmaxsel1, ptmaxsel2));
  
  EXPECT_FALSE(dijetcore::SelectorPtMinGreaterThan(compound_sel, ptminsel2));
  EXPECT_TRUE(dijetcore::SelectorPtMaxGreaterThan(compound_sel, ptmaxsel1));
  EXPECT_TRUE(dijetcore::SelectorPtMinLessThan(compound_sel, ptminsel2));
  EXPECT_FALSE(dijetcore::SelectorPtMaxLessThan(compound_sel, ptmaxsel1));
  
  EXPECT_TRUE(dijetcore::SelectorPtMaxLessEqualThan(compound_sel, ptmaxsel2));
  EXPECT_TRUE(dijetcore::SelectorPtMinLessEqualThan(compound_sel, ptminsel1));
  EXPECT_FALSE(dijetcore::SelectorPtMaxLessThan(compound_sel, ptmaxsel2));
  EXPECT_FALSE(dijetcore::SelectorPtMinLessThan(compound_sel, ptminsel1));
}
