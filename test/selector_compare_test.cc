#include "jet_analysis/util/selector_compare.hh"

#include "fastjet/Selector.hh"

#include <iostream>
#include <string>
#include <sstream>

int main() {
  
  fastjet::Selector ptminsel1 = fastjet::SelectorPtMin(0.2);
  fastjet::Selector ptminsel2 = fastjet::SelectorPtMin(2.0);
  fastjet::Selector ptmaxsel1 = fastjet::SelectorPtMax(0.2);
  fastjet::Selector ptmaxsel2 = fastjet::SelectorPtMax(2.0);
  
  if (!SelectorPtMinLessThan(ptminsel1, ptminsel2) ||
      SelectorPtMinGreaterThan(ptminsel1, ptminsel2))
    return 1;
  if (!SelectorPtMaxLessThan(ptmaxsel1, ptmaxsel2) ||
      SelectorPtMaxGreaterThan(ptmaxsel1, ptmaxsel2))
    return 1;
  
  if (SelectorPtMinLessThan(ptminsel2, ptminsel1) ||
      !SelectorPtMinGreaterThan(ptminsel2, ptminsel1))
    return 1;
  if (SelectorPtMaxLessThan(ptmaxsel2, ptmaxsel1) ||
      !SelectorPtMaxGreaterThan(ptmaxsel2, ptmaxsel1))
    return 1;

  if (SelectorPtMinLessThan(ptminsel2, ptminsel2) ||
      SelectorPtMinGreaterThan(ptminsel2, ptminsel2))
    return 1;
  if (SelectorPtMaxLessThan(ptmaxsel2, ptmaxsel2) ||
      SelectorPtMaxGreaterThan(ptmaxsel2, ptmaxsel2))
    return 1;

  if (!SelectorPtMinLessEqualThan(ptminsel2, ptminsel2) ||
      !SelectorPtMinGreaterEqualThan(ptminsel2, ptminsel2))
    return 1;
  if (!SelectorPtMaxLessEqualThan(ptmaxsel2, ptmaxsel2) ||
      !SelectorPtMaxGreaterEqualThan(ptmaxsel2, ptmaxsel2))
    return 1;

  
  fastjet::Selector compound_sel = ptminsel1 && ptmaxsel2;
  if (SelectorPtMinGreaterThan(compound_sel, ptminsel2) ||
      !SelectorPtMaxGreaterThan(compound_sel, ptmaxsel1) ||
      !SelectorPtMinLessThan(compound_sel, ptminsel2) ||
      SelectorPtMaxLessThan(compound_sel, ptmaxsel1))
    return 1;
  
  if (!SelectorPtMaxLessEqualThan(compound_sel, ptmaxsel2) ||
      !SelectorPtMinLessEqualThan(compound_sel, ptminsel1) ||
      SelectorPtMaxLessThan(compound_sel, ptmaxsel2) ||
      SelectorPtMinLessThan(compound_sel, ptminsel1))
    return 1;
  
  return 0;
}
