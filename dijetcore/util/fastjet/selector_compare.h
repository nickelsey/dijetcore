#ifndef DIJETCORE_UTIL_FASTJET_SELECTOR_COMPARE_H
#define DIJETCORE_UTIL_FASTJET_SELECTOR_COMPARE_H

// used to compare selectors based on their string description
// to keep things as simple as possible, I have decided that
// when applicable, these will act like operators
// (eg, SelectorPtMinGreaterThan(lhs, rhs) will return true if lhs' pt cut
// is greater than rhs's).
// To avoid throwing exceptions, it will be assumed that if a selector
// does not contain a specific cut, then its value is taken to be zero

#include "fastjet/Selector.hh"

#include "dijetcore/lib/types.h"

namespace dijetcore {
  
  double ExtractDoubleFromSelector(const fastjet::Selector& sel, string descriptor);
  
  bool SelectorPtMinGreaterThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs);
  
  bool SelectorPtMinGreaterEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs);
  
  bool SelectorPtMinLessThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs);
  
  bool SelectorPtMinLessEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs);
  
  bool SelectorPtMaxGreaterThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs);
  
  bool SelectorPtMaxGreaterEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs);
  
  bool SelectorPtMaxLessThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs);
  
  bool SelectorPtMaxLessEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs);
  
} // namespace dijetcore

#endif //DIJETCORE_UTIL_FASTJET_SELECTOR_COMPARE_H
