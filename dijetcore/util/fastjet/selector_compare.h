#ifndef DIJETCORE_UTIL_FASTJET_SELECTOR_COMPARE_H
#define DIJETCORE_UTIL_FASTJET_SELECTOR_COMPARE_H

// used to compare selectors based on their string description
// to keep things as simple as possible, I have decided that
// when applicable, these will act like operators
// (eg, SelectorPtMinGreaterThan(lhs, rhs) will return true if lhs' pt cut
// is greater than rhs's).
// To avoid throwing exceptions, it will be assumed that if a selector
// does not contain a specific cut, then its value is taken to be zero

#include <sstream>

#include "fastjet/Selector.hh"

#include "dijetcore/lib/types.h"

namespace dijetcore {
  
  double ExtractDoubleFromSelector(const fastjet::Selector& sel, string descriptor) {
    
    string sel_str = sel.description();
    std::size_t length = descriptor.size();
    std::size_t pos    = sel_str.find(descriptor);
    
    if ( pos == string::npos)
      return 0.0;
    
    std::istringstream stream(string(sel_str.begin() + pos + length, sel_str.end()));
    double ret;
    stream >> std::skipws >> ret;
    return ret;
  }
  
  bool SelectorPtMinGreaterThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
    string descriptor = "pt >=";
    
    double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
    double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
    
    return lhs_pt > rhs_pt;
  }
  
  bool SelectorPtMinGreaterEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
    string descriptor = "pt >=";
    
    double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
    double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
    
    return lhs_pt >= rhs_pt;
  }
  
  bool SelectorPtMinLessThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
    string descriptor = "pt >=";
    
    double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
    double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
    
    return lhs_pt < rhs_pt;
  }
  
  bool SelectorPtMinLessEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
    string descriptor = "pt >=";
    
    double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
    double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
    
    return lhs_pt <= rhs_pt;
  }
  
  bool SelectorPtMaxGreaterThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
    string descriptor = "pt <=";
    
    double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
    double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
    
    return lhs_pt > rhs_pt;
  }
  
  bool SelectorPtMaxGreaterEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
    string descriptor = "pt <=";
    
    double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
    double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
    
    return lhs_pt >= rhs_pt;
  }
  
  bool SelectorPtMaxLessThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
    string descriptor = "pt <=";
    
    double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
    double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
    
    return lhs_pt < rhs_pt;
  }
  
  bool SelectorPtMaxLessEqualThan(const fastjet::Selector& lhs, const fastjet::Selector& rhs) {
    string descriptor = "pt <=";
    
    double lhs_pt = ExtractDoubleFromSelector(lhs, descriptor);
    double rhs_pt = ExtractDoubleFromSelector(rhs, descriptor);
    
    return lhs_pt <= rhs_pt;
  }
  
} // namespace dijetcore

#endif //DIJETCORE_UTIL_FASTJET_SELECTOR_COMPARE_H
