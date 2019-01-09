#ifndef DIJETCORE_LIB_STRING_STRING_CAST_H
#define DIJETCORE_LIB_STRING_STRING_CAST_H

#include <string>
#include <sstream>

namespace dijetcore {
  
  template<typename T>
  bool CanCast(std::string s) {
    std::istringstream iss(s);
    T dummy;
    iss >> std::skipws >> dummy;
    return iss && iss.eof();
  }
  
  template<typename T>
  T CastTo(std::string s) {
    std::istringstream iss(s);
    T dummy;
    iss >> std::skipws >> dummy;
    return dummy;
  }
  
} // namespace dijetcore

#endif // DIJETCORE_LIB_STRING_STRING_CAST_H
