#include "dijetcore/lib/assert.h"

#include <algorithm>
#include <numeric>
#include <cstring>

#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/string/path_utils.h"

namespace dijetcore {
  
  AssertionFailure::AssertionFailure(const char* file,
                                     int line,
                                     const char* failure,
                                     const string& msg) :
    msg_stack_({MakeString("[assertion failure: ", GetFileName(string(file)), "::", line, "] "),
               string(failure), " ", msg}), msg_(message()) {}
  
  string AssertionFailure::message() {
    return std::accumulate(msg_stack_.begin(), msg_stack_.end(), string(" "));
  }
  
  void AssertionFailure::append(const string& str) {
    msg_stack_.push_back(str);
    msg_ = message();
  }
  
  const char* AssertionFailure::what() const noexcept {
    return msg_.c_str();
  }
  
} // namespace dijetcore
