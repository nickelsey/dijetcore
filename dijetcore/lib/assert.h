#ifndef DIJETCORE_LIB_ASSERT_H
#define DIJETCORE_LIB_ASSERT_H

// verbose exception handling for dijetcore

#include <exception>
#include <vector>

#include "dijetcore/lib/types.h"

namespace dijetcore {
    // assertion can be anything that will be evaluated to a boolean
    // failure will throw an AssertionFailure, which provides some verbose
    // error logging (where it was thrown from in the source code, for instance)
 #define DIJETCORE_ASSERT(assertion, ...)                                    \
if (!(assertion)) {                                                        \
  throw ::dijetcore::AssertionFailure(                                      \
    __FILE__, __LINE__, #assertion , ::dijetcore::MakeString(__VA_ARGS__)); \
}

    // throws an AssertionFailure, as described above
#define DIJETCORE_THROW(...)                                     \
{ throw ::dijetcore::AssertionFailure(                           \
  __FILE__, __LINE__, "", ::dijetcore::MakeString(__VA_ARGS__)); \
}
    
    // can collect verbose messages recursively through a call stack 
    // of try/catch blocks, using Append
  class AssertionFailure : public std::exception {
  public:
    AssertionFailure(const char* file,
                     int line,
                     const char* failure,
                     const string& msg);
    
    const std::vector<string>& MessageStack() const {return msg_stack_;}
    
    void Append(const string& str);
    
    string Msg();
    
    const char* what() const noexcept override;
    
  private:
    std::vector<string> msg_stack_;
    string msg_;
  };
    
} // namespace dijetcore

#endif // DIJETCORE_LIB_ASSERT_H
