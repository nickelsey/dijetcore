#include "dijetcore/lib/string/string_utils.h"

namespace dijetcore {
  
  template<>
  void SplitString(const string& str, std::set<string>& cont, char delim) {
    std::stringstream ss(str);
    string token;
    while (std::getline(ss, token, delim)) {
      cont.insert(token);
    }
  }
  
  template<>
  std::string MakeString(const std::string& s) {
    return s;
  }
  
  std::string MakeString(const char* chr) {
    return std::string(chr);
  }
  
  bool HasEnding(const string& full_string, const string& ending) {
    if (full_string.length() >= ending.length()) {
      return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
    } else {
      return false;
    }
  }
  
  bool Consume(std::string& arg, const std::string& seq) {
    if (seq.size() > arg.size()) return false;
    if (arg.compare( 0, seq.length(), seq ) == 0) {
      arg = std::string(arg.begin() + seq.length(), arg.end());
      return true;
    }
    return false;
  }
  
  bool Consume(std::string& arg, char c) {
    if (arg.size() == 0)
      return false;
    if (arg[0] == c) {
      arg = std::string(arg.begin()+1, arg.end());
      return true;
    }
    return false;
  }
  
} // namespace dijetcore

