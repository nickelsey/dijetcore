#include "dijetcore/lib/filesystem.h"

#include "boost/filesystem.hpp"

namespace dijetcore {

  bool CreateDirectory(const string& path) {
    boost::filesystem::path dir(path.c_str());
    return boost::filesystem::create_directories(dir);
  }

} // namespace dijetcore