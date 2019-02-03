#include "dijetcore/lib/filesystem.h"

#include "boost/filesystem.hpp"

namespace dijetcore {

  bool CreateDirectory(const string& path) {
    boost::filesystem::path dir(path.c_str());
    boost::filesystem::create_directories(dir);

    return boost::filesystem::is_directory(dir);
  }

  string ConcatenatePath(const string& path1, const string& path2) {
    boost::filesystem::path concat_path(path1.c_str());
    concat_path /= path2;
    return concat_path.string();
  }

} // namespace dijetcore