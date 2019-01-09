#ifndef DIJETCORE_LIB_FILESYSTEM_H
#define DIJETCORE_LIB_FILESYSTEM_H

// functions pertaining to filesystem manipulation
// such as creating directories

#include "dijetcore/lib/types.h"

namespace dijetcore {

  bool CreateDirectory(const string& path);

  // can concatenate two paths or path & file
  string ConcatenatePath(const string& path1, const string& path2);

} // namespace dijetcore

#endif // DIJETCORE_LIB_FILESYSTEM_H