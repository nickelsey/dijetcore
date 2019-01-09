#ifndef DIJETCORE_LIB_FILESYSTEM_H
#define DIJETCORE_LIB_FILESYSTEM_H

// functions pertaining to filesystem manipulation
// such as creating directories

#include "boost/filesystem.hpp"

#include "dijetcore/lib/types.h"

namespace dijetcore {

  bool CreateDirectory(const string& path);

} // namespace dijetcore

#endif // DIJETCORE_LIB_FILESYSTEM_H