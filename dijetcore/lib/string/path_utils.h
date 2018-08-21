#ifndef DIJETCORE_LIB_STRING_PATH_UTILS_H
#define DIJETCORE_LIB_STRING_PATH_UTILS_H

// functionality to manipulate strings representing paths
// and file names

#include "dijetcore/lib/types.h"

namespace dijetcore {

    // strips the path and returns the file name
    string GetFileName(const string& path);
    
    // strips the file name and returns the path
    string GetPath(const string& path);

} // namespace dijetcore


#endif // DIJETCORE_LIB_STRING_PATH_UTILS_H
