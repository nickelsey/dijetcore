#ifndef DIJETCORE_LIB_MEMORY_H
#define DIJETCORE_LIB_MEMORY_H

// memory includes for libdijetcore's internal use

#include <memory>

// In case we're compiling with c++11, we will emulate c++14's std::make_unique                                                                                                                                                                                                 
#if __cplusplus < 201402L && (!defined __cpp_lib_make_unique)
  #include "dijetcore/lib/memory/make_unique.h"
#endif

namespace dijetcore {
    // commonly used memory functionality
    using std::unique_ptr;
    using std::shared_ptr;
    using std::make_unique;
    using std::make_shared;
}

#endif // DIJETCORE_LIB_MEMORY_H
