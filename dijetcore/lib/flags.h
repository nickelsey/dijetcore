
#ifndef DIJETCORE_CORE_FLAGS_HH
#define DIJETCORE_CORE_FLAGS_HH

/* implementation of command line argument parser
 * using gflags
 */

#include "dijetcore/lib/types.h"

namespace dijetcore {
  
  // returns the usage message
  const char* ProgramUsage();
  
  // sets the usage message that is printed when the executable
  // is called with flag --help. Do not add lines for individual
  // flags, this is handled by gflags
  void SetUsageMessage(const string& str);
  
  // Parses command line options. Follows gflags convention,
  // where argv and argc are modified to remove any flags that
  // have been handled, leaving any unhandled flags
  bool ParseCommandLineFlags(int* argc, char** argv);
  
} // namespace dijetcore

#include "gflags/gflags.h"

// including gflags declared variables in the dijetcore namespace

#define DIJETCORE_GFLAGS_DEFINE(type, name, default_value, help_str)             \
  DEFINE_##type(name, default_value, help_str);                                 \
  namespace dijetcore {                                                          \
    using ::FLAGS_##name;                                                       \
  } // namespace dijetcore

#define DIJETCORE_DEFINE_int(name, default_value, help_str)                      \
  DIJETCORE_GFLAGS_DEFINE(int32, name, default_value, help_str)
#define DIJETCORE_DEFINE_int64(name, default_value, help_str)                    \
  DIJETCORE_GFLAGS_DEFINE(int64, name, default_value, help_str)
#define DIJETCORE_DEFINE_double(name, default_value, help_str)                   \
  DIJETCORE_GFLAGS_DEFINE(double, name, default_value, help_str)
#define DIJETCORE_DEFINE_bool(name, default_value, help_str)                     \
  DIJETCORE_GFLAGS_DEFINE(bool, name, default_value, help_str)
#define DIJETCORE_DEFINE_string(name, default_value, help_str)                   \
  DIJETCORE_GFLAGS_DEFINE(string, name, default_value, help_str)

#define DIJETCORE_GFLAGS_DECLARE(type, name)            \
  DECLARE_##type(name);                                \
  namespace dijetcore {                                 \
    using ::FLAGS_##name;                              \
  } // namespace dijetcore


#define DIJETCORE_DECLARE_int(name)                    \
  DIJETCORE_GFLAGS_DECLARE(int32, name)
#define DIJETCORE_DECLARE_int64(name)                  \
  DIJETCORE_GFLAGS_DECLARE(int64, name)
#define DIJETCORE_DECLARE_double(name)                 \
  DIJETCORE_GFLAGS_DECLARE(double, name)
#define DIJETCORE_DECLARE_bool(name)                   \
  DIJETCORE_GFLAGS_DECLARE(bool, name)
#define DIJETCORE_DECLARE_string(name)                 \
  DIJETCORE_GFLAGS_DECLARE(string, name)
  

#endif // DIJETCORE_CORE_FLAGS_HH
