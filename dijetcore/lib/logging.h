
#ifndef DIJETCORE_CORE_LOGGING_HH
#define DIJETCORE_CORE_LOGGING_HH

/* dijetcore uses GLog for logging, style is 
 * LOG(INFO) <<... LOG(ERROR) << ... etc
 * for more info on using glog, look in 
 * https://github.com/google/glog
 */

#include "glog/stl_logging.h"
#include "glog/logging.h"

#include "dijetcore/lib/flags.h"

// used to control levels of output using glog
// levels: INFO = 0, WARNING = 1, ERROR = 2,
// FATAL = 3. Fatal will cause exit if NDEBUG is turned on
// for dijetcoreloglevel, negative values represent glog's
// verbosity parameter, for extra output, mainly during debugging.
DIJETCORE_DECLARE_int(dijetcoreloglevel);
DIJETCORE_DECLARE_int(minloglevel);
// logs all output to stderr
DIJETCORE_DECLARE_bool(logtostderr);
  
namespace dijetcore {
  
  bool InitLogging(int* argc, char** argv);
  bool InitLogging(char* argv);
  
  void LogInfoToStdErr();
  
  namespace internal {
    void SetGlogVerbosity();
  }
  
} // namespace dijetcore

#endif // DIJETCORE_CORE_LOGGING_HH
