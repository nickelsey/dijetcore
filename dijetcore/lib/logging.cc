#include "dijetcore/lib/logging.h"

// initial setting for logging level
DIJETCORE_DEFINE_int(dijetcoreloglevel, google::ERROR,
                    "Minimum severity level to output.");

namespace dijetcore {
  
  bool InitLogging(int* argc, char** argv) {
    if (*argc == 0)
      return true;
    ::google::InitGoogleLogging(argv[0]);

    internal::SetGlogVerbosity();
  
    return true;
  }
  
  bool InitLogging(char* argv) {
    if (argv == nullptr)
      return true;
    ::google::InitGoogleLogging(argv);
    
    internal::SetGlogVerbosity();
    
    return true;
  }
  
  void LogInfoToStdErr() {
    FLAGS_logtostderr = true;
    FLAGS_minloglevel = std::min(FLAGS_minloglevel, google::GLOG_INFO);
  }
  
  namespace internal {
    void SetGlogVerbosity() {
      // choose the lower of dijetcore log level & minloglevel
      // dijetcoreloglevel < 0 for verbose output from glog (FLAGS_v)
      FLAGS_minloglevel = std::min(FLAGS_dijetcoreloglevel, FLAGS_minloglevel);
      if (FLAGS_dijetcoreloglevel <=0) {
        FLAGS_logtostderr = true;
        FLAGS_v = std::max(-FLAGS_dijetcoreloglevel, FLAGS_v);
      }
    }
  }
} // namespace dijetcore
