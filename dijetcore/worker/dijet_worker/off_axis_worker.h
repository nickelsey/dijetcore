#ifndef DIJETCORE_WORKER_DIJET_WORKER_OFF_AXIS_WORKER_H
#define DIJETCORE_WORKER_DIJET_WORKER_OFF_AXIS_WORKER_H

#include "dijetcore/worker/dijet_worker/dijet_worker.h"

namespace dijetcore {
  
  struct OffAxisOutput {
    
  };
  
  class OffAxisWorker {
  public:
    OffAxisWorker();
    ~OffAxisWorker();
    
    std::unordered_map<std::string, std::unique_ptr<OffAxisOutput>>& Run(DijetWorker& worker);
    
    std::unordered_map<std::string, std::unique_ptr<OffAxisOutput>>& OffAxisResult() {return off_axis_result;}
    
  private:
    std::unordered_map<std::string, std::unique_ptr<OffAxisOutput>> off_axis_result;
    
  };
  
} // namespace dijetcore

#endif // DIJETCORE_WORKER_DIJET_WORKER_OFF_AXIS_WORKER_H
