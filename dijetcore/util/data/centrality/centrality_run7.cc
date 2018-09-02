#include "dijetcore/util/data/centrality/centrality_run7.h"

namespace dijetcore {
  
  CentralityRun7::CentralityRun7() {
    cent_bin_9_ = std::vector<unsigned>{485, 399, 269, 178, 114, 69, 39, 21, 10};
  }
  
  CentralityRun7::~CentralityRun7() {}
  
  int CentralityRun7::Centrality9(unsigned grefmult) {
    for (unsigned i = 0; i < cent_bin_9_.size(); ++i) {
      if (grefmult >= cent_bin_9_[i])
      return i;
    }
    return -1;
  }
  
} // namespace dijetcore
