#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"

namespace dijetcore {
  
  bool MatchDef::equivalentCluster(const MatchDef& rhs) const {
    return rhs.initialJetDef().equivalentCluster(initialJetDef()) &&
    rhs.matchedJetDef().equivalentCluster(matchedJetDef());
  }
  
} // namespace dijetcore
