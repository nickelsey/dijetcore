#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"

namespace dijetcore {
  
  bool MatchDef::EquivalentCluster(const MatchDef& rhs) const {
    return rhs.InitialJetDef().EquivalentCluster(InitialJetDef()) &&
    rhs.MatchedJetDef().EquivalentCluster(InitialJetDef());
  }
  
} // namespace dijetcore
