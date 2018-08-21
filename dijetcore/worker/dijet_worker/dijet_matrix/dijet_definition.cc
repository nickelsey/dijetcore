#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"

namespace dijetcore {
  
  bool DijetDefinition::EquivalentCluster(const DijetDefinition& rhs) const {
    return lead->EquivalentCluster(*rhs.lead) && sub->EquivalentCluster(*rhs.sub);
  }
  
} // namespace dijetcore
