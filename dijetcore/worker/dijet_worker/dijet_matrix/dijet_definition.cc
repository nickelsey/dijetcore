#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"

namespace dijetcore {
  
  bool DijetDefinition::equivalentCluster(const DijetDefinition& rhs) const {
    return lead->equivalentCluster(*rhs.lead) && sub->equivalentCluster(*rhs.sub);
  }
  
} // namespace dijetcore
