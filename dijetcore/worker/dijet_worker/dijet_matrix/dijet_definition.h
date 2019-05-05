#ifndef DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_DIJET_DEFINITION_H
#define DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_DIJET_DEFINITION_H

// A DijetDefinition contains a MatchDef for
// both the leading & subleading jet

// direct access is given to the MatchDefs, since
// this is just a container

#include "dijetcore/lib/memory.h"

#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"

namespace dijetcore {
  
  class DijetDefinition {
  public:
    DijetDefinition() : lead(nullptr), sub(nullptr), dPhi(0.4) { };
    DijetDefinition(const DijetDefinition& rhs) :
    lead(rhs.lead), sub(rhs.sub), dPhi(rhs.dPhi) { };
    DijetDefinition(MatchDef* lead, MatchDef* sub, double phi = 0.4) :
    lead(lead), sub(sub), dPhi(phi) { };
    DijetDefinition(unique_ptr<MatchDef> lead, unique_ptr<MatchDef>sub, double phi = 0.4) :
    lead(lead.get()), sub(sub.get()), dPhi(phi) { };
    
    virtual ~DijetDefinition() { };
    
    inline bool isValid()    const {return lead && sub && lead->isValid() && sub->isValid();}
    inline bool doMatching() const {return lead && sub && lead->canMatch() && sub->canMatch();}
    
    // checks if two dijet definitions can use equivalent
    // clustersequences
    bool equivalentCluster(const DijetDefinition& rhs) const;
    
    MatchDef* lead;
    MatchDef* sub;
    
    double dPhi;
  };
  
  // tests for logical inequality, not pointer
  // inequality - so if the pointers differ, but
  // the underlying MatchDefs are equivalent, will
  // return true
  inline bool operator==(const DijetDefinition& lhs, const DijetDefinition& rhs) {
    if (lhs.dPhi != rhs.dPhi)
      return false;
    if (lhs.lead == nullptr && lhs.sub == nullptr &&
        rhs.lead == nullptr && rhs.sub == nullptr)
      return true;
    if ((lhs.lead == nullptr && rhs.lead != nullptr) ||
        (lhs.lead != nullptr && rhs.lead == nullptr) ||
        (lhs.sub == nullptr  && rhs.sub != nullptr) ||
        (lhs.sub != nullptr  && rhs.sub == nullptr))
      return false;
    if (lhs.lead != nullptr &&
        *lhs.lead != *rhs.lead)
      return false;
    if (lhs.sub != nullptr &&
        *lhs.sub != *rhs.sub)
      return false;
    return true;
  }
  inline bool operator!=(const DijetDefinition& lhs, const DijetDefinition& rhs) {
    return !(lhs == rhs);
  }
  
} // namespace dijetcore

#endif // DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_DIJET_DEFINITION_H
