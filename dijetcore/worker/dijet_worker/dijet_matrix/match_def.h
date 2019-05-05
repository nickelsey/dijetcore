#ifndef DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_MATCH_DEF_H
#define DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_MATCH_DEF_H

// a set of jet definitions intended to be used to match jets
// between cluster sequences - for instance, in the STAR A_J,
// where a hard pT const cut is applied first, then matched to
// jets clustered with the full event

#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"

namespace dijetcore {
  
  class MatchDef {
  public:
    
    // By default, both initial & match definitions
    // are in invalid states to prevent misunderstanding
    MatchDef() : initial_(), matched_(), R_(0.0) { };
    
    // if not specified, dR is set to the minimum of the
    // two JetDef radii
    MatchDef(const JetDef& init, const JetDef& match) :
    initial_(init), matched_(match) {
      R_ = initial_.R() < matched_.R() ? initial_.R() : matched_.R();
    }
    
    // can specify a non-default radius
    MatchDef(const JetDef& init, const JetDef& match,
             double dR) : initial_(init), matched_(match), R_(dR) { };
    
    virtual ~MatchDef() { }
    
    inline JetDef initialJetDef() const {return initial_;}
    inline JetDef matchedJetDef() const {return matched_;}
    
    inline double  dR() const {return R_;}
    
    void setInitialJetDef(const JetDef& def) {initial_ = def;}
    void setMatchedJetDef(const JetDef& def) {matched_ = def;}
    
    void setdR(double dR_in) {R_ = dR_in;}
    
    // matching is done when a both JetDefs are in valid states
    inline bool isValid()  const  {return initial_.isValid();}
    inline bool canMatch() const  {return isValid() && matched_.isValid();}
    
    // decide if two matchdefs can use the same cluster sequences
    // for both initial and matched jetfinding
    bool equivalentCluster(const MatchDef& rhs) const;
    
  private:
    
    JetDef initial_;
    JetDef matched_;
    
    double R_;
    
  };
  
  inline bool operator==(const MatchDef& lhs, const MatchDef& rhs) {
    if (lhs.initialJetDef() != rhs.initialJetDef() ||
        lhs.matchedJetDef() != rhs.matchedJetDef() ||
        lhs.dR() != rhs.dR())
      return false;
    return true;
  }
  
  inline bool operator!=(const MatchDef& lhs, const MatchDef& rhs) {
    return !(lhs == rhs);
  }

} // namespace dijetcore
  
#endif // DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_MATCH_DEF_H
