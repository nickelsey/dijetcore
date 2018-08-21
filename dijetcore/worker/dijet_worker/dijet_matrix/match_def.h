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
    
    inline JetDef InitialJetDef() const {return initial_;}
    inline JetDef MatchedJetDef() const {return matched_;}
    
    inline double  dR() const {return R_;}
    
    void SetInitialJetDef(const JetDef& def) {initial_ = def;}
    void SetMatchedJetDef(const JetDef& def) {matched_ = def;}
    
    void SetdR(double dR_in) {R_ = dR_in;}
    
    // matching is done when a both JetDefs are in valid states
    inline bool IsValid()  const  {return initial_.IsValid();}
    inline bool CanMatch() const  {return IsValid() && matched_.IsValid();}
    
    // decide if two matchdefs can use the same cluster sequences
    // for both initial and matched jetfinding
    bool EquivalentCluster(const MatchDef& rhs) const;
    
  private:
    
    JetDef initial_;
    JetDef matched_;
    
    double R_;
    
  };
  
  inline bool operator==(const MatchDef& lhs, const MatchDef& rhs) {
    if (lhs.InitialJetDef() != rhs.InitialJetDef() ||
        lhs.MatchedJetDef() != rhs.MatchedJetDef() ||
        lhs.dR() != rhs.dR())
      return false;
    return true;
  }
  
  inline bool operator!=(const MatchDef& lhs, const MatchDef& rhs) {
    return !(lhs == rhs);
  }

} // namespace dijetcore
  
#endif // DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_MATCH_DEF_H
