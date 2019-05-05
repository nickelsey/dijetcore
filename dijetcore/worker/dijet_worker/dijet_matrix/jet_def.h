#ifndef DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_JET_DEF_H
#define DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_JET_DEF_H

// extension of the fastjet::JetDefinition class
// to encapsulate background subtraction options

#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/Selector.hh"

namespace dijetcore {
  
  class JetDef : public fastjet::JetDefinition {
  public:
    // default constructor has invalid settings
    // to prevent jetfinding without explicitly
    // choosing an algorithm and radius
    JetDef();
    
    // initializes the fastjet::JetDefinition,
    // by default,invalid area estimation &
    // invalid bkg subtraction settings
    JetDef(fastjet::JetAlgorithm alg, double R,
           fastjet::AreaDefinition area_def =
           fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec()),
           fastjet::JetDefinition bkg_jet_def = fastjet::JetDefinition(),
           fastjet::AreaDefinition bkg_area_def =
           fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec()),
           fastjet::Selector bkg_selector = !fastjet::SelectorNHardest(2));
    
    JetDef(fastjet::JetAlgorithm alg, double R,
           double xtra_param,
           fastjet::AreaDefinition area_def =
           fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec()),
           fastjet::JetDefinition bkg_jet_def = fastjet::JetDefinition(),
           fastjet::AreaDefinition bkg_area_def_ =
           fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec()),
           fastjet::Selector bkg_selector = !fastjet::SelectorNHardest(2));
    
    JetDef(fastjet::JetDefinition jet_def_in,
           fastjet::AreaDefinition area_def =
           fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec()),
           fastjet::JetDefinition bkg_jet_def = fastjet::JetDefinition(),
           fastjet::AreaDefinition bkg_area_def =
           fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec()),
           fastjet::Selector bkg_selector = !fastjet::SelectorNHardest(2));
    
    JetDef(const JetDef& rhs);
    
    virtual ~JetDef() { }
    
    // selector used on constituents pre-clustering
    inline fastjet::Selector constituentSelector()   const {return const_selector_;}
    
    // selector used on jets post-clustering
    // (and post bkg subtraction, if it is being applied)
    inline fastjet::Selector jetSelector()           const {return jet_selector_;}
    
    // area definition used for area estimation &
    // background subtraction
    inline fastjet::AreaDefinition areaDefinition()  const {return area_def_;}
    
    // background jet definition
    inline fastjet::JetDefinition backgroundJetDef() const {return bkg_jet_def_;}
    
    // background area definition
    inline fastjet::AreaDefinition backgroundAreaDef() const {return bkg_area_def_;}
    
    // selector used to estimate background. By default,
    // set to !SelectorNHardest(2), but should also include
    // a user defined rapidity cut
    inline fastjet::Selector backgroundSelector()    const {return bkg_selector_;}
    
    // set the various fastjet definitions & selectors
    inline void setConstituentSelector(fastjet::Selector sel)     {const_selector_ = sel;}
    inline void setJetSelector(fastjet::Selector sel)             {jet_selector_ = sel;}
    inline void setAreaDefinition(fastjet::AreaDefinition def)    {area_def_ = def;}
    inline void setBackgroundJetDef(fastjet::JetDefinition def)   {bkg_jet_def_ = def;}
    inline void setBackgroundAreaDef(fastjet::AreaDefinition def) {bkg_area_def_ = def;}
    inline void setBackgroundSelector(fastjet::Selector sel)      {bkg_selector_ = sel;}
    
    // returns true if the jetdefinition is valid
    inline bool isValid() const          {return jet_algorithm() != fastjet::undefined_jet_algorithm;}
    
    // returns true if a valid area definition is present
    inline bool canEstimateArea() const  {return area_def_.area_type() != fastjet::invalid_area;}
    
    // returns true if a valid area definition & valid
    // background jet definition are defined
    inline bool canBackgroundSub() {return canEstimateArea() &&
      (bkg_jet_def_.jet_algorithm() != fastjet::undefined_jet_algorithm) &&
      (bkg_area_def_.area_type() != fastjet::invalid_area);}
    
    // returns true if two JetDefs use an equivalent
    // clustersequence & input particles, which means that
    // clustering can be done only once
    bool equivalentCluster(const JetDef& rhs) const;
    
  private:
    
    // Selector for constituents pre-clustering
    // by default is an identity function
    fastjet::Selector const_selector_;
    
    // Selector for jets post-clustering
    // by default is an identity function
    fastjet::Selector jet_selector_;
    
    // area estimation tools
    
    // the area definition used for area estimation.
    // don't use FastJet's default GhostedAreaSpec for
    // STAR, since the default ghosted area is much larger
    // in eta than necessary, and will increase run time.
    fastjet::AreaDefinition area_def_;
    
    // background subtraction tools
    
    // jet definition used for background subtraction.
    // a good starting point (per FastJet) is to use
    // the kt algorithm, with R = 0.4-0.6
    fastjet::JetDefinition bkg_jet_def_;
    
    // an area definition that can be used with the background
    // jet definition - really only needs to differ from the area
    // definition of the other area definition if the two
    // JetDefinitions have very different R
    fastjet::AreaDefinition bkg_area_def_;
    
    // selector used by the JetMedianBackgroundEstimator,
    // by default, set to !SelectorNHardest(2)
    fastjet::Selector bkg_selector_;
    
  };
  
  // we use the textual description of selectors to compare, since fastjet
  // does not provide an equality operator. Therefore, compound Selectors
  // may give equivalent outputs but still not be considered equal.
  inline bool operator==(const JetDef& lhs, const JetDef& rhs) {
    if (lhs.jet_algorithm() != rhs.jet_algorithm() ||
        lhs.R() != rhs.R() ||
        lhs.recombination_scheme() != rhs.recombination_scheme() ||
        lhs.strategy() != rhs.strategy() ||
        lhs.constituentSelector().description() != rhs.constituentSelector().description() ||
        lhs.jetSelector().description() != rhs.jetSelector().description() ||
        lhs.areaDefinition().area_type() != rhs.areaDefinition().area_type() ||
        lhs.areaDefinition().ghost_spec().ghost_maxrap() != rhs.areaDefinition().ghost_spec().ghost_maxrap() ||
        lhs.areaDefinition().ghost_spec().ghost_area() != rhs.areaDefinition().ghost_spec().ghost_area() ||
        lhs.areaDefinition().ghost_spec().grid_scatter() != rhs.areaDefinition().ghost_spec().grid_scatter() ||
        lhs.areaDefinition().ghost_spec().pt_scatter() != rhs.areaDefinition().ghost_spec().pt_scatter() ||
        lhs.areaDefinition().ghost_spec().mean_ghost_pt() != rhs.areaDefinition().ghost_spec().mean_ghost_pt() ||
        lhs.areaDefinition().ghost_spec().repeat() != rhs.areaDefinition().ghost_spec().repeat() ||
        lhs.backgroundAreaDef().area_type() != rhs.backgroundAreaDef().area_type() ||
        lhs.backgroundAreaDef().ghost_spec().ghost_maxrap() != rhs.backgroundAreaDef().ghost_spec().ghost_maxrap() ||
        lhs.backgroundAreaDef().ghost_spec().ghost_area() != rhs.backgroundAreaDef().ghost_spec().ghost_area() ||
        lhs.backgroundAreaDef().ghost_spec().grid_scatter() != rhs.backgroundAreaDef().ghost_spec().grid_scatter() ||
        lhs.backgroundAreaDef().ghost_spec().pt_scatter() != rhs.backgroundAreaDef().ghost_spec().pt_scatter() ||
        lhs.backgroundAreaDef().ghost_spec().mean_ghost_pt() != rhs.backgroundAreaDef().ghost_spec().mean_ghost_pt() ||
        lhs.backgroundAreaDef().ghost_spec().repeat() != rhs.backgroundAreaDef().ghost_spec().repeat() ||
        lhs.backgroundJetDef().jet_algorithm() != rhs.backgroundJetDef().jet_algorithm() ||
        lhs.backgroundJetDef().R() != rhs.backgroundJetDef().R() ||
        lhs.backgroundJetDef().recombination_scheme() != rhs.backgroundJetDef().recombination_scheme() ||
        lhs.backgroundJetDef().strategy() != rhs.backgroundJetDef().strategy() ||
        lhs.backgroundSelector().description() != rhs.backgroundSelector().description())
      return false;
    return true;
  }
  
  inline bool operator!=(const JetDef& lhs, const JetDef& rhs) {
    return !(lhs == rhs);
  }
  
} // namespace dijetcore

#endif  // DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_JET_DEF_H
