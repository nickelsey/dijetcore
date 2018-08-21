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
    inline fastjet::Selector ConstituentSelector()   const {return const_selector_;}
    
    // selector used on jets post-clustering
    // (and post bkg subtraction, if it is being applied)
    inline fastjet::Selector JetSelector()           const {return jet_selector_;}
    
    // area definition used for area estimation &
    // background subtraction
    inline fastjet::AreaDefinition AreaDefinition()  const {return area_def_;}
    
    // background jet definition
    inline fastjet::JetDefinition BackgroundJetDef() const {return bkg_jet_def_;}
    
    // background area definition
    inline fastjet::AreaDefinition BackgroundAreaDef() const {return bkg_area_def_;}
    
    // selector used to estimate background. By default,
    // set to !SelectorNHardest(2), but should also include
    // a user defined rapidity cut
    inline fastjet::Selector BackgroundSelector()    const {return bkg_selector_;}
    
    // set the various fastjet definitions & selectors
    inline void SetConstituentSelector(fastjet::Selector sel)     {const_selector_ = sel;}
    inline void SetJetSelector(fastjet::Selector sel)             {jet_selector_ = sel;}
    inline void SetAreaDefinition(fastjet::AreaDefinition def)    {area_def_ = def;}
    inline void SetBackgroundJetDef(fastjet::JetDefinition def)   {bkg_jet_def_ = def;}
    inline void SetBackgroundAreaDef(fastjet::AreaDefinition def) {bkg_area_def_ = def;}
    inline void SetBackgroundSelector(fastjet::Selector sel)      {bkg_selector_ = sel;}
    
    // returns true if the jetdefinition is valid
    inline bool IsValid() const          {return jet_algorithm() != fastjet::undefined_jet_algorithm;}
    
    // returns true if a valid area definition is present
    inline bool CanEstimateArea() const  {return area_def_.area_type() != fastjet::invalid_area;}
    
    // returns true if a valid area definition & valid
    // background jet definition are defined
    inline bool CanBackgroundSub() {return CanEstimateArea() &&
      (bkg_jet_def_.jet_algorithm() != fastjet::undefined_jet_algorithm) &&
      (bkg_area_def_.area_type() != fastjet::invalid_area);}
    
    // returns true if two JetDefs use an equivalent
    // clustersequence & input particles, which means that
    // clustering can be done only once
    bool EquivalentCluster(const JetDef& rhs) const;
    
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
        lhs.ConstituentSelector().description() != rhs.ConstituentSelector().description() ||
        lhs.JetSelector().description() != rhs.JetSelector().description() ||
        lhs.AreaDefinition().area_type() != rhs.AreaDefinition().area_type() ||
        lhs.AreaDefinition().ghost_spec().ghost_maxrap() != rhs.AreaDefinition().ghost_spec().ghost_maxrap() ||
        lhs.AreaDefinition().ghost_spec().ghost_area() != rhs.AreaDefinition().ghost_spec().ghost_area() ||
        lhs.AreaDefinition().ghost_spec().grid_scatter() != rhs.AreaDefinition().ghost_spec().grid_scatter() ||
        lhs.AreaDefinition().ghost_spec().pt_scatter() != rhs.AreaDefinition().ghost_spec().pt_scatter() ||
        lhs.AreaDefinition().ghost_spec().mean_ghost_pt() != rhs.AreaDefinition().ghost_spec().mean_ghost_pt() ||
        lhs.AreaDefinition().ghost_spec().repeat() != rhs.AreaDefinition().ghost_spec().repeat() ||
        lhs.BackgroundAreaDef().area_type() != rhs.BackgroundAreaDef().area_type() ||
        lhs.BackgroundAreaDef().ghost_spec().ghost_maxrap() != rhs.BackgroundAreaDef().ghost_spec().ghost_maxrap() ||
        lhs.BackgroundAreaDef().ghost_spec().ghost_area() != rhs.BackgroundAreaDef().ghost_spec().ghost_area() ||
        lhs.BackgroundAreaDef().ghost_spec().grid_scatter() != rhs.BackgroundAreaDef().ghost_spec().grid_scatter() ||
        lhs.BackgroundAreaDef().ghost_spec().pt_scatter() != rhs.BackgroundAreaDef().ghost_spec().pt_scatter() ||
        lhs.BackgroundAreaDef().ghost_spec().mean_ghost_pt() != rhs.BackgroundAreaDef().ghost_spec().mean_ghost_pt() ||
        lhs.BackgroundAreaDef().ghost_spec().repeat() != rhs.BackgroundAreaDef().ghost_spec().repeat() ||
        lhs.BackgroundJetDef().jet_algorithm() != rhs.BackgroundJetDef().jet_algorithm() ||
        lhs.BackgroundJetDef().R() != rhs.BackgroundJetDef().R() ||
        lhs.BackgroundJetDef().recombination_scheme() != rhs.BackgroundJetDef().recombination_scheme() ||
        lhs.BackgroundJetDef().strategy() != rhs.BackgroundJetDef().strategy() ||
        lhs.BackgroundSelector().description() != rhs.BackgroundSelector().description())
      return false;
    return true;
  }
  
  inline bool operator!=(const JetDef& lhs, const JetDef& rhs) {
    return !(lhs == rhs);
  }
  
} // namespace dijetcore

#endif  // DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_JET_DEF_H
