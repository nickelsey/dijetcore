#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"

namespace dijetcore {

  JetDef::JetDef() :
  fastjet::JetDefinition(),
  const_selector_(fastjet::SelectorIdentity()),
  jet_selector_(fastjet::SelectorIdentity()),
  area_def_(fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec())),
  bkg_jet_def_(fastjet::JetDefinition()),
  bkg_area_def_(fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec())),
  bkg_selector_(!fastjet::SelectorNHardest(2)) { }
  
  JetDef::JetDef(fastjet::JetAlgorithm alg, double R,
                 fastjet::AreaDefinition area_def,
                 fastjet::JetDefinition bkg_jet_def,
                 fastjet::AreaDefinition bkg_area_def,
                 fastjet::Selector bkg_selector) :
  fastjet::JetDefinition(alg, R),
  const_selector_(fastjet::SelectorIdentity()),
  jet_selector_(fastjet::SelectorIdentity()),
  area_def_(area_def),
  bkg_jet_def_(bkg_jet_def),
  bkg_area_def_(bkg_area_def),
  bkg_selector_(bkg_selector) { }
  
  JetDef::JetDef(fastjet::JetAlgorithm alg, double R,
                 double xtra_param,
                 fastjet::AreaDefinition area_def,
                 fastjet::JetDefinition bkg_jet_def,
                 fastjet::AreaDefinition bkg_area_def,
                 fastjet::Selector bkg_selector) :
  fastjet::JetDefinition(alg, R, xtra_param),
  const_selector_(fastjet::SelectorIdentity()),
  jet_selector_(fastjet::SelectorIdentity()),
  area_def_(area_def),
  bkg_jet_def_(bkg_jet_def),
  bkg_area_def_(bkg_area_def),
  bkg_selector_(bkg_selector) { }
  
  JetDef::JetDef(fastjet::JetDefinition jet_def_in,
                 fastjet::AreaDefinition area_def,
                 fastjet::JetDefinition bkg_jet_def,
                 fastjet::AreaDefinition bkg_area_def,
                 fastjet::Selector bkg_selector) :
  fastjet::JetDefinition(jet_def_in),
  const_selector_(fastjet::SelectorIdentity()),
  jet_selector_(fastjet::SelectorIdentity()),
  area_def_(area_def),
  bkg_jet_def_(bkg_jet_def),
  bkg_area_def_(bkg_area_def),
  bkg_selector_(bkg_selector) { }
  
  JetDef::JetDef(const JetDef& rhs) :
  fastjet::JetDefinition(rhs.jet_algorithm(), rhs.R(),
                         rhs.recombination_scheme(),
                         rhs.strategy()),
  const_selector_(rhs.const_selector_),
  jet_selector_(rhs.jet_selector_),
  area_def_(rhs.area_def_),
  bkg_jet_def_(rhs.bkg_jet_def_),
  bkg_area_def_(rhs.bkg_area_def_),
  bkg_selector_(rhs.bkg_selector_) { }
  
  bool JetDef::equivalentCluster(const JetDef& rhs) const {
    // essentially == operator, without jet selector comparison
    if (jet_algorithm() != rhs.jet_algorithm() ||
        R() != rhs.R() ||
        recombination_scheme() != rhs.recombination_scheme() ||
        strategy() != rhs.strategy() ||
        constituentSelector().description() != rhs.constituentSelector().description() ||
        areaDefinition().area_type() != rhs.areaDefinition().area_type() ||
        areaDefinition().ghost_spec().ghost_maxrap() != rhs.areaDefinition().ghost_spec().ghost_maxrap() ||
        areaDefinition().ghost_spec().ghost_area() != rhs.areaDefinition().ghost_spec().ghost_area() ||
        areaDefinition().ghost_spec().grid_scatter() != rhs.areaDefinition().ghost_spec().grid_scatter() ||
        areaDefinition().ghost_spec().pt_scatter() != rhs.areaDefinition().ghost_spec().pt_scatter() ||
        areaDefinition().ghost_spec().mean_ghost_pt() != rhs.areaDefinition().ghost_spec().mean_ghost_pt() ||
        areaDefinition().ghost_spec().repeat() != rhs.areaDefinition().ghost_spec().repeat() ||
        backgroundAreaDef().area_type() != rhs.backgroundAreaDef().area_type() ||
        backgroundAreaDef().ghost_spec().ghost_maxrap() != rhs.backgroundAreaDef().ghost_spec().ghost_maxrap() ||
        backgroundAreaDef().ghost_spec().ghost_area() != rhs.backgroundAreaDef().ghost_spec().ghost_area() ||
        backgroundAreaDef().ghost_spec().grid_scatter() != rhs.backgroundAreaDef().ghost_spec().grid_scatter() ||
        backgroundAreaDef().ghost_spec().pt_scatter() != rhs.backgroundAreaDef().ghost_spec().pt_scatter() ||
        backgroundAreaDef().ghost_spec().mean_ghost_pt() != rhs.backgroundAreaDef().ghost_spec().mean_ghost_pt() ||
        backgroundAreaDef().ghost_spec().repeat() != rhs.backgroundAreaDef().ghost_spec().repeat() ||
        backgroundJetDef().jet_algorithm() != rhs.backgroundJetDef().jet_algorithm() ||
        backgroundJetDef().R() != rhs.backgroundJetDef().R() ||
        backgroundJetDef().recombination_scheme() != rhs.backgroundJetDef().recombination_scheme() ||
        backgroundJetDef().strategy() != rhs.backgroundJetDef().strategy() ||
        backgroundSelector().description() != rhs.backgroundSelector().description())
      return false;
    return true;
  }

} // namespace dijetcore


