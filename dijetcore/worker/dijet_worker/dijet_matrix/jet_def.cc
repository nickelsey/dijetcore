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
  
  bool JetDef::EquivalentCluster(const JetDef& rhs) const {
    // essentially == operator, without jet selector comparison
    if (jet_algorithm() != rhs.jet_algorithm() ||
        R() != rhs.R() ||
        recombination_scheme() != rhs.recombination_scheme() ||
        strategy() != rhs.strategy() ||
        ConstituentSelector().description() != rhs.ConstituentSelector().description() ||
        AreaDefinition().area_type() != rhs.AreaDefinition().area_type() ||
        AreaDefinition().ghost_spec().ghost_maxrap() != rhs.AreaDefinition().ghost_spec().ghost_maxrap() ||
        AreaDefinition().ghost_spec().ghost_area() != rhs.AreaDefinition().ghost_spec().ghost_area() ||
        AreaDefinition().ghost_spec().grid_scatter() != rhs.AreaDefinition().ghost_spec().grid_scatter() ||
        AreaDefinition().ghost_spec().pt_scatter() != rhs.AreaDefinition().ghost_spec().pt_scatter() ||
        AreaDefinition().ghost_spec().mean_ghost_pt() != rhs.AreaDefinition().ghost_spec().mean_ghost_pt() ||
        AreaDefinition().ghost_spec().repeat() != rhs.AreaDefinition().ghost_spec().repeat() ||
        BackgroundAreaDef().area_type() != rhs.BackgroundAreaDef().area_type() ||
        BackgroundAreaDef().ghost_spec().ghost_maxrap() != rhs.BackgroundAreaDef().ghost_spec().ghost_maxrap() ||
        BackgroundAreaDef().ghost_spec().ghost_area() != rhs.BackgroundAreaDef().ghost_spec().ghost_area() ||
        BackgroundAreaDef().ghost_spec().grid_scatter() != rhs.BackgroundAreaDef().ghost_spec().grid_scatter() ||
        BackgroundAreaDef().ghost_spec().pt_scatter() != rhs.BackgroundAreaDef().ghost_spec().pt_scatter() ||
        BackgroundAreaDef().ghost_spec().mean_ghost_pt() != rhs.BackgroundAreaDef().ghost_spec().mean_ghost_pt() ||
        BackgroundAreaDef().ghost_spec().repeat() != rhs.BackgroundAreaDef().ghost_spec().repeat() ||
        BackgroundJetDef().jet_algorithm() != rhs.BackgroundJetDef().jet_algorithm() ||
        BackgroundJetDef().R() != rhs.BackgroundJetDef().R() ||
        BackgroundJetDef().recombination_scheme() != rhs.BackgroundJetDef().recombination_scheme() ||
        BackgroundJetDef().strategy() != rhs.BackgroundJetDef().strategy() ||
        BackgroundSelector().description() != rhs.BackgroundSelector().description())
      return false;
    return true;
  }

} // namespace dijetcore


