#include <iostream>
#include <math.h>
#include <string>

#include "jet_analysis/dijet_worker/dijet_matrix/jet_def.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/Selector.hh"

bool AreaDefEquality(fastjet::AreaDefinition a, fastjet::AreaDefinition b) {
  if (a.area_type() != b.area_type() ||
      a.ghost_spec().ghost_maxrap() != b.ghost_spec().ghost_maxrap() ||
      a.ghost_spec().repeat() != b.ghost_spec().repeat() ||
      a.ghost_spec().ghost_area() != b.ghost_spec().ghost_area() ||
      a.ghost_spec().grid_scatter() != b.ghost_spec().grid_scatter() ||
      a.ghost_spec().pt_scatter() != b.ghost_spec().pt_scatter() ||
      a.ghost_spec().mean_ghost_pt() != b.ghost_spec().mean_ghost_pt())
    return false;
  return true;
}

bool JetDefinitionEquality(fastjet::JetDefinition a, fastjet::JetDefinition b) {
  if (a.jet_algorithm() != b.jet_algorithm() ||
      a.R() != b.R() ||
      a.recombination_scheme() != b.recombination_scheme() ||
      a.extra_param() != b.extra_param() ||
      a.strategy() != b.strategy())
    return false;
  return true;
}

// selector comparison is done by descriptive strings,
// since I was too lazy to come up with a better way to
// compare compound selectors
bool SelectorEquality(fastjet::Selector a, fastjet::Selector b) {
  return a.description() == b.description();
}

bool JetDefEquality(JetDef a, JetDef b) {
  if (!JetDefinitionEquality(a, b) ||
      !JetDefinitionEquality(a.BackgroundJetDef(), b.BackgroundJetDef()) ||
      !SelectorEquality(a.ConstituentSelector(), b.ConstituentSelector()) ||
      !SelectorEquality(a.JetSelector(), b.JetSelector()) ||
      !SelectorEquality(a.BackgroundSelector(), b.BackgroundSelector()) ||
      !AreaDefEquality(a.AreaDefinition(), b.AreaDefinition()) ||
      !AreaDefEquality(a.BackgroundAreaDef(), b.BackgroundAreaDef()))
    return false;
  return true;
}

int main() {
  
  // check defaults
  JetDef default_def;
  
  if (JetDefinitionEquality(default_def, fastjet::JetDefinition()) != true ||
      JetDefinitionEquality(default_def.BackgroundJetDef(), fastjet::JetDefinition()) != true ||
      SelectorEquality(default_def.ConstituentSelector(), fastjet::SelectorIdentity()) != true ||
      SelectorEquality(default_def.JetSelector(), fastjet::SelectorIdentity()) != true ||
      SelectorEquality(default_def.BackgroundSelector(), (!fastjet::SelectorNHardest(2))) != true ||
      AreaDefEquality(default_def.AreaDefinition(), fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec())) != true)
    return 1;
  
  if (default_def.IsValid() == true ||
      default_def.CanEstimateArea() == true ||
      default_def.CanBackgroundSub() == true)
    return 1;
  
  // check non-default constructor
  JetDef simple_def(fastjet::antikt_algorithm, 1.0);
  default_def.set_jet_algorithm(fastjet::antikt_algorithm);
  if (JetDefEquality(simple_def, default_def) != true)
    return 1;
  
  if (simple_def.IsValid() != true ||
      simple_def.CanEstimateArea() != false ||
      simple_def.CanBackgroundSub() != false)
    return 1;
  
  fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
  fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
  fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);
  
  JetDef full_def(fastjet::antikt_algorithm, 1.0,
                  area_def,
                  bkg_def,
                  area_def,
                  bkg_sel);
  
  default_def.SetAreaDefinition(area_def);
  default_def.SetBackgroundJetDef(bkg_def);
  default_def.SetBackgroundSelector(bkg_sel);
  default_def.SetBackgroundAreaDef(area_def);
  
  if (JetDefEquality(full_def, default_def) != true)
    return 1;

  if (full_def.IsValid() != true ||
      full_def.CanEstimateArea() != true ||
      full_def.CanBackgroundSub() != true)
    return 1;
  
  // now test our EquivalentCluster algorithm
  // which should be false for any difference between
  // two jetdefs EXCEPT for jet selectors
  JetDef equivalent1;
  JetDef equivalent2;
  
  if (!equivalent1.EquivalentCluster(equivalent2))
    return 1;
  
  JetDef equivalent3(fastjet::antikt_algorithm, 0.4);
  JetDef equivalent4(fastjet::antikt_algorithm, 0.5);
  JetDef equivalent5(fastjet::kt_algorithm, 0.4);
  
  if (equivalent3.EquivalentCluster(equivalent4))
    return 1;
  if (equivalent3.EquivalentCluster(equivalent5))
    return 1;
  if (equivalent4.EquivalentCluster(equivalent5))
    return 1;
  
  return 0;
}

