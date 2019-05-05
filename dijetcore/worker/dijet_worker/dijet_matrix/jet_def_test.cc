#include "gtest/gtest.h"

#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"

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

bool JetDefEquality(dijetcore::JetDef a, dijetcore::JetDef b) {
  if (!JetDefinitionEquality(a, b) ||
      !JetDefinitionEquality(a.backgroundJetDef(), b.backgroundJetDef()) ||
      !SelectorEquality(a.constituentSelector(), b.constituentSelector()) ||
      !SelectorEquality(a.jetSelector(), b.jetSelector()) ||
      !SelectorEquality(a.backgroundSelector(), b.backgroundSelector()) ||
      !AreaDefEquality(a.areaDefinition(), b.areaDefinition()) ||
      !AreaDefEquality(a.backgroundAreaDef(), b.backgroundAreaDef()))
    return false;
  return true;
}


TEST(JetDef, DefaultSettings) {
    dijetcore::JetDef default_def;

    EXPECT_TRUE(JetDefinitionEquality(default_def, fastjet::JetDefinition()));
    EXPECT_TRUE(JetDefinitionEquality(default_def.backgroundJetDef(), fastjet::JetDefinition()));
    EXPECT_TRUE(SelectorEquality(default_def.constituentSelector(), fastjet::SelectorIdentity()));
    EXPECT_TRUE(SelectorEquality(default_def.backgroundSelector(), (!fastjet::SelectorNHardest(2))));
    EXPECT_TRUE(SelectorEquality(default_def.jetSelector(), fastjet::SelectorIdentity()));
    EXPECT_TRUE(AreaDefEquality(default_def.areaDefinition(), fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec())));

    EXPECT_FALSE(default_def.isValid());
    EXPECT_FALSE(default_def.canEstimateArea());
    EXPECT_FALSE(default_def.canBackgroundSub());
}

TEST(JetDef, SimpleDefinition) {
    dijetcore::JetDef simple_def(fastjet::antikt_algorithm, 1.0);
    dijetcore::JetDef default_def;
    default_def.set_jet_algorithm(fastjet::antikt_algorithm);
    EXPECT_TRUE(JetDefEquality(simple_def, default_def));
    
    EXPECT_TRUE(simple_def.isValid());
    EXPECT_FALSE(simple_def.canEstimateArea());
    EXPECT_FALSE(simple_def.canBackgroundSub());
}

TEST(JetDef, FullDefinition) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);

    dijetcore::JetDef full_def(fastjet::antikt_algorithm, 1.0, area_def, bkg_def, area_def, bkg_sel);

    dijetcore::JetDef default_def;
    default_def.set_jet_algorithm(fastjet::antikt_algorithm);
    default_def.setAreaDefinition(area_def);
    default_def.setBackgroundJetDef(bkg_def);
    default_def.setBackgroundSelector(bkg_sel);
    default_def.setBackgroundAreaDef(area_def);

    EXPECT_TRUE(JetDefEquality(full_def, default_def));

    EXPECT_TRUE(full_def.isValid());
    EXPECT_TRUE(full_def.canEstimateArea());
    EXPECT_TRUE(full_def.canBackgroundSub());
}

TEST(JetDef, EquivalantClustering) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);

    dijetcore::JetDef full_def(fastjet::antikt_algorithm, 1.0, area_def, bkg_def, area_def, bkg_sel);

    dijetcore::JetDef default_def;
    default_def.set_jet_algorithm(fastjet::antikt_algorithm);
    default_def.setAreaDefinition(area_def);
    default_def.setBackgroundJetDef(bkg_def);
    default_def.setBackgroundSelector(bkg_sel);
    default_def.setBackgroundAreaDef(area_def);

    EXPECT_TRUE(full_def.equivalentCluster(default_def));
}
