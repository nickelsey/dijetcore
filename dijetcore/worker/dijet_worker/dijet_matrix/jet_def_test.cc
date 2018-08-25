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
      !JetDefinitionEquality(a.BackgroundJetDef(), b.BackgroundJetDef()) ||
      !SelectorEquality(a.ConstituentSelector(), b.ConstituentSelector()) ||
      !SelectorEquality(a.JetSelector(), b.JetSelector()) ||
      !SelectorEquality(a.BackgroundSelector(), b.BackgroundSelector()) ||
      !AreaDefEquality(a.AreaDefinition(), b.AreaDefinition()) ||
      !AreaDefEquality(a.BackgroundAreaDef(), b.BackgroundAreaDef()))
    return false;
  return true;
}


TEST(JetDef, DefaultSettings) {
    dijetcore::JetDef default_def;

    EXPECT_TRUE(JetDefinitionEquality(default_def, fastjet::JetDefinition()));
    EXPECT_TRUE(JetDefinitionEquality(default_def.BackgroundJetDef(), fastjet::JetDefinition()));
    EXPECT_TRUE(SelectorEquality(default_def.ConstituentSelector(), fastjet::SelectorIdentity()));
    EXPECT_TRUE(SelectorEquality(default_def.BackgroundSelector(), (!fastjet::SelectorNHardest(2))));
    EXPECT_TRUE(SelectorEquality(default_def.JetSelector(), fastjet::SelectorIdentity()));
    EXPECT_TRUE(AreaDefEquality(default_def.AreaDefinition(), fastjet::AreaDefinition(fastjet::invalid_area, fastjet::GhostedAreaSpec())));

    EXPECT_FALSE(default_def.IsValid());
    EXPECT_FALSE(default_def.CanEstimateArea());
    EXPECT_FALSE(default_def.CanBackgroundSub());
}

TEST(JetDef, SimpleDefinition) {
    dijetcore::JetDef simple_def(fastjet::antikt_algorithm, 1.0);
    dijetcore::JetDef default_def;
    default_def.set_jet_algorithm(fastjet::antikt_algorithm);
    EXPECT_TRUE(JetDefEquality(simple_def, default_def));
    
    EXPECT_TRUE(simple_def.IsValid());
    EXPECT_FALSE(simple_def.CanEstimateArea());
    EXPECT_FALSE(simple_def.CanBackgroundSub());
}

TEST(JetDef, FullDefinition) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);

    dijetcore::JetDef full_def(fastjet::antikt_algorithm, 1.0, area_def, bkg_def, area_def, bkg_sel);

    dijetcore::JetDef default_def;
    default_def.set_jet_algorithm(fastjet::antikt_algorithm);
    default_def.SetAreaDefinition(area_def);
    default_def.SetBackgroundJetDef(bkg_def);
    default_def.SetBackgroundSelector(bkg_sel);
    default_def.SetBackgroundAreaDef(area_def);

    EXPECT_TRUE(JetDefEquality(full_def, default_def));

    EXPECT_TRUE(full_def.IsValid());
    EXPECT_TRUE(full_def.CanEstimateArea());
    EXPECT_TRUE(full_def.CanBackgroundSub());
}

TEST(JetDef, EquivalantClustering) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);

    dijetcore::JetDef full_def(fastjet::antikt_algorithm, 1.0, area_def, bkg_def, area_def, bkg_sel);

    dijetcore::JetDef default_def;
    default_def.set_jet_algorithm(fastjet::antikt_algorithm);
    default_def.SetAreaDefinition(area_def);
    default_def.SetBackgroundJetDef(bkg_def);
    default_def.SetBackgroundSelector(bkg_sel);
    default_def.SetBackgroundAreaDef(area_def);

    EXPECT_TRUE(full_def.EquivalentCluster(default_def));
}
