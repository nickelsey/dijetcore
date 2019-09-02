#include "dijetcore/test/dijetcore_test_helper.h"

namespace dijetcore {
namespace testing {

void CheckJetDef(dijetcore::JetDef def, fastjet::JetAlgorithm alg,
                 double radius, double jet_pt, double const_pt,
                 double const_eta, fastjet::RecombinationScheme scheme,
                 fastjet::Strategy strategy, fastjet::AreaType area_type,
                 int ghost_repeat, double ghost_area, double grid_scatter,
                 double pt_scatter, double mean_ghost_pt,
                 fastjet::JetDefinition bkg_def,
                 fastjet::AreaDefinition bkg_area) {
  EXPECT_EQ(def.R(), radius);
  EXPECT_EQ(def.jet_algorithm(), alg);
  EXPECT_EQ(
      dijetcore::ExtractDoubleFromSelector(def.constituentSelector(), "pt >="),
      const_pt);
  if (jet_pt > 0.0)
    EXPECT_EQ(dijetcore::ExtractDoubleFromSelector(def.jetSelector(), "pt >="),
              jet_pt);
  EXPECT_EQ(def.recombination_scheme(), scheme);
  EXPECT_EQ(def.strategy(), strategy);
  EXPECT_EQ(def.areaDefinition().area_type(), area_type);
  EXPECT_EQ(def.areaDefinition().ghost_spec().ghost_maxrap(),
            const_eta + radius);
  EXPECT_EQ(def.areaDefinition().ghost_spec().repeat(), ghost_repeat);
  EXPECT_EQ(def.areaDefinition().ghost_spec().ghost_area(), ghost_area);
  EXPECT_EQ(def.areaDefinition().ghost_spec().grid_scatter(), grid_scatter);
  EXPECT_EQ(def.areaDefinition().ghost_spec().pt_scatter(), pt_scatter);
  EXPECT_EQ(def.areaDefinition().ghost_spec().mean_ghost_pt(), mean_ghost_pt);
  EXPECT_EQ(def.backgroundJetDef().R(), bkg_def.R());
  EXPECT_EQ(def.backgroundJetDef().jet_algorithm(), bkg_def.jet_algorithm());
  EXPECT_EQ(def.backgroundJetDef().strategy(), bkg_def.strategy());
  EXPECT_EQ(def.backgroundJetDef().recombination_scheme(),
            bkg_def.recombination_scheme());
  EXPECT_EQ(def.backgroundAreaDef().area_type(), bkg_area.area_type());
  EXPECT_EQ(def.backgroundAreaDef().ghost_spec().ghost_maxrap(),
            bkg_area.ghost_spec().ghost_maxrap());
  EXPECT_EQ(def.backgroundAreaDef().ghost_spec().repeat(),
            bkg_area.ghost_spec().repeat());
  EXPECT_EQ(def.backgroundAreaDef().ghost_spec().ghost_area(),
            bkg_area.ghost_spec().ghost_area());
  EXPECT_EQ(def.backgroundAreaDef().ghost_spec().grid_scatter(),
            bkg_area.ghost_spec().grid_scatter());
  EXPECT_EQ(def.backgroundAreaDef().ghost_spec().pt_scatter(),
            bkg_area.ghost_spec().pt_scatter());
  EXPECT_EQ(def.backgroundAreaDef().ghost_spec().mean_ghost_pt(),
            bkg_area.ghost_spec().mean_ghost_pt());
  return;
}

void CheckMatchDef(dijetcore::MatchDef *def, fastjet::JetAlgorithm alg,
                   double radius, double match_radius, double jet_pt,
                   double const_init_pt, double const_match_pt,
                   double const_eta, fastjet::RecombinationScheme scheme,
                   fastjet::Strategy strategy, fastjet::AreaType area_type,
                   int ghost_repeat, double ghost_area, double grid_scatter,
                   double pt_scatter, double mean_ghost_pt,
                   fastjet::JetDefinition bkg_def,
                   fastjet::AreaDefinition bkg_area) {

  CheckJetDef(def->initialJetDef(), alg, radius, jet_pt, const_init_pt,
              const_eta, scheme, strategy, area_type, ghost_repeat, ghost_area,
              grid_scatter, pt_scatter, mean_ghost_pt, bkg_def, bkg_area);

  CheckJetDef(def->matchedJetDef(), alg, match_radius, 0.0, const_match_pt,
              const_eta, scheme, strategy, area_type, ghost_repeat, ghost_area,
              grid_scatter, pt_scatter, mean_ghost_pt, bkg_def, bkg_area);
  return;
}

void CheckDijetDefinition(
    dijetcore::DijetDefinition *def, fastjet::JetAlgorithm lead_alg,
    fastjet::JetAlgorithm sub_alg, double lead_R, double lead_match_R,
    double sub_R, double sub_match_R, double lead_pt, double sub_pt,
    double lead_const_init_pt, double lead_const_match_pt,
    double sub_const_init_pt, double sub_const_match_pt, double const_eta,
    fastjet::RecombinationScheme scheme, fastjet::Strategy strategy,
    fastjet::AreaType area_type, int ghost_repeat, double ghost_area,
    double grid_scatter, double pt_scatter, double mean_ghost_pt,
    fastjet::JetDefinition bkg_def, fastjet::AreaDefinition bkg_area_lead,
    fastjet::AreaDefinition bkg_area_sub) {

  CheckMatchDef(def->lead, lead_alg, lead_R, lead_match_R, lead_pt,
                lead_const_init_pt, lead_const_match_pt, const_eta, scheme,
                strategy, area_type, ghost_repeat, ghost_area, grid_scatter,
                pt_scatter, mean_ghost_pt, bkg_def, bkg_area_lead);
  CheckMatchDef(def->sub, sub_alg, sub_R, sub_match_R, sub_pt,
                sub_const_init_pt, sub_const_match_pt, const_eta, scheme,
                strategy, area_type, ghost_repeat, ghost_area, grid_scatter,
                pt_scatter, mean_ghost_pt, bkg_def, bkg_area_sub);
  return;
}

} // namespace testing
} // namespace dijetcore