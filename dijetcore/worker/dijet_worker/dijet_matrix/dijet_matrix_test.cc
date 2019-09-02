#include "gtest/gtest.h"

#include <string>

#include "dijetcore/lib/logging.h"
#include "dijetcore/test/dijetcore_test_helper.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/fastjet/selector_compare.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_matrix.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"

TEST(DijetMatrix, Default) {
  dijetcore::DijetMatrix default_matrix;

  default_matrix.initialize();
  std::string expected_default_string =
      "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_"
      "0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_"
      "pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2";

  std::set<std::string> keys = default_matrix.keys();

  EXPECT_NE(keys.find(expected_default_string), keys.end());
}

TEST(DijetMatrix, Initialization) {
  dijetcore::DijetMatrix matrix;
  matrix.forceConstituentPtEquality(false);
  matrix.forceConstituentEtaEquality(false);
  matrix.forceJetResolutionEquality(false);
  matrix.addJetAlgorithm({fastjet::antikt_algorithm, fastjet::kt_algorithm});
  matrix.addLeadJetR({0.4, 0.5});
  matrix.addSubJetR({0.4, 0.5});
  matrix.addLeadJetPt({20, 18});
  matrix.addSubJetPt({10, 9});
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 32);

  matrix.addSubJetPt(40);
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 32);

  matrix.addLeadJetPt({2, 60});
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 48 + 8);

  matrix.clear();
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 1);
}

TEST(DijetMatrix, ForceConstistuentEta) {
  dijetcore::DijetMatrix matrix;
  matrix.addConstituentEta({1.0, 1.5});
  matrix.forceConstituentEtaEquality(true);
  matrix.forceJetResolutionEquality(false);
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 2);

  matrix.forceConstituentEtaEquality(false);
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 4);
}

TEST(DijetMatrix, ForceConstituentPt) {
  dijetcore::DijetMatrix matrix;

  matrix.forceConstituentPtEquality(true);
  matrix.forceJetResolutionEquality(false);

  matrix.addConstituentLeadInitialPt(2.0);
  matrix.addConstituentSubInitialPt(2.1);

  matrix.initialize();
  EXPECT_EQ(matrix.size(), 0);

  matrix.forceConstituentPtEquality(false);

  matrix.initialize();
  EXPECT_EQ(matrix.size(), 1);

  matrix.clear();

  matrix.addConstituentLeadInitialPt({2.0, 2.2});
  matrix.addConstituentSubInitialPt({2.1, 2.2});
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 4);

  matrix.forceConstituentPtEquality(true);
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 1);

  matrix.addConstituentLeadMatchPt({0.2, 0.3});
  matrix.addConstituentSubMatchPt(0.2);
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 1);

  matrix.forceConstituentPtEquality(false);
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 8);
}

TEST(DijetMatrix, SingleCaseFullTest) {
  dijetcore::DijetMatrix matrix;
  matrix.forceConstituentPtEquality(true);
  matrix.forceConstituentEtaEquality(true);
  matrix.forceJetResolutionEquality(false);
  matrix.forceMatchJetResolutionEquality(true);
  matrix.addLeadJetPt(18);
  matrix.addLeadJetR(0.4);
  matrix.addLeadJetR(0.5);
  matrix.addSubJetR(0.4);
  matrix.addSubJetR(0.6);
  matrix.addSubJetPt(9);
  matrix.addConstituentEta(1.0);
  matrix.addConstituentLeadInitialPt(2.5);
  matrix.addConstituentLeadMatchPt(0.5);
  matrix.addConstituentSubInitialPt(2.3);
  matrix.addConstituentSubMatchPt(0.4);
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 0);

  matrix.forceConstituentPtEquality(false);
  matrix.forceConstituentEtaEquality(false);
  matrix.forceMatchJetResolutionEquality(false);
  matrix.initialize();

  EXPECT_EQ(matrix.size(), 4);

  std::vector<std::pair<double, double>> lead_jet_pairs{{0.4, 0.4}, {0.5, 0.4}};
  std::vector<std::pair<double, double>> sub_jet_pairs{{0.4, 0.4}, {0.6, 0.4}};

  for (auto key : matrix.keys()) {
    dijetcore::DijetKey parsed_key = dijetcore::ParseStringToDijetKey(key);

    if (std::find(lead_jet_pairs.begin(), lead_jet_pairs.end(),
                  std::pair<double, double>{parsed_key.lead_init_r,
                                            parsed_key.lead_match_r}) ==
        lead_jet_pairs.end()) {
      EXPECT_STREQ("Found bad DijetDefinition in DijetMatrix", "");
      continue;
    }
    double lead_R = parsed_key.lead_init_r;
    double lead_match_R = parsed_key.lead_match_r;

    if (std::find(sub_jet_pairs.begin(), sub_jet_pairs.end(),
                  std::pair<double, double>{parsed_key.sub_init_r,
                                            parsed_key.sub_match_r}) ==
        sub_jet_pairs.end()) {
      EXPECT_STREQ("Found bad DijetDefinition in DijetMatrix", "");
      continue;
    }
    double sub_R = parsed_key.sub_init_r;
    double sub_match_R = parsed_key.sub_match_r;
    // build the background area def

    fastjet::GhostedAreaSpec ghost_spec_lead(1.0 + 0.4, 1, 0.01, 1, 0.1,
                                             1e-100);
    fastjet::AreaDefinition bkg_area_def_lead(
        fastjet::active_area_explicit_ghosts, ghost_spec_lead);

    fastjet::GhostedAreaSpec ghost_spec_sub(1.0 + 0.4, 1, 0.01, 1, 0.1, 1e-100);
    fastjet::AreaDefinition bkg_area_def_sub(
        fastjet::active_area_explicit_ghosts, ghost_spec_sub);

    dijetcore::testing::CheckDijetDefinition(
        matrix.dijetDefinitions()[key].get(), fastjet::antikt_algorithm,
        fastjet::antikt_algorithm, lead_R, lead_match_R, sub_R, sub_match_R, 18,
        9, 2.5, 0.5, 2.3, 0.4, 1.0, fastjet::E_scheme, fastjet::Best,
        fastjet::active_area_explicit_ghosts, 1, 0.01, 1, 0.1, 1e-100,
        fastjet::JetDefinition(fastjet::kt_algorithm, 0.4), bkg_area_def_lead,
        bkg_area_def_sub);
  }
}

TEST(DijetMatrix, TestClusterSequenceSets) {
  dijetcore::DijetMatrix matrix;
  matrix.forceJetResolutionEquality(false);
  matrix.addLeadJetPt({20.0, 18.0});
  matrix.addSubJetPt({10.0, 9.0});
  matrix.addLeadJetR(0.4);
  matrix.addSubJetR(0.4);
  matrix.addConstituentEta(0.4);
  matrix.addConstituentLeadInitialPt(2.0);
  matrix.addConstituentSubInitialPt(2.0);
  matrix.addConstituentLeadMatchPt(0.2);
  matrix.addConstituentSubMatchPt(0.2);
  matrix.initialize();

  ASSERT_EQ(matrix.sortedDefinitions().size(), 1);
  EXPECT_EQ(matrix.sortedDefinitions()[0].size(), 4);

  matrix.addLeadJetR(0.5);
  matrix.addSubJetR(0.5);
  matrix.initialize();

  ASSERT_EQ(matrix.sortedDefinitions().size(), 4);

  for (int i = 0; i < matrix.sortedDefinitions().size(); ++i) {
    EXPECT_EQ(matrix.sortedDefinitions()[i].size(), 4);

    auto &defs = matrix.sortedDefinitions()[i];

    std::set<double> lead_init_r;
    std::set<double> lead_match_r;
    std::set<double> sub_init_r;
    std::set<double> sub_match_r;

    for (auto &e : defs) {
      lead_init_r.insert(e.second->lead->initialJetDef().R());
      lead_match_r.insert(e.second->lead->matchedJetDef().R());
      sub_init_r.insert(e.second->sub->initialJetDef().R());
      sub_match_r.insert(e.second->sub->matchedJetDef().R());
    }

    EXPECT_EQ(lead_init_r.size(), 1);
    EXPECT_EQ(lead_match_r.size(), 1);
    EXPECT_EQ(sub_init_r.size(), 1);
    EXPECT_EQ(sub_match_r.size(), 1);
  }

  matrix.addConstituentLeadInitialPt({1.0, 3.0});
  matrix.addConstituentLeadMatchPt({0.3, 0.4});
  matrix.addConstituentSubInitialPt({1.0, 3.0});
  matrix.addConstituentSubMatchPt({0.3, 0.4});
  matrix.initialize();

  for (int i = 0; i < matrix.sortedDefinitions().size(); ++i) {
    EXPECT_EQ(matrix.sortedDefinitions()[i].size(), 4);

    auto &defs = matrix.sortedDefinitions()[i];

    std::set<double> lead_init_r;
    std::set<double> lead_match_r;
    std::set<double> sub_init_r;
    std::set<double> sub_match_r;
    std::set<double> lead_init_const_pt;
    std::set<double> lead_match_const_pt;
    std::set<double> sub_init_const_pt;
    std::set<double> sub_match_const_pt;

    for (auto &e : defs) {
      lead_init_r.insert(e.second->lead->initialJetDef().R());
      lead_match_r.insert(e.second->lead->matchedJetDef().R());
      sub_init_r.insert(e.second->sub->initialJetDef().R());
      sub_match_r.insert(e.second->sub->matchedJetDef().R());
      lead_init_const_pt.insert(dijetcore::ExtractDoubleFromSelector(
          e.second->lead->initialJetDef().constituentSelector(), "pt >="));
      lead_match_const_pt.insert(dijetcore::ExtractDoubleFromSelector(
          e.second->lead->matchedJetDef().constituentSelector(), "pt >="));
      sub_init_const_pt.insert(dijetcore::ExtractDoubleFromSelector(
          e.second->sub->initialJetDef().constituentSelector(), "pt >="));
      sub_match_const_pt.insert(dijetcore::ExtractDoubleFromSelector(
          e.second->sub->matchedJetDef().constituentSelector(), "pt >="));
    }

    EXPECT_EQ(lead_init_r.size(), 1);
    EXPECT_EQ(lead_match_r.size(), 1);
    EXPECT_EQ(sub_init_r.size(), 1);
    EXPECT_EQ(sub_match_r.size(), 1);
    EXPECT_EQ(lead_init_const_pt.size(), 1);
    EXPECT_EQ(lead_match_const_pt.size(), 1);
    EXPECT_EQ(sub_init_const_pt.size(), 1);
    EXPECT_EQ(sub_match_const_pt.size(), 1);
  }
}
