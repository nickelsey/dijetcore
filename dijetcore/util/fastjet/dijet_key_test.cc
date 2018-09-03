#include "gtest/gtest.h"

#include <string>

#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

using std::string;

TEST(DijetKey, CompleteKey) {
  
  string expected = "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.5_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.6_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.7_alg_2_pt_0_const_eta_1_const_pt_0.2";
  
  fastjet::Selector const_sel = fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  fastjet::Selector const_sel_2 = fastjet::SelectorPtMin(2.0) && fastjet::SelectorAbsRapMax(1.0);
  fastjet::Selector jet_sel = fastjet::SelectorPtMin(10.0);
  fastjet::Selector jet_sel_2 = fastjet::SelectorPtMin(20.0);
  
  
  dijetcore::JetDef def_A(fastjet::antikt_algorithm, 0.4);
  def_A.SetConstituentSelector(const_sel_2);
  def_A.SetJetSelector(jet_sel_2);
  dijetcore::JetDef def_B(fastjet::antikt_algorithm, 0.5);
  def_B.SetConstituentSelector(const_sel);
  dijetcore::JetDef def_C(fastjet::antikt_algorithm, 0.6);
  def_C.SetConstituentSelector(const_sel_2);
  def_C.SetJetSelector(jet_sel_2);
  dijetcore::JetDef def_D(fastjet::antikt_algorithm, 0.7);
  def_D.SetConstituentSelector(const_sel);
  
  dijetcore::MatchDef* match_A = new dijetcore::MatchDef(def_A, def_B, 0.4);
  dijetcore::MatchDef* match_B = new dijetcore::MatchDef(def_C, def_D, 0.6);
  
  dijetcore::DijetDefinition dijet_A(match_A, match_B);
  
  string key = dijetcore::MakeKeyFromDijetDefinition(dijet_A);
  
    EXPECT_EQ(key, expected);
}

TEST(DijetKey, Conversion) {
  string expected = "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.5_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.6_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.7_alg_2_pt_0_const_eta_1_const_pt_0.2";
  dijetcore::DijetKey key = dijetcore::ParseStringToDijetKey(expected);
  string output = dijetcore::ParseDijetKeyToString(key);
  
  EXPECT_EQ(expected, output);
  
  EXPECT_EQ(key.lead_init_r, 0.4);
  EXPECT_EQ(key.lead_init_alg, 2);
  EXPECT_EQ(key.lead_init_pt, 20);
  EXPECT_EQ(key.lead_match_r, 0.5);
  EXPECT_EQ(key.lead_match_alg, 2);
  EXPECT_EQ(key.lead_match_pt, 0);
  EXPECT_EQ(key.lead_init_const_eta, 1);
  EXPECT_EQ(key.lead_init_const_pt, 2);
  EXPECT_EQ(key.lead_match_const_eta, 1);
  EXPECT_EQ(key.lead_match_const_pt, 0.2);
  
  EXPECT_EQ(key.sub_init_r, 0.6);
  EXPECT_EQ(key.sub_init_alg, 2);
  EXPECT_EQ(key.sub_init_pt, 20);
  EXPECT_EQ(key.sub_match_r, 0.7);
  EXPECT_EQ(key.sub_match_alg, 2);
  EXPECT_EQ(key.sub_match_pt, 0);
  EXPECT_EQ(key.sub_init_const_eta, 1);
  EXPECT_EQ(key.sub_init_const_pt, 2);
  EXPECT_EQ(key.sub_match_const_eta, 1);
  EXPECT_EQ(key.sub_match_const_pt, 0.2);
}
