#include "jet_analysis/util/dijet_key.hh"
#include "jet_analysis/dijet_worker/dijet_matrix/jet_def.hh"
#include "jet_analysis/dijet_worker/dijet_matrix/match_def.hh"
#include "jet_analysis/dijet_worker/dijet_matrix/dijet_definition.hh"

#include <iostream>

int main() {
  
  std::string expected = "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.5_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.6_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.7_alg_2_pt_0_const_eta_1_const_pt_0.2";
  
  fastjet::Selector const_sel = fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  fastjet::Selector const_sel_2 = fastjet::SelectorPtMin(2.0) && fastjet::SelectorAbsRapMax(1.0);
  fastjet::Selector jet_sel = fastjet::SelectorPtMin(10.0);
  fastjet::Selector jet_sel_2 = fastjet::SelectorPtMin(20.0);
  
  
  JetDef def_A(fastjet::antikt_algorithm, 0.4);
  def_A.SetConstituentSelector(const_sel_2);
  def_A.SetJetSelector(jet_sel_2);
  JetDef def_B(fastjet::antikt_algorithm, 0.5);
  def_B.SetConstituentSelector(const_sel);
  JetDef def_C(fastjet::antikt_algorithm, 0.6);
  def_C.SetConstituentSelector(const_sel_2);
  def_C.SetJetSelector(jet_sel_2);
  JetDef def_D(fastjet::antikt_algorithm, 0.7);
  def_D.SetConstituentSelector(const_sel);
  
  MatchDef* match_A = new MatchDef(def_A, def_B, 0.4);
  MatchDef* match_B = new MatchDef(def_C, def_D, 0.6);
  
  DijetDefinition dijet_A(match_A, match_B);
  
  std::string key = MakeKeyFromDijetDefinition(dijet_A);
  
  if (key != expected)
    return 1;
  
  return 0;
}
