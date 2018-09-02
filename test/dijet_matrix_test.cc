#include "jet_analysis/dijet_worker/dijet_matrix/dijet_matrix.hh"
#include "jet_analysis/util/selector_compare.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"

#include <string>
#include <set>

bool CheckDijetDefinition(DijetDefinition def, fastjet::JetAlgorithm lead_alg, fastjet::JetAlgorithm sub_alg,
                          double lead_R, double sub_R, double lead_pt, double sub_pt, double lead_const_init_pt,
                          double lead_const_match_pt, double sub_const_init_pt, double sub_const_match_pt,
                          double const_eta, fastjet::RecombinationScheme scheme, fastjet::Strategy strategy,
                          fastjet::AreaType area_type, int ghost_repeat, double ghost_area, double grid_scatter,
                          double pt_scatter, double mean_ghost_pt, fastjet::JetDefinition bkg_def,
                          fastjet::AreaDefinition bkg_area_lead, fastjet::AreaDefinition bkg_area_sub) {
  MatchDef* lead = def.lead;
  MatchDef* sub = def.sub;
  
  // check the leading jet
  if (lead->InitialJetDef().R() != lead_R ||
      lead->MatchedJetDef().R() != lead_R ||
      lead->InitialJetDef().jet_algorithm() != lead_alg ||
      lead->MatchedJetDef().jet_algorithm() != lead_alg ||
      ExtractDoubleFromSelector(lead->InitialJetDef().ConstituentSelector(), "pt >=") != lead_const_init_pt ||
      ExtractDoubleFromSelector(lead->MatchedJetDef().ConstituentSelector(), "pt >=") != lead_const_match_pt ||
      ExtractDoubleFromSelector(lead->InitialJetDef().JetSelector(), "pt >=") != lead_pt ||
      ExtractDoubleFromSelector(lead->MatchedJetDef().JetSelector(), "pt >=") != 0 ||
      lead->InitialJetDef().recombination_scheme() != scheme ||
      lead->MatchedJetDef().recombination_scheme() != scheme ||
      lead->InitialJetDef().strategy() != strategy ||
      lead->MatchedJetDef().strategy() != strategy ||
      lead->InitialJetDef().AreaDefinition().area_type() != area_type ||
      lead->MatchedJetDef().AreaDefinition().area_type() != area_type ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().ghost_maxrap() != (const_eta + lead_R) ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().ghost_maxrap() != (const_eta + lead_R) ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().repeat() != ghost_repeat ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().repeat() != ghost_repeat ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().ghost_area() != ghost_area ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().ghost_area() != ghost_area ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      lead->InitialJetDef().BackgroundJetDef().R() != bkg_def.R() ||
      lead->MatchedJetDef().BackgroundJetDef().R() != bkg_def.R() ||
      lead->InitialJetDef().BackgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      lead->MatchedJetDef().BackgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      lead->InitialJetDef().BackgroundJetDef().strategy() != bkg_def.strategy() ||
      lead->MatchedJetDef().BackgroundJetDef().strategy() != bkg_def.strategy() ||
      lead->InitialJetDef().BackgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      lead->MatchedJetDef().BackgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      lead->InitialJetDef().BackgroundAreaDef().area_type() != bkg_area_lead.area_type() ||
      lead->MatchedJetDef().BackgroundAreaDef().area_type() != bkg_area_lead.area_type() ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_lead.ghost_spec().ghost_maxrap() ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_lead.ghost_spec().ghost_maxrap() ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().repeat() != bkg_area_lead.ghost_spec().repeat() ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().repeat() != bkg_area_lead.ghost_spec().repeat() ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().ghost_area() != bkg_area_lead.ghost_spec().ghost_area() ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().ghost_area() != bkg_area_lead.ghost_spec().ghost_area() ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_lead.ghost_spec().grid_scatter()  ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_lead.ghost_spec().grid_scatter()  ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_lead.ghost_spec().pt_scatter()  ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_lead.ghost_spec().pt_scatter()  ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_lead.ghost_spec().mean_ghost_pt()  ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_lead.ghost_spec().mean_ghost_pt() )
    return false;
  
  // and check subleading jet
  if (sub->InitialJetDef().R() != sub_R ||
      sub->MatchedJetDef().R() != sub_R ||
      sub->InitialJetDef().jet_algorithm() != sub_alg ||
      sub->MatchedJetDef().jet_algorithm() != sub_alg ||
      ExtractDoubleFromSelector(sub->InitialJetDef().ConstituentSelector(), "pt >=") != sub_const_init_pt ||
      ExtractDoubleFromSelector(sub->MatchedJetDef().ConstituentSelector(), "pt >=") != sub_const_match_pt ||
      ExtractDoubleFromSelector(sub->InitialJetDef().JetSelector(), "pt >=") != sub_pt ||
      ExtractDoubleFromSelector(sub->MatchedJetDef().JetSelector(), "pt >=") != 0 ||
      sub->InitialJetDef().recombination_scheme() != scheme ||
      sub->MatchedJetDef().recombination_scheme() != scheme ||
      sub->InitialJetDef().strategy() != strategy ||
      sub->MatchedJetDef().strategy() != strategy ||
      sub->InitialJetDef().AreaDefinition().area_type() != area_type ||
      sub->MatchedJetDef().AreaDefinition().area_type() != area_type ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().ghost_maxrap() != (const_eta + sub_R) ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().ghost_maxrap() != (const_eta + sub_R) ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().repeat() != ghost_repeat ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().repeat() != ghost_repeat ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().ghost_area() != ghost_area ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().ghost_area() != ghost_area ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      sub->InitialJetDef().BackgroundJetDef().R() != bkg_def.R() ||
      sub->MatchedJetDef().BackgroundJetDef().R() != bkg_def.R() ||
      sub->InitialJetDef().BackgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      sub->MatchedJetDef().BackgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      sub->InitialJetDef().BackgroundJetDef().strategy() != bkg_def.strategy() ||
      sub->MatchedJetDef().BackgroundJetDef().strategy() != bkg_def.strategy() ||
      sub->InitialJetDef().BackgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      sub->MatchedJetDef().BackgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      sub->InitialJetDef().BackgroundAreaDef().area_type() != bkg_area_sub.area_type() ||
      sub->MatchedJetDef().BackgroundAreaDef().area_type() != bkg_area_sub.area_type() ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_sub.ghost_spec().ghost_maxrap() ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_sub.ghost_spec().ghost_maxrap() ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().repeat() != bkg_area_sub.ghost_spec().repeat() ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().repeat() != bkg_area_sub.ghost_spec().repeat() ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().ghost_area() != bkg_area_sub.ghost_spec().ghost_area() ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().ghost_area() != bkg_area_sub.ghost_spec().ghost_area() ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_sub.ghost_spec().grid_scatter()  ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_sub.ghost_spec().grid_scatter()  ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_sub.ghost_spec().pt_scatter()  ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_sub.ghost_spec().pt_scatter()  ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_sub.ghost_spec().mean_ghost_pt()  ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_sub.ghost_spec().mean_ghost_pt() )
    return false;
  
  return true;
}

int main() {
  
  DijetMatrix default_matrix;
  
  default_matrix.Initialize();
  std::string expected_default_string = "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2";
  
  std::set<std::string> keys = default_matrix.Keys();
  
  if (keys.find(expected_default_string) == keys.end())
    return 1;
  
  try {
    default_matrix.DijetDefinitions().at(expected_default_string);
  } catch(...) {
    return 1;
  }
  
  DijetMatrix matrix;
  matrix.ForceConstituentPtEquality(false);
  matrix.ForceConstituentEtaEquality(false);
  matrix.AddJetAlgorithm({fastjet::antikt_algorithm, fastjet::kt_algorithm});
  matrix.AddLeadJetR({0.4, 0.5});
  matrix.AddSubJetR({0.4, 0.5});
  matrix.AddLeadJetPt({20.0, 18.0});
  matrix.AddSubJetPt({10.0, 8.0});
  matrix.Initialize();
  if (matrix.Size() != 32)
    return 1;
  
  matrix.AddSubJetPt(40);
  matrix.Initialize();
  if (matrix.Size() != 32)
    return 1;
  
  matrix.AddLeadJetPt(60);
  matrix.Initialize();
  if (matrix.Size() != 56)
    return 1;
  
  matrix.Clear();
  if (matrix.Size() != 0)
    return 1;
  
  matrix.Initialize();
  if (matrix.Size() != 1)
    return 1;
  
  matrix.Clear();
  // build a "simple" test case to do a full check on
  // including background jet definitions, area definitions,
  // etc
  matrix.AddLeadJetPt(18);
  matrix.AddLeadJetR(0.4);
  matrix.AddLeadJetR(0.5);
  matrix.AddSubJetPt(9);
  matrix.AddSubJetR(0.6);
  matrix.AddConstituentEta(1.0);
  matrix.AddConstituentLeadInitialPt(2.5);
  matrix.AddConstituentLeadMatchPt(0.5);
  matrix.AddConstituentSubInitialPt(2.3);
  matrix.AddConstituentSubMatchPt(0.4);
  
  matrix.Initialize();
  
  
  
  // there should be two dijet definitions
  if (matrix.Size() != 2)
    return 1;
  
  for (auto key : matrix.Keys()) {
    double lead_R;
    if ((key.find("INIT_R_0.4") != std::string::npos) &&
        (key.find("MATCH_R_0.4") != std::string::npos)) {
      lead_R = 0.4;
    }
    else if ((key.find("INIT_R_0.5") != std::string::npos) &&
             (key.find("MATCH_R_0.5") != std::string::npos)) {
      lead_R = 0.5;
    }
    else return 1;
    
    // build the background area def
  
    fastjet::GhostedAreaSpec ghost_spec_lead(1.0 + 0.4, 1, 0.01, 1, 0.1, 1e-100);
    fastjet::AreaDefinition bkg_area_def_lead(fastjet::active_area_explicit_ghosts, ghost_spec_lead);
    
    fastjet::GhostedAreaSpec ghost_spec_sub(1.0 + 0.4, 1, 0.01, 1, 0.1, 1e-100);
    fastjet::AreaDefinition bkg_area_def_sub(fastjet::active_area_explicit_ghosts, ghost_spec_sub);
    
    if (!CheckDijetDefinition(*matrix.DijetDefinitions()[key], fastjet::antikt_algorithm, fastjet::antikt_algorithm,
                              lead_R, 0.6, 18, 9, 2.5, 0.5, 2.3, 0.4, 1.0, fastjet::E_scheme, fastjet::Best,
                              fastjet::active_area_explicit_ghosts, 1, 0.01, 1, 0.1, 1e-100,
                              fastjet::JetDefinition(fastjet::kt_algorithm, 0.4), bkg_area_def_lead,
                              bkg_area_def_sub))
      return 1;
  }
  
  // now test that forcing equality of constituent pt & eta in
  // leading & subleading definitions works as intended
  matrix.Clear();
  matrix.ForceConstituentPtEquality(true);
  matrix.ForceConstituentEtaEquality(true);
  matrix.AddConstituentLeadInitialPt({2.0, 3.0});
  matrix.AddConstituentSubInitialPt({2.0, 3.0});
  matrix.AddConstituentEta({1.0, 2.0});
  matrix.Initialize();
  if (matrix.Size() != 4)
    return 1;
  
  matrix.ForceConstituentEtaEquality(false);
  matrix.Initialize();
  if(matrix.Size() != 8)
    return 1;
  
  matrix.Clear();
  matrix.ForceConstituentPtEquality(true);
  matrix.ForceConstituentEtaEquality(true);
  matrix.AddConstituentLeadInitialPt(2.0);
  matrix.AddConstituentSubInitialPt(2.0);
  matrix.AddConstituentEta(1.0);
  matrix.AddConstituentLeadMatchPt({0.2, 0.3});
  matrix.AddConstituentSubMatchPt({0.2, 0.3});
  matrix.Initialize();
  
  if (matrix.Size() != 2)
    return 1;
  
  // test to make sure that a "lonely" pt is not used
  matrix.Clear();
  matrix.ForceConstituentPtEquality(true);
  matrix.ForceConstituentEtaEquality(true);
  matrix.AddConstituentLeadInitialPt({2.0, 3.0});
  matrix.AddConstituentSubInitialPt(2.0);
  matrix.AddConstituentEta(1.0);
  matrix.AddConstituentLeadMatchPt({0.2, 0.3});
  matrix.AddConstituentSubMatchPt({0.2, 0.3});
  matrix.Initialize();
  
  if (matrix.Size() != 2)
    return 1;
  
  return 0;
}
