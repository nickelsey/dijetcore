#include "gtest/gtest.h"

#include "dijetcore/worker/dijet_worker/dijet_worker.h"
#include "dijetcore/util/fastjet/selector_compare.h"

bool CheckDijetDefinition(dijetcore::DijetDefinition def, fastjet::JetAlgorithm lead_alg, fastjet::JetAlgorithm sub_alg,
                          double lead_R, double sub_R, double lead_pt, double sub_pt, double lead_const_init_pt,
                          double lead_const_match_pt, double sub_const_init_pt, double sub_const_match_pt,
                          double const_eta, fastjet::RecombinationScheme scheme, fastjet::Strategy strategy,
                          fastjet::AreaType area_type, int ghost_repeat, double ghost_area, double grid_scatter,
                          double pt_scatter, double mean_ghost_pt, fastjet::JetDefinition bkg_def,
                          fastjet::AreaDefinition bkg_area_lead, fastjet::AreaDefinition bkg_area_sub) {
  dijetcore::MatchDef* lead = def.lead;
  dijetcore::MatchDef* sub = def.sub;
  
  // check the leading jet
  if (lead->InitialJetDef().R() != lead_R ||
      lead->MatchedJetDef().R() != lead_R ||
      lead->InitialJetDef().jet_algorithm() != lead_alg ||
      lead->MatchedJetDef().jet_algorithm() != lead_alg ||
      dijetcore::ExtractDoubleFromSelector(lead->InitialJetDef().ConstituentSelector(), "pt >=") != lead_const_init_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->MatchedJetDef().ConstituentSelector(), "pt >=") != lead_const_match_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->InitialJetDef().JetSelector(), "pt >=") != lead_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->MatchedJetDef().JetSelector(), "pt >=") != 0 ||
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
      dijetcore::ExtractDoubleFromSelector(sub->InitialJetDef().ConstituentSelector(), "pt >=") != sub_const_init_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->MatchedJetDef().ConstituentSelector(), "pt >=") != sub_const_match_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->InitialJetDef().JetSelector(), "pt >=") != sub_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->MatchedJetDef().JetSelector(), "pt >=") != 0 ||
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

TEST(DijetWorker, Default) {
  // fast default settings test
  dijetcore::DijetWorker worker(fastjet::antikt_algorithm);
  worker.Initialize();
  
  EXPECT_EQ(worker.Size(), 1);
    
  
  fastjet::JetDefinition bkg_def(fastjet::kt_algorithm, 0.4);
  fastjet::GhostedAreaSpec ghost_spec(1.4, 1, 0.01, 1, 0.1, 1e-100);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
  
  EXPECT_TRUE(CheckDijetDefinition(*(worker.DijetDefinitions().begin()->second), fastjet::antikt_algorithm, fastjet::antikt_algorithm,
                                   0.4, 0.4, 20.0, 10.0, 2.0, 0.2, 2.0, 0.2, 1.0, fastjet::E_scheme, fastjet::Best,
                                   fastjet::active_area_explicit_ghosts, 1, 0.01, 1, 0.1, 1e-100, bkg_def, area_def, area_def));
}

TEST(DijetWorker, SimpleCluster) {
  // fast default settings test
  dijetcore::DijetWorker worker(fastjet::antikt_algorithm);
  worker.Initialize();
  
  EXPECT_EQ(worker.Size(), 1);
  
  //do a full clustering routine - generate two "hard" jets
  double lead_pt = 30.0;
  double lead_eta = 0.4;
  double lead_phi = 0.4;
  double lead_m = 0;
  double sub_pt = 15.0;
  double sub_eta = -0.4;
  double sub_phi = 3.54;
  double sub_m = 0;
  
  fastjet::PseudoJet leading_jet_in;
  leading_jet_in.reset_PtYPhiM(lead_pt, lead_eta, lead_phi, lead_m);
  fastjet::PseudoJet subleading_jet_in;
  subleading_jet_in.reset_PtYPhiM(sub_pt, sub_eta, sub_phi, sub_m);
  
  std::vector<fastjet::PseudoJet> event_particles{leading_jet_in, subleading_jet_in};
  
  auto& worker_output = worker.Run(event_particles);
  
  ASSERT_EQ(worker_output.size(), 1);
  
  dijetcore::ClusterOutput dijet_worker_output = *(*worker_output.begin()).second.get();
  
}

TEST(DijetWorker, LessSimpleCluster) {
  // fast default settings test
  dijetcore::DijetWorker worker(fastjet::antikt_algorithm);
  
  worker.AddSubJetPt({10.0, 16.0});
  
  worker.Initialize();
  
  EXPECT_EQ(worker.Size(), 2);
  
  //do a full clustering routine - generate two "hard" jets
  double lead_pt = 30.0;
  double lead_eta = 0.4;
  double lead_phi = 0.4;
  double lead_m = 0;
  double sub_pt = 15.0;
  double sub_eta = -0.4;
  double sub_phi = 3.54;
  double sub_m = 0;
  
  fastjet::PseudoJet leading_jet_in;
  leading_jet_in.reset_PtYPhiM(lead_pt, lead_eta, lead_phi, lead_m);
  fastjet::PseudoJet subleading_jet_in;
  subleading_jet_in.reset_PtYPhiM(sub_pt, sub_eta, sub_phi, sub_m);
  
  std::vector<fastjet::PseudoJet> event_particles{leading_jet_in, subleading_jet_in};
  
  auto& worker_output = worker.Run(event_particles);
  
  ASSERT_EQ(worker_output.size(), 2);
  
  int found_sub = 0;
  int found_lead = 0;
  int found_match = 0;
  for (auto& out : worker_output) {
    if (out.second->found_lead)
      found_lead++;
    if (out.second->found_sublead)
      found_sub++;
    if (out.second->found_match)
      found_match++;
  }
  
  EXPECT_EQ(found_sub, 1);
  EXPECT_EQ(found_lead, 2);
  EXPECT_EQ(found_match, 1);
  
  for (auto& out : worker_output) {
    if (out.second->found_lead)
      EXPECT_NEAR(out.second->lead_hard.pt(), lead_pt, 1e-1);
    if (out.second->found_sublead)
      EXPECT_NEAR(out.second->sublead_hard.pt(), sub_pt, 1e-1);
    if (out.second->found_match) {
      EXPECT_NEAR(out.second->lead_match.pt(), lead_pt, 1e-1);
      EXPECT_NEAR(out.second->sublead_match.pt(), sub_pt, 1e-1);
    }
  }
}
