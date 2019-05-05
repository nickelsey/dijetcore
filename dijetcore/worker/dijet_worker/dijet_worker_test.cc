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
  if (lead->initialJetDef().R() != lead_R ||
      lead->matchedJetDef().R() != lead_R ||
      lead->initialJetDef().jet_algorithm() != lead_alg ||
      lead->matchedJetDef().jet_algorithm() != lead_alg ||
      dijetcore::ExtractDoubleFromSelector(lead->initialJetDef().constituentSelector(), "pt >=") != lead_const_init_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->matchedJetDef().constituentSelector(), "pt >=") != lead_const_match_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->initialJetDef().jetSelector(), "pt >=") != lead_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->matchedJetDef().jetSelector(), "pt >=") != 0 ||
      lead->initialJetDef().recombination_scheme() != scheme ||
      lead->matchedJetDef().recombination_scheme() != scheme ||
      lead->initialJetDef().strategy() != strategy ||
      lead->matchedJetDef().strategy() != strategy ||
      lead->initialJetDef().areaDefinition().area_type() != area_type ||
      lead->matchedJetDef().areaDefinition().area_type() != area_type ||
      lead->initialJetDef().areaDefinition().ghost_spec().ghost_maxrap() != (const_eta + lead_R) ||
      lead->matchedJetDef().areaDefinition().ghost_spec().ghost_maxrap() != (const_eta + lead_R) ||
      lead->initialJetDef().areaDefinition().ghost_spec().repeat() != ghost_repeat ||
      lead->matchedJetDef().areaDefinition().ghost_spec().repeat() != ghost_repeat ||
      lead->initialJetDef().areaDefinition().ghost_spec().ghost_area() != ghost_area ||
      lead->matchedJetDef().areaDefinition().ghost_spec().ghost_area() != ghost_area ||
      lead->initialJetDef().areaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      lead->matchedJetDef().areaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      lead->initialJetDef().areaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      lead->matchedJetDef().areaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      lead->initialJetDef().areaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      lead->matchedJetDef().areaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      lead->initialJetDef().backgroundJetDef().R() != bkg_def.R() ||
      lead->matchedJetDef().backgroundJetDef().R() != bkg_def.R() ||
      lead->initialJetDef().backgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      lead->matchedJetDef().backgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      lead->initialJetDef().backgroundJetDef().strategy() != bkg_def.strategy() ||
      lead->matchedJetDef().backgroundJetDef().strategy() != bkg_def.strategy() ||
      lead->initialJetDef().backgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      lead->matchedJetDef().backgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      lead->initialJetDef().backgroundAreaDef().area_type() != bkg_area_lead.area_type() ||
      lead->matchedJetDef().backgroundAreaDef().area_type() != bkg_area_lead.area_type() ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_lead.ghost_spec().ghost_maxrap() ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_lead.ghost_spec().ghost_maxrap() ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().repeat() != bkg_area_lead.ghost_spec().repeat() ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().repeat() != bkg_area_lead.ghost_spec().repeat() ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().ghost_area() != bkg_area_lead.ghost_spec().ghost_area() ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().ghost_area() != bkg_area_lead.ghost_spec().ghost_area() ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_lead.ghost_spec().grid_scatter()  ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_lead.ghost_spec().grid_scatter()  ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_lead.ghost_spec().pt_scatter()  ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_lead.ghost_spec().pt_scatter()  ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_lead.ghost_spec().mean_ghost_pt()  ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_lead.ghost_spec().mean_ghost_pt() )
    return false;
  
  // and check subleading jet
  if (sub->initialJetDef().R() != sub_R ||
      sub->matchedJetDef().R() != sub_R ||
      sub->initialJetDef().jet_algorithm() != sub_alg ||
      sub->matchedJetDef().jet_algorithm() != sub_alg ||
      dijetcore::ExtractDoubleFromSelector(sub->initialJetDef().constituentSelector(), "pt >=") != sub_const_init_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->matchedJetDef().constituentSelector(), "pt >=") != sub_const_match_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->initialJetDef().jetSelector(), "pt >=") != sub_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->matchedJetDef().jetSelector(), "pt >=") != 0 ||
      sub->initialJetDef().recombination_scheme() != scheme ||
      sub->matchedJetDef().recombination_scheme() != scheme ||
      sub->initialJetDef().strategy() != strategy ||
      sub->matchedJetDef().strategy() != strategy ||
      sub->initialJetDef().areaDefinition().area_type() != area_type ||
      sub->matchedJetDef().areaDefinition().area_type() != area_type ||
      sub->initialJetDef().areaDefinition().ghost_spec().ghost_maxrap() != (const_eta + sub_R) ||
      sub->matchedJetDef().areaDefinition().ghost_spec().ghost_maxrap() != (const_eta + sub_R) ||
      sub->initialJetDef().areaDefinition().ghost_spec().repeat() != ghost_repeat ||
      sub->matchedJetDef().areaDefinition().ghost_spec().repeat() != ghost_repeat ||
      sub->initialJetDef().areaDefinition().ghost_spec().ghost_area() != ghost_area ||
      sub->matchedJetDef().areaDefinition().ghost_spec().ghost_area() != ghost_area ||
      sub->initialJetDef().areaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      sub->matchedJetDef().areaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      sub->initialJetDef().areaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      sub->matchedJetDef().areaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      sub->initialJetDef().areaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      sub->matchedJetDef().areaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      sub->initialJetDef().backgroundJetDef().R() != bkg_def.R() ||
      sub->matchedJetDef().backgroundJetDef().R() != bkg_def.R() ||
      sub->initialJetDef().backgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      sub->matchedJetDef().backgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      sub->initialJetDef().backgroundJetDef().strategy() != bkg_def.strategy() ||
      sub->matchedJetDef().backgroundJetDef().strategy() != bkg_def.strategy() ||
      sub->initialJetDef().backgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      sub->matchedJetDef().backgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      sub->initialJetDef().backgroundAreaDef().area_type() != bkg_area_sub.area_type() ||
      sub->matchedJetDef().backgroundAreaDef().area_type() != bkg_area_sub.area_type() ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_sub.ghost_spec().ghost_maxrap() ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_sub.ghost_spec().ghost_maxrap() ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().repeat() != bkg_area_sub.ghost_spec().repeat() ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().repeat() != bkg_area_sub.ghost_spec().repeat() ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().ghost_area() != bkg_area_sub.ghost_spec().ghost_area() ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().ghost_area() != bkg_area_sub.ghost_spec().ghost_area() ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_sub.ghost_spec().grid_scatter()  ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_sub.ghost_spec().grid_scatter()  ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_sub.ghost_spec().pt_scatter()  ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_sub.ghost_spec().pt_scatter()  ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_sub.ghost_spec().mean_ghost_pt()  ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_sub.ghost_spec().mean_ghost_pt() )
    return false;
  
  return true;
}

TEST(DijetWorker, Default) {
  // fast default settings test
  dijetcore::DijetWorker worker(fastjet::antikt_algorithm);
  worker.initialize();
  
  EXPECT_EQ(worker.size(), 1);
    
  
  fastjet::JetDefinition bkg_def(fastjet::kt_algorithm, 0.4);
  fastjet::GhostedAreaSpec ghost_spec(1.4, 1, 0.01, 1, 0.1, 1e-100);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
  
  EXPECT_TRUE(CheckDijetDefinition(*(worker.dijetDefinitions().begin()->second), fastjet::antikt_algorithm, fastjet::antikt_algorithm,
                                   0.4, 0.4, 20.0, 10.0, 2.0, 0.2, 2.0, 0.2, 1.0, fastjet::E_scheme, fastjet::Best,
                                   fastjet::active_area_explicit_ghosts, 1, 0.01, 1, 0.1, 1e-100, bkg_def, area_def, area_def));
}

TEST(DijetWorker, SimpleCluster) {
  // fast default settings test
  dijetcore::DijetWorker worker(fastjet::antikt_algorithm);
  worker.initialize();
  
  EXPECT_EQ(worker.size(), 1);
  
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
  
  auto& worker_output = worker.run(event_particles);
  
  ASSERT_EQ(worker_output.size(), 1);
  
  dijetcore::ClusterOutput dijet_worker_output = *(*worker_output.begin()).second.get();
  
}

TEST(DijetWorker, LessSimpleCluster) {
  // fast default settings test
  dijetcore::DijetWorker worker(fastjet::antikt_algorithm);
  
  worker.addSubJetPt({10.0, 16.0});
  
  worker.initialize();
  
  EXPECT_EQ(worker.size(), 2);
  
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
  
  auto& worker_output = worker.run(event_particles);
  
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
