#include "gtest/gtest.h"

#include "dijetcore/worker/dijet_worker/dijet_worker.h"
#include "dijetcore/util/fastjet/selector_compare.h"
#include "dijetcore/test/dijetcore_test_helper.h"

TEST(DijetWorker, Default) {
  // fast default settings test
  dijetcore::DijetWorker worker(fastjet::antikt_algorithm);
  worker.initialize();
  
  EXPECT_EQ(worker.size(), 1);
    
  
  fastjet::JetDefinition bkg_def(fastjet::kt_algorithm, 0.4);
  fastjet::GhostedAreaSpec ghost_spec(1.4, 1, 0.01, 1, 0.1, 1e-100);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
  
  dijetcore::testing::CheckDijetDefinition(worker.dijetDefinitions().begin()->second.get(), fastjet::antikt_algorithm, fastjet::antikt_algorithm,
                                   0.4, 0.4, 0.4, 0.4, 20.0, 10.0, 2.0, 0.2, 2.0, 0.2, 1.0, fastjet::E_scheme, fastjet::Best,
                                   fastjet::active_area_explicit_ghosts, 1, 0.01, 1, 0.1, 1e-100, bkg_def, area_def, area_def);
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
