// dijet_worker_test.cc

// no tests I write are really "unit tests"
// but these are really not unit tests
// I just use them as a sanity check
// for my core analysis routine

#include "jet_analysis/dijet_worker/dijet_worker.hh"
#include "jet_analysis/util/selector_compare.hh"

#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include <vector>
#include <random>
#include <iostream>

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


int main () {
  
  // fast default settings test
  DijetWorker a(fastjet::antikt_algorithm);
  a.Initialize();
  
  if (a.Size() != 1)
    return 1;
  
  fastjet::JetDefinition bkg_def(fastjet::kt_algorithm, 0.4);
  fastjet::GhostedAreaSpec ghost_spec(1.4, 1, 0.01, 1, 0.1, 1e-100);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
  
  if (!CheckDijetDefinition(*(a.DijetDefinitions().begin()->second), fastjet::antikt_algorithm, fastjet::antikt_algorithm,
                            0.4, 0.4, 20.0, 10.0, 2.0, 0.2, 2.0, 0.2, 1.0, fastjet::E_scheme, fastjet::Best,
                            fastjet::active_area_explicit_ghosts, 1, 0.01, 1, 0.1, 1e-100, bkg_def, area_def, area_def))
    return 1;
  
  
  //do a full clustering routine - first, generate two "hard" jets, with a large background
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
  
  // define distributions for pt, eta, phi
  double eta_max = 1.0;
  double eta_min = -eta_max;
  std::uniform_real_distribution<> eta_dist(eta_min, eta_max);
  
  double phi_max = 2.0 * TMath::Pi();
  double phi_min = 0.0;
  std::uniform_real_distribution<> phi_dist(phi_min, phi_max);
  
  double pt_max = 2.0;
  double pt_min = 0.0;
  std::uniform_real_distribution<> pt_dist(pt_min, pt_max);
  
  // create the RNG
  // standard mersenne_twister_engine seeded with a constant,
  // so that it is reproducible
  std::mt19937 gen(14342);
  
  int n_bkg = 500;
  event_particles.resize(n_bkg + 2);
  // generate n_bkg particle background
  for (int i = 0; i < n_bkg; ++i) {
    double pt_in = pt_dist(gen);
    double eta_in = eta_dist(gen);
    double phi_in = phi_dist(gen);
    double m_in = 0.0;
    fastjet::PseudoJet tmp;
    tmp.reset_PtYPhiM(pt_in, eta_in, phi_in, m_in);
    event_particles[i+2] = tmp;
  }
  
  // define our jetfinding criteria
  fastjet::JetAlgorithm alg = fastjet::antikt_algorithm;
  double lead_jet_pt_cut = 25.0;
  double sublead_jet_pt_cut = 10.0;
  double lead_jet_R = 0.3;
  double sublead_jet_R = 0.25;
  double const_pt_lead_hard = 2.0;
  double const_pt_sub_hard = 2.0;
  double const_pt_lead_match = 0.2;
  double const_pt_sub_match = 0.2;
  
  DijetWorker test_worker(alg, lead_jet_pt_cut, lead_jet_R, sublead_jet_pt_cut, sublead_jet_R,
                          const_pt_lead_hard, const_pt_lead_match, const_pt_sub_hard,
                          const_pt_sub_match);
  
  test_worker.Initialize();
  auto worker_output = test_worker.Run(event_particles);
  if (worker_output.size() != 1)
    return 1;
  
  ClusterOutput dijet_worker_output = (*worker_output.begin()).second;
  
  // now perform the the same calculation without the Dijet Worker
  // to test for equality
  
  // jet definitions
  fastjet::JetDefinition lead_jet_def(alg, lead_jet_R);
  fastjet::JetDefinition sublead_jet_def(alg, sublead_jet_R);
  fastjet::JetDefinition bkg_jet_def(fastjet::kt_algorithm, 0.4);
  // build the selectors we need
  fastjet::Selector constituent_selector_hard = fastjet::SelectorAbsRapMax(1.0) && fastjet::SelectorPtMin(2.0);
  fastjet::Selector constituent_selector_match = fastjet::SelectorAbsRapMax(1.0) && fastjet::SelectorPtMin(0.2);
  fastjet::Selector jet_selector_lead_hard = fastjet::SelectorAbsRapMax(1.0 - lead_jet_R) && fastjet::SelectorPtMin(lead_jet_pt_cut);
  fastjet::Selector jet_selector_sublead_hard = fastjet::SelectorAbsRapMax(1.0 - sublead_jet_R) && fastjet::SelectorPtMin(sublead_jet_pt_cut);
  fastjet::Selector bkg_selector = !fastjet::SelectorNHardest(2) && fastjet::SelectorAbsRapMax(0.6);
  
  // area definitions
  fastjet::GhostedAreaSpec ghost_spec_lead(1.0 + lead_jet_R, 1, 0.01);
  fastjet::GhostedAreaSpec ghost_spec_sub(1.0 + sublead_jet_R, 1, 0.01);
  fastjet::GhostedAreaSpec ghost_spec_bkg(1.0 + 0.4, 1, 0.01);
  fastjet::AreaDefinition area_def_lead(fastjet::active_area_explicit_ghosts, ghost_spec_lead);
  fastjet::AreaDefinition area_def_sub(fastjet::active_area_explicit_ghosts, ghost_spec_sub);
  fastjet::AreaDefinition area_def_bkg(fastjet::active_area_explicit_ghosts, ghost_spec_bkg);
  
  // start by selecting our input particles
  std::vector<fastjet::PseudoJet> cluster_input = constituent_selector_hard(event_particles);
  
  //initial clustering
  fastjet::ClusterSequenceArea cl_lead_hard(cluster_input, lead_jet_def, area_def_lead);
  fastjet::ClusterSequenceArea cl_sublead_hard(cluster_input, sublead_jet_def, area_def_sub);
  
  std::vector<fastjet::PseudoJet> lead_hard_jets = fastjet::sorted_by_pt(jet_selector_lead_hard(cl_lead_hard.inclusive_jets()));
  std::vector<fastjet::PseudoJet> sublead_hard_jets = fastjet::sorted_by_pt(jet_selector_sublead_hard(cl_sublead_hard.inclusive_jets()));
  
  if (lead_hard_jets.size() == 0 || sublead_hard_jets.size() == 0)
    return 1;
  
  // select our leading and subleading hard jets
  fastjet::PseudoJet leading_jet = lead_hard_jets[0];
  
  fastjet::Selector recoil_selector = !fastjet::SelectorRectangle(2.1, TMath::Pi() - dijet_worker_output.dijet_def->dPhi);
  recoil_selector.set_reference(leading_jet);
  std::vector<fastjet::PseudoJet> dphi_cut_jets = fastjet::sorted_by_pt(recoil_selector(sublead_hard_jets));
  
  fastjet::PseudoJet subleading_jet = dphi_cut_jets[0];
  
  // check equality of leading jets
  if (dijet_worker_output.lead_hard.pt() != leading_jet.pt() ||
      dijet_worker_output.lead_hard.eta() != leading_jet.eta() ||
      dijet_worker_output.lead_hard.phi() != leading_jet.phi() ||
      dijet_worker_output.lead_hard.m() != leading_jet.m())
    return 1;
  
  if (dijet_worker_output.sublead_hard.pt() != subleading_jet.pt() ||
      dijet_worker_output.sublead_hard.eta() != subleading_jet.eta() ||
      dijet_worker_output.sublead_hard.phi() != subleading_jet.phi() ||
      dijet_worker_output.sublead_hard.m() != subleading_jet.m())
    return 1;
  
  if (fabs(dijet_worker_output.lead_hard.area() - leading_jet.area()) > 0.001 ||
      fabs(dijet_worker_output.sublead_hard.area() != subleading_jet.area()) > 0.001 )
    return 1;
  
  if (fabs(dijet_worker_output.lead_hard_rho) > 1e-10 ||
      fabs(dijet_worker_output.lead_hard_sigma) > 1e-10 ||
      fabs(dijet_worker_output.sublead_hard_rho) > 1e-10 ||
      fabs(dijet_worker_output.sublead_hard_sigma) > 1e-10)
    return 1;
  
  // now get the matched set of input particles
  cluster_input = constituent_selector_match(event_particles);
  
  fastjet::ClusterSequenceArea cl_lead_match(cluster_input, lead_jet_def, area_def_lead);
  fastjet::ClusterSequenceArea cl_sublead_match(cluster_input, sublead_jet_def, area_def_sub);
  
  std::vector<fastjet::PseudoJet> lead_match_jets = fastjet::sorted_by_pt(cl_lead_match.inclusive_jets());
  std::vector<fastjet::PseudoJet> sublead_match_jets = fastjet::sorted_by_pt(cl_sublead_match.inclusive_jets());
  
  if (lead_match_jets.size() == 0 || sublead_match_jets.size() == 0)
    return 1;
  
  // do background subtraction
  fastjet::JetMedianBackgroundEstimator lead_bkg_est(bkg_selector, bkg_def, area_def_bkg);
  lead_bkg_est.set_particles(cluster_input);
  fastjet::JetMedianBackgroundEstimator sublead_bkg_est(bkg_selector, bkg_def, area_def_bkg);
  sublead_bkg_est.set_particles(cluster_input);
  
  fastjet::Subtractor lead_subtractor(&lead_bkg_est);
  fastjet::Subtractor sublead_subtractor(&sublead_bkg_est);
  
  std::vector<fastjet::PseudoJet> lead_subtracted_jets = fastjet::sorted_by_pt(lead_subtractor(lead_match_jets));
  std::vector<fastjet::PseudoJet> sublead_subtracted_jets = fastjet::sorted_by_pt(sublead_subtractor(sublead_match_jets));
  
  // match to the hard jets
  fastjet::Selector lead_circle = fastjet::SelectorCircle(lead_jet_R);
  lead_circle.set_reference(leading_jet);
  fastjet::Selector sublead_circle = fastjet::SelectorCircle(lead_jet_R);
  sublead_circle.set_reference(subleading_jet);
  
  std::vector<fastjet::PseudoJet> matched_to_lead = fastjet::sorted_by_pt(lead_circle(lead_subtracted_jets));
  std::vector<fastjet::PseudoJet> matched_to_sublead = fastjet::sorted_by_pt(sublead_circle(sublead_subtracted_jets));
  
  if (matched_to_lead.size() == 0 || matched_to_sublead.size() == 0)
    return 1;
  
  
  fastjet::PseudoJet matched_lead_jet =  matched_to_lead[0];
  fastjet::PseudoJet matched_sublead_jet =  matched_to_sublead[0];
  
  // now check the matched jets agree with the dijetworker results
  if (fabs(matched_lead_jet.eta() - dijet_worker_output.lead_match.eta()) > 0.05 ||
      fabs(matched_sublead_jet.eta() - dijet_worker_output.sublead_match.eta()) > 0.05 ||
      fabs(matched_lead_jet.phi() - dijet_worker_output.lead_match.phi()) > 0.05 ||
      fabs(matched_sublead_jet.phi() - dijet_worker_output.sublead_match.phi()) > 0.05)
    return 1;
  
  if (fabs((dijet_worker_output.lead_match.pt() - matched_lead_jet.pt())/dijet_worker_output.lead_match.pt()) > 0.1 ||
      fabs((dijet_worker_output.sublead_match.pt() - matched_sublead_jet.pt())/dijet_worker_output.sublead_match.pt()) > 0.1)
    return 1;
  
  // **************************
  // next test: see if the jets bkg is being calculated properly, by finding the average
  // over an ensemble of events
  double lead_pt_total = 0.0;
  double sublead_pt_total = 0.0;
  
  int n_tries = 200;
  for (int i = 0; i < n_tries; ++i) {
    // set the inital jets
    event_particles[0] = leading_jet_in;
    event_particles[1] = subleading_jet_in;
    // generate n_bkg particle background
    for (int i = 0; i < n_bkg; ++i) {
      double pt_in = pt_dist(gen);
      double eta_in = eta_dist(gen);
      double phi_in = phi_dist(gen);
      double m_in = 0.0;
      fastjet::PseudoJet tmp;
      tmp.reset_PtYPhiM(pt_in, eta_in, phi_in, m_in);
      event_particles[i+2] = tmp;
    }
    
    auto result = test_worker.Run(event_particles);
    ClusterOutput out = (*result.begin()).second;
    lead_pt_total += out.lead_match.pt();
    sublead_pt_total += out.sublead_match.pt();
  }
  
  if (fabs((lead_pt_total/n_tries-lead_pt)/lead_pt) > 0.05 ||
      fabs((sublead_pt_total/n_tries-sub_pt)/sub_pt) > 0.05)
    return 1;
  
  // **************************
  // next test: make sure the clustering output is successfully finding events
  // that have the leading jet, the subleading jet, and the matched jets
  cluster_input.clear();
  leading_jet_in.reset_PtYPhiM(21, 0.1, 3.14, 0);
  subleading_jet_in.reset_PtYPhiM(9, -0.1, 0.01, 0);
  cluster_input.push_back(leading_jet_in);
  cluster_input.push_back(subleading_jet_in);
  
  DijetWorker b(fastjet::antikt_algorithm);
  b.Initialize();
  worker_output = b.Run(cluster_input);
  
  if (worker_output.size() != 1)
    return 1;
  
  auto only_output = worker_output.begin();
  if ((*only_output).second.found_lead != true ||
      (*only_output).second.found_sublead != false)
    return 1;
  
  // now make sure it can find a subleading jet & match
  cluster_input.clear();
  subleading_jet_in.reset_PtYPhiM(11, -0.1, 0.01, 0);
  cluster_input.push_back(leading_jet_in);
  cluster_input.push_back(subleading_jet_in);
  
  worker_output = b.Run(cluster_input);
  if (worker_output.size() != 1)
    return 1;
  
  only_output = worker_output.begin();
  if ((*only_output).second.found_lead != true ||
      (*only_output).second.found_sublead != true ||
      (*only_output).second.found_match != true)
    return 1;
  
  return 0;
}
