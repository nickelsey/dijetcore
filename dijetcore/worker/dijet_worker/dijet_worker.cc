#include "dijetcore/worker/dijet_worker/dijet_worker.h"
#include "dijetcore/lib/logging.h"

#include "TMath.h"

#include <iostream>

namespace dijetcore {
  
  // default constructor
  DijetWorker::DijetWorker() : DijetMatrix() { }
  
  // single entry construction
  DijetWorker::DijetWorker(fastjet::JetAlgorithm jet_alg_in,
                           double lead_pt_in,
                           double lead_R_in,
                           double lead_match_R_in,
                           double sub_pt_in,
                           double sub_R_in,
                           double sub_match_R_in,
                           double const_lead_pt_init_in,
                           double const_lead_pt_match_in,
                           double const_sub_pt_init_in,
                           double const_sub_pt_match_in,
                           double eta_in) :
  DijetMatrix(jet_alg_in, lead_pt_in, lead_R_in, lead_match_R_in, sub_pt_in, sub_R_in,
              sub_match_R_in, const_lead_pt_init_in, const_lead_pt_match_in,
              const_sub_pt_init_in, const_sub_pt_match_in, eta_in) { };
  
  // set construction
  DijetWorker::DijetWorker(std::set<fastjet::JetAlgorithm> jet_alg_in,
                           std::set<double> lead_pt_in,
                           std::set<double> lead_R_in,
                           std::set<double> lead_match_R_in,
                           std::set<double> sub_pt_in,
                           std::set<double> sub_R_in,
                           std::set<double> sub_match_R_in,
                           std::set<double> const_lead_pt_init_in,
                           std::set<double> const_lead_pt_match_in,
                           std::set<double> const_sub_pt_init_in,
                           std::set<double> const_sub_pt_match_in,
                           std::set<double> eta_in) :
  DijetMatrix(jet_alg_in, lead_pt_in, lead_R_in, lead_match_R_in, sub_pt_in, sub_R_in,
              sub_match_R_in, const_lead_pt_init_in, const_lead_pt_match_in,
              const_sub_pt_init_in, const_sub_pt_match_in, eta_in) { };
  
  DijetWorker::DijetWorker(const DijetMatrix& rhs) :
  DijetMatrix(rhs) { }
  
  std::unordered_map<std::string, unique_ptr<ClusterOutput>>&
  DijetWorker::run(const std::vector<fastjet::PseudoJet>& input) {
    // clear the last event
    clearResults();

    // make sure the DijetMatrix is initialized,
    // if not, do so
    if (size() == 0) {
      initialize();
    }
    
    // get the sorted list of all dijet definitions
    auto& sorted_dijet_definitions = sortedDefinitions();
    auto& sorted_dijet_keys = sortedKeys();
    auto& sorted_dijet_min_pt = sortedDefinitionsMinPt();
    
    for (int i = 0; i < sorted_dijet_definitions.size(); ++i) {

      auto& dijet_set = sorted_dijet_definitions.at(i);
      const string& dijet_set_key = sorted_dijet_keys.at(i);
      
      double min_pt_lead = sorted_dijet_min_pt[i].first;
      double min_pt_sub  = sorted_dijet_min_pt[i].second;
      
      if (dijet_set.size() == 0) {
        LOG(ERROR) << "error in DijetMatrix: size of sorted DijetDefinition container is zero";
        continue;
      }
      
      // now lets loop over all dijet definitions in the set
      for (auto&& dijet_def : dijet_set) {

        // create output container
        unique_ptr<ClusterOutput> dijet_container = make_unique<ClusterOutput>(dijet_def.second);
        
        // load the dijet key & Matchdefs
        const string& key = dijet_def.first;
        MatchDef* lead = dijet_def.second->lead;
        MatchDef* sub = dijet_def.second->sub;
        
        // run the clustering if its necessary
        if (!runClustering(input, *dijet_def.second, dijet_set_key, min_pt_lead, min_pt_sub)) {
          break;
        }

        // load proper cluster sequence
        auto cl_sequences = loadProperClusterSequence(dijet_set_key);
        bool fail = false;
        for (auto& e : cl_sequences)
          if (e == nullptr)
            fail = true;

        if (fail) {
          LOG(ERROR) << "failed to load cluster sequences: " << key;
          continue;
        }

        fastjet::ClusterSequenceArea* cl_lead_hard = cl_sequences[0];
        fastjet::ClusterSequenceArea* cl_lead_match = cl_sequences[2];
        fastjet::ClusterSequenceArea* cl_sub_hard = cl_sequences[1];
        fastjet::ClusterSequenceArea* cl_sub_match = cl_sequences[3];

        // now we can decide if we have hard core jets that satisfy our dijet criteria
        // get output jets, and make sure neither are zero length
        std::vector<fastjet::PseudoJet> lead_hard_jets = fastjet::sorted_by_pt(lead->initialJetDef().jetSelector()(cl_lead_hard->inclusive_jets()));
        std::vector<fastjet::PseudoJet> sublead_hard_jets = fastjet::sorted_by_pt(sub->initialJetDef().jetSelector()(cl_sub_hard->inclusive_jets()));
        
        if (lead_hard_jets.size() == 0)
          continue;
        
        // now do background subtraction for these jets
        std::vector<fastjet::PseudoJet> lead_subtracted_hard;
        std::vector<fastjet::PseudoJet> sub_subtracted_hard;
        
        if (equivalentBkgEstimationInput(lead->initialJetDef(), sub->initialJetDef())) {
          auto bkg_subtractor = getBackgroundSubtractor(input, lead->initialJetDef());
          lead_subtracted_hard = lead->initialJetDef().jetSelector()((*bkg_subtractor.second)(lead_hard_jets));
          sub_subtracted_hard = sub->initialJetDef().jetSelector()((*bkg_subtractor.second)(sublead_hard_jets));
          dijet_container->lead_hard_rho = bkg_subtractor.first->rho();
          dijet_container->sublead_hard_rho = bkg_subtractor.first->rho();
          dijet_container->lead_hard_sigma = bkg_subtractor.first->sigma();
          dijet_container->sublead_hard_sigma = bkg_subtractor.first->sigma();
        }
        else {
          auto bkg_subtractor_lead = getBackgroundSubtractor(input, lead->initialJetDef());
          auto bkg_subtractor_sub = getBackgroundSubtractor(input, sub->initialJetDef());
          lead_subtracted_hard = lead->initialJetDef().jetSelector()((*bkg_subtractor_lead.second)(lead_hard_jets));
          sublead_hard_jets = sub->initialJetDef().jetSelector()((*bkg_subtractor_sub.second)(sublead_hard_jets));
          dijet_container->lead_hard_rho = bkg_subtractor_lead.first->rho();
          dijet_container->sublead_hard_rho = bkg_subtractor_sub.first->rho();
          dijet_container->lead_hard_sigma = bkg_subtractor_lead.first->sigma();
          dijet_container->sublead_hard_sigma = bkg_subtractor_sub.first->sigma();
        }
        
        // check to make sure we still have a leading jet
        if (lead_subtracted_hard.size() == 0)
          continue;
        
        // select our trigger jet & recoil jet
        fastjet::PseudoJet leading_hard_jet = selectLeadingHardJet(lead_subtracted_hard);
        fastjet::PseudoJet subleading_hard_jet = selectSubLeadingHardJet(leading_hard_jet,
                                                                         sub_subtracted_hard,
                                                                         dijet_def.second->dPhi);

        // we have a trigger jet at least
        dijet_container->found_lead = true;
        dijet_container->lead_hard = leading_hard_jet;

        // make sure the subleading hard jet isn't a failed match
        if (!subleading_hard_jet.has_associated_cluster_sequence()) {
          cluster_result[key] = std::move(dijet_container);
          continue;
        }

        // fill in the subleading jet
        dijet_container->found_sublead = true;
        dijet_container->sublead_hard = subleading_hard_jet;

        // Now for matched jets
        // get the resulting jets and make sure neither are zero length
        std::vector<fastjet::PseudoJet> lead_match_jets = fastjet::sorted_by_pt(cl_lead_match->inclusive_jets());
        std::vector<fastjet::PseudoJet> sublead_match_jets = fastjet::sorted_by_pt(cl_sub_match->inclusive_jets());

        // if there are no jets (should not happen) then continue
        if (lead_match_jets.size() == 0 || sublead_match_jets.size() == 0) {
          cluster_result[key] = std::move(dijet_container);
          continue;
        }

        // do background estimation & subtraction for the matched jets
        std::vector<fastjet::PseudoJet> lead_subtracted_matched;
        std::vector<fastjet::PseudoJet> sub_subtracted_matched;

        std::pair<double, double> lead_match_bkg = subtractBackgroundFromJets(input, lead->matchedJetDef(),
                                                                              lead_match_jets, lead_subtracted_matched);
        std::pair<double, double> sub_match_bkg = subtractBackgroundFromJets(input, sub->matchedJetDef(),
                                                                             sublead_match_jets, sub_subtracted_matched);

        // get the resulting jets and make sure neither are zero length
        lead_subtracted_matched = fastjet::sorted_by_pt(lead->matchedJetDef().jetSelector()(lead_subtracted_matched));
        sub_subtracted_matched = fastjet::sorted_by_pt(sub->matchedJetDef().jetSelector()(sub_subtracted_matched));

        if (lead_subtracted_matched.size() == 0 ||
            sub_subtracted_matched.size() == 0) {
          cluster_result[key] = std::move(dijet_container);
          continue;
        }

        // now we perform matching
        fastjet::PseudoJet leading_matched_jet = matchJets(leading_hard_jet, *lead, lead_subtracted_matched);
        fastjet::PseudoJet subleading_matched_jet = matchJets(subleading_hard_jet, *sub, sub_subtracted_matched);

        if (!leading_matched_jet.has_associated_cluster_sequence() ||
            !subleading_matched_jet.has_associated_cluster_sequence()) {
          cluster_result[key] = std::move(dijet_container);
          continue;
        }

        // success - fill in the rest of the ClusterOutput
        dijet_container->found_match = true;
        dijet_container->lead_match = leading_matched_jet;
        dijet_container->lead_match_rho = lead_match_bkg.first;
        dijet_container->lead_match_sigma = lead_match_bkg.second;
        dijet_container->sublead_match = subleading_matched_jet;
        dijet_container->sublead_match_rho = sub_match_bkg.first;
        dijet_container->sublead_match_sigma = sub_match_bkg.second;
        cluster_result[key] = std::move(dijet_container);
      }
    }
    return cluster_result;
  }
  
  void DijetWorker::clearResults() {
    cluster_seq_lead_hard.clear();
    cluster_seq_sub_hard.clear();
    cluster_seq_lead_match.clear();
    cluster_seq_sub_match.clear();
    cluster_result.clear();
  }
  
  bool DijetWorker::runClustering(const std::vector<fastjet::PseudoJet>& input, const DijetDefinition& def,
                                  string cluster_identifier, double min_pt_lead, double min_pt_sub) {
    // check if the clustering has been run yet -
    // if not, run the clustering for this dijet set
    if (cluster_seq_lead_hard.find(cluster_identifier) == cluster_seq_lead_hard.end()) {
      MatchDef* lead = def.lead;
      MatchDef* sub = def.sub;
      
      unique_ptr<fastjet::ClusterSequenceArea> lead_hard_seq(nullptr);
      unique_ptr<fastjet::ClusterSequenceArea> lead_match_seq(nullptr);
      unique_ptr<fastjet::ClusterSequenceArea> sub_hard_seq(nullptr);
      unique_ptr<fastjet::ClusterSequenceArea> sub_match_seq(nullptr);
      
      // first, the hard core jets
      
      // if the same cluster sequence can be used for both leading and subleading, do so
      if (equivalentClusterInput(lead->initialJetDef(), sub->initialJetDef())) {
        lead_hard_seq = make_unique<fastjet::ClusterSequenceArea>(lead->initialJetDef().constituentSelector()(input),
                                                                  lead->initialJetDef(),
                                                                  lead->initialJetDef().areaDefinition());
      }
      else {
        lead_hard_seq = make_unique<fastjet::ClusterSequenceArea>(lead->initialJetDef().constituentSelector()(input),
                                                                  lead->initialJetDef(),
                                                                  lead->initialJetDef().areaDefinition());
        sub_hard_seq = make_unique<fastjet::ClusterSequenceArea>(sub->initialJetDef().constituentSelector()(input),
                                                                 sub->initialJetDef(),
                                                                 sub->initialJetDef().areaDefinition());
      }
      
      // insert into dictionary
      if (sub_hard_seq.get() == nullptr) {
        cluster_seq_lead_hard[cluster_identifier] = std::move(lead_hard_seq);
      }
      else {
        cluster_seq_lead_hard[cluster_identifier] = std::move(lead_hard_seq);
        cluster_seq_sub_hard[cluster_identifier] = std::move(sub_hard_seq);
      }
      
      // if there are no jets we can continue -
      if (cluster_seq_sub_hard.find(cluster_identifier) != cluster_seq_sub_hard.end() &&
          cluster_seq_sub_hard[cluster_identifier].get() != nullptr) {
        if (cluster_seq_sub_hard[cluster_identifier]->inclusive_jets().size() < 1 ||
            cluster_seq_lead_hard[cluster_identifier]->inclusive_jets().size() < 1)
          return false;
      }
      else {
        if (cluster_seq_lead_hard[cluster_identifier]->inclusive_jets().size() < 2)
          return false;
      }
      
      // if the jets do not satisfy our min pT cuts, we can continue - we will
      // identify potential leading & subleading jets, and check against min pT
      fastjet::PseudoJet tmp_lead = fastjet::sorted_by_pt(cluster_seq_lead_hard[cluster_identifier]->inclusive_jets())[0];
      fastjet::Selector tmp_sel = !fastjet::SelectorCircle(lead->initialJetDef().R());
      tmp_sel.set_reference(tmp_lead);
      
      fastjet::PseudoJet tmp_sub;
      if (cluster_seq_sub_hard.find(cluster_identifier) != cluster_seq_sub_hard.end()) {
        if (tmp_sel(cluster_seq_sub_hard[cluster_identifier]->inclusive_jets()).size() == 0)
          return false;
        tmp_sub = fastjet::sorted_by_pt(tmp_sel(cluster_seq_sub_hard[cluster_identifier]->inclusive_jets()))[0];
      }
      
      else {
        if (tmp_sel(cluster_seq_lead_hard[cluster_identifier]->inclusive_jets()).size() == 0)
          return false;
        tmp_sub = fastjet::sorted_by_pt(tmp_sel(cluster_seq_lead_hard[cluster_identifier]->inclusive_jets()))[0];
      }
      
      if (tmp_lead.pt() < min_pt_lead || tmp_sub.pt() < min_pt_sub)
        return false;
      
      // now we will perform the matching jet finding
      // for both leading & subleading jets
      
      if (equivalentClusterInput(lead->matchedJetDef(), sub->matchedJetDef())) {
        lead_match_seq = make_unique<fastjet::ClusterSequenceArea>(lead->matchedJetDef().constituentSelector()(input),
                                                                   lead->matchedJetDef(),
                                                                   lead->matchedJetDef().areaDefinition());
      }
      else {
        lead_match_seq = make_unique<fastjet::ClusterSequenceArea>(lead->matchedJetDef().constituentSelector()(input),
                                                                   lead->matchedJetDef(),
                                                                   lead->matchedJetDef().areaDefinition());
        sub_match_seq = make_unique<fastjet::ClusterSequenceArea>(sub->matchedJetDef().constituentSelector()(input),
                                                                  sub->matchedJetDef(),
                                                                  sub->matchedJetDef().areaDefinition());
      }
      
      // insert into dictionary
      if (sub_match_seq.get() == nullptr) {
        cluster_seq_lead_match[cluster_identifier] = std::move(lead_match_seq);
      }
      else {
        cluster_seq_lead_match[cluster_identifier] = std::move(lead_match_seq);
        cluster_seq_sub_match[cluster_identifier] = std::move(sub_match_seq);
      }
    }
    return true;
  }
  
  std::array<fastjet::ClusterSequenceArea*, 4> DijetWorker::loadProperClusterSequence(string cluster_identifier) {
    std::array<fastjet::ClusterSequenceArea*, 4> ret{nullptr, nullptr, nullptr, nullptr};
    // now load the proper cluster sequences
    if (cluster_seq_lead_hard.find(cluster_identifier) == cluster_seq_lead_hard.end() ||
        cluster_seq_lead_match.find(cluster_identifier) == cluster_seq_lead_match.end()) {
      LOG(ERROR) << "could not load cluster sequence for key: " << cluster_identifier;
      return ret;
    }
    
    ret[0] = cluster_seq_lead_hard[cluster_identifier].get();
    ret[2] = cluster_seq_lead_match[cluster_identifier].get();
    
    if (cluster_seq_sub_hard.find(cluster_identifier) == cluster_seq_sub_hard.end())
      ret[1] = cluster_seq_lead_hard[cluster_identifier].get();
    else
      ret[1] = cluster_seq_sub_hard[cluster_identifier].get();
    
    if (cluster_seq_sub_match.find(cluster_identifier) == cluster_seq_sub_match.end())
      ret[3] = cluster_seq_lead_match[cluster_identifier].get();
    else
      ret[3] = cluster_seq_sub_match[cluster_identifier].get();
    
    return ret;
  }
  
  fastjet::PseudoJet DijetWorker::selectLeadingHardJet(const std::vector<fastjet::PseudoJet>& clustered_jets) {
    // for now, only returns the highest momentum jet
    return fastjet::sorted_by_pt(clustered_jets)[0];
  }
  
  
  fastjet::PseudoJet DijetWorker::selectSubLeadingHardJet(const fastjet::PseudoJet& trigger_jet,
                                                          const std::vector<fastjet::PseudoJet>& sublead_candidates,
                                                          double dPhi_range) {
    // select the highest momentum jet that is dPhi > pi - dPhi_cut
    fastjet::Selector recoil_selector = !fastjet::SelectorRectangle(2.1, TMath::Pi() - dPhi_range);
    recoil_selector.set_reference(trigger_jet);
    std::vector<fastjet::PseudoJet> dphi_selected_recoil = fastjet::sorted_by_pt(recoil_selector(sublead_candidates));
    
    if (dphi_selected_recoil.size() == 0)
      return fastjet::PseudoJet();
    
    return dphi_selected_recoil[0];
  }
  
  std::pair<double, double> DijetWorker::estimateBackgroundDensity(const std::vector<fastjet::PseudoJet>& input,
                                                      const JetDef& jet_def) {
    fastjet::JetMedianBackgroundEstimator bkg_est(jet_def.backgroundSelector(),
                                                  jet_def.backgroundJetDef(),
                                                  jet_def.backgroundAreaDef());
    bkg_est.set_particles(jet_def.constituentSelector()(input));
    return {bkg_est.rho(), bkg_est.sigma()};
  }
  
  std::pair<double, double> DijetWorker::subtractBackgroundFromJets(const std::vector<fastjet::PseudoJet>& input,
                                                                    const JetDef& jet_def,
                                                                    const std::vector<fastjet::PseudoJet>& jets,
                                                                    std::vector<fastjet::PseudoJet>& subtracted_jets) {
    fastjet::JetMedianBackgroundEstimator bkg_est(jet_def.backgroundSelector(),
                                                  jet_def.backgroundJetDef(),
                                                  jet_def.backgroundAreaDef());
    bkg_est.set_particles(jet_def.constituentSelector()(input));
    fastjet::Subtractor bkgdSubtractor(&bkg_est);
    bkgdSubtractor.set_use_rho_m(true);
    bkgdSubtractor.set_safe_mass(true);
    subtracted_jets = fastjet::sorted_by_pt(jet_def.jetSelector()(bkgdSubtractor(jets)));
    
    return {bkg_est.rho(), bkg_est.sigma()};
  }
  
  std::pair<unique_ptr<fastjet::JetMedianBackgroundEstimator>,
  unique_ptr<fastjet::Subtractor>> DijetWorker::getBackgroundSubtractor(const std::vector<fastjet::PseudoJet>& input,
                                                                        const JetDef& jet_def) {
    unique_ptr<fastjet::JetMedianBackgroundEstimator> bkg_est =
    make_unique<fastjet::JetMedianBackgroundEstimator>(jet_def.backgroundSelector(),
                                                       jet_def.backgroundJetDef(),
                                                       jet_def.backgroundAreaDef());
    bkg_est->set_particles(jet_def.constituentSelector()(input));
    unique_ptr<fastjet::Subtractor> bkgdSubtractor = make_unique<fastjet::Subtractor>(bkg_est.get());
    bkgdSubtractor->set_use_rho_m(true);
    bkgdSubtractor->set_safe_mass(true);
    return {std::move(bkg_est), std::move(bkgdSubtractor)};
  }
  
  fastjet::PseudoJet DijetWorker::matchJets(const fastjet::PseudoJet& target,
                                            const MatchDef& matchdef,
                                            const std::vector<fastjet::PseudoJet>& candidates) {
    // now match to the hard jets
    fastjet::Selector circle_selector = fastjet::SelectorCircle(matchdef.dR());
    circle_selector.set_reference(target);
    std::vector<fastjet::PseudoJet> tmp = fastjet::sorted_by_pt(circle_selector(candidates));
    if (tmp.size() ==0)
      return fastjet::PseudoJet();
    return tmp[0];
  }
  
  bool DijetWorker::equivalentAreaDefinition(const fastjet::AreaDefinition& a1, const fastjet::AreaDefinition& a2) {
    
    // compare area definitions
    if (a1.area_type() != a2.area_type() ||
        a1.ghost_spec().ghost_maxrap() != a2.ghost_spec().ghost_maxrap() ||
        a1.ghost_spec().repeat() != a2.ghost_spec().repeat() ||
        a1.ghost_spec().ghost_area() != a2.ghost_spec().ghost_area() ||
        a1.ghost_spec().grid_scatter() != a2.ghost_spec().grid_scatter() ||
        a1.ghost_spec().pt_scatter() != a2.ghost_spec().pt_scatter() ||
        a1.ghost_spec().mean_ghost_pt() != a2.ghost_spec().mean_ghost_pt())
      return false;
    
    return true;
  }
  
  bool DijetWorker::equivalentClusterInput(const JetDef& c1, const JetDef& c2) {
    // compare the jet definition
    if (c1.R() != c2.R() ||
        c1.jet_algorithm() != c2.jet_algorithm() ||
        c1.strategy() != c2.strategy() ||
        c1.recombination_scheme() != c2.recombination_scheme())
      return false;
    
    // compare AreaDefinitions
    if (!equivalentAreaDefinition(c1.areaDefinition(), c2.areaDefinition()))
      return false;
    
    // compare constituent selector since I dont want to
    // compare all possible selectors, will only compare
    // description strings, so if they werent constructed
    // identically this will fail, even if they give
    // identical output
    if (c1.constituentSelector().description() != c2.constituentSelector().description())
      return false;
    
    return true;
  }
  
  bool DijetWorker::equivalentBkgEstimationInput(const JetDef& c1, const JetDef& c2) {
    
    // first. compare constituent and bkg jet selector
    if (c1.constituentSelector().description() != c2.constituentSelector().description())
      return false;
    
    if (c1.backgroundSelector().description() != c2.backgroundSelector().description())
      return false;
    
    // compare background area definitions
    if (!equivalentAreaDefinition(c1.backgroundAreaDef(), c2.backgroundAreaDef()))
      return false;
    
    // finally compare background jet definitions
    if (c1.backgroundJetDef().R() != c2.backgroundJetDef().R() ||
        c1.backgroundJetDef().jet_algorithm() != c2.backgroundJetDef().jet_algorithm() ||
        c1.backgroundJetDef().recombination_scheme() != c2.backgroundJetDef().recombination_scheme() ||
        c1.backgroundJetDef().strategy() != c2.backgroundJetDef().strategy())
      return false;
    
    return true;
  }
  
} // namespace dijetcore
