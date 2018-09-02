#include "dijetcore/worker/dijet_worker/off_axis_worker.h"

#include "dijetcore/util/data/vector_conversion.h"
#include "dijetcore/util/fastjet/selector_compare.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"

#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

namespace dijetcore {
  
  OffAxisWorker::OffAxisWorker() {
    
  }
  
  OffAxisWorker::~OffAxisWorker() {
    
  }
  
  std::unordered_map<std::string, unique_ptr<OffAxisOutput>>& OffAxisWorker::Run(DijetWorker& worker, int centrality) {
    
    // clear our output
    off_axis_result_.clear();
    cluster_sequence_.clear();
    
    // find an event with the proper centrality
    if (!NextEvent())
      ReadEvent(0);
    
    int cent_bin = GetCentrality();
    
    while (cent_bin != centrality) {
      if (!NextEvent())
        ReadEvent(0);
      cent_bin = GetCentrality();
    }
    
    // load an event in the proper centrality class
    
    // get the input container
    auto& dijets = worker.Dijets();
    
    for (auto& dijet : dijets) {
      
      if (!dijet.second->found_match)
        continue;
      
      auto& key = dijet.first;
      auto& def = dijet.second->dijet_def;
      auto& lead_hard_jet = dijet.second->lead_hard;
      auto& sub_hard_jet = dijet.second->sublead_hard;
      
      // build our new input list
      std::vector<fastjet::PseudoJet> primary_particles;
      dijetcore::ConvertTStarJetVector(GetOutputContainer(), primary_particles);
      primary_particles.push_back(lead_hard_jet);
      primary_particles.push_back(sub_hard_jet);
      
      auto lead_matched_jet = RunCluster(def->lead, primary_particles, lead_hard_jet);
      auto sub_matched_jet = RunCluster(def->sub, primary_particles, sub_hard_jet);
      
      unique_ptr<OffAxisOutput> res = make_unique<OffAxisOutput>(lead_matched_jet.first, sub_matched_jet.first);
      res->lead_jet_rho = lead_matched_jet.second.first;
      res->lead_jet_sigma = lead_matched_jet.second.second;
      res->sub_jet_rho = sub_matched_jet.second.first;
      res->sub_jet_sigma = sub_matched_jet.second.second;
      
      off_axis_result_[key] = std::move(res);
    }
    
    return off_axis_result_;
  }
  
  int OffAxisWorker::GetCentrality() {
    static const std::vector<int> cent_bins_ = {485, 399, 269, 178, 114, 69, 39, 21, 10};
    int cent_bin = -1;
    for (int i = 0; i < cent_bins_.size(); ++i)
      if (GetEvent()->GetHeader()->GetProperReferenceMultiplicity() > cent_bins_[i]) {
        cent_bin = i;
        break;
      }
    return cent_bin;
  }
  
  std::pair<fastjet::PseudoJet, std::pair<double, double>>
  OffAxisWorker::RunCluster(MatchDef* def, const std::vector<fastjet::PseudoJet>& input, const fastjet::PseudoJet& reference) {
  
    // have to remove bkg contribution from input tracks that are in the same kinematic window as the reference jet
    double max_pt = ExtractDoubleFromSelector(def->InitialJetDef().ConstituentSelector(), "pt >=");
    fastjet::Selector max_pt_sel = fastjet::SelectorPtMax(max_pt);
    
    std::unique_ptr<fastjet::ClusterSequenceArea>
    cluster = make_unique<fastjet::ClusterSequenceArea>(def->MatchedJetDef().ConstituentSelector()(max_pt_sel(input)),
                                                        def->MatchedJetDef(),
                                                        def->MatchedJetDef().AreaDefinition());
    
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(def->MatchedJetDef().JetSelector()(cluster->inclusive_jets()));
  
    // do bkg subtraction
    fastjet::JetMedianBackgroundEstimator bkg_est(def->MatchedJetDef().BackgroundSelector(),
                                                  def->MatchedJetDef().BackgroundJetDef(),
                                                  def->MatchedJetDef().BackgroundAreaDef());
    bkg_est.set_particles(def->MatchedJetDef().ConstituentSelector()(input));
    fastjet::Subtractor bkgdSubtractor(&bkg_est);
    std::vector<fastjet::PseudoJet> subtracted_jets = fastjet::sorted_by_pt(bkgdSubtractor(jets));
  
    std::pair<double, double> bkg_stats = {bkg_est.rho(), bkg_est.sigma()};
    
    // find matched jet
    fastjet::Selector circle_sel = fastjet::SelectorCircle(def->InitialJetDef().R());
    circle_sel.set_reference(reference);
    
    subtracted_jets = fastjet::sorted_by_pt(circle_sel(subtracted_jets));
  
    // insert cluster sequence
    cluster_sequence_.push_back(std::move(cluster));
    
    if (subtracted_jets.size() != 0)
      return {subtracted_jets[0], bkg_stats};
    else
      return {fastjet::PseudoJet(), bkg_stats};
    
  }
  
} // namespace dijetcore
