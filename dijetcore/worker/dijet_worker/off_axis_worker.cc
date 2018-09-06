#include "dijetcore/worker/dijet_worker/off_axis_worker.h"

#include "dijetcore/util/data/vector_conversion.h"
#include "dijetcore/util/fastjet/selector_compare.h"
#include "dijetcore/lib/logging.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"

#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

namespace dijetcore {
  
  OffAxisWorker::OffAxisWorker() {
    cent_min_ = 0;
    cent_max_ = 8;
  }
  
  OffAxisWorker::~OffAxisWorker() {
    
  }
  
  void OffAxisWorker::SetCentralityRange(int centLow, int centHigh) {
    LOG(INFO) << "Initializing centrality range for OffAxisWorker";
    if (centLow >= 0 && centLow < 9)
      centLow = centLow;
    else
      centLow = 0;
    if (centHigh < 0 || centHigh >= 9)
      centHigh = 8;
    if (centHigh < centLow) {
      LOG(ERROR) << "centrality high bound is less than centrality low bound: exiting";
      return;
    }
    
    cent_min_ = centLow;
    cent_max_ = centHigh;
    int nbins = centHigh - centLow + 1;
    pre_selected_events_ = std::vector<std::vector<unsigned>>(nbins, std::vector<unsigned>());
    current_event_ = std::vector<unsigned>(nbins, 0);
    
    ReadEvent(0);
    int centrality = GetCentrality();
    if (centrality >= cent_min_ && centrality <= cent_max_) {
      int idx =  centrality - cent_min_;
      pre_selected_events_[idx].push_back(GetNOfCurrentEvent());
    }
    while (NextEvent()) {
      int centrality = GetCentrality();
      if (centrality >= cent_min_ && centrality <= cent_max_) {
        int idx =  centrality - cent_min_;
        pre_selected_events_[idx].push_back(GetNOfCurrentEvent());
      }
    }
  }
  
  std::unordered_map<std::string, unique_ptr<OffAxisOutput>>& OffAxisWorker::Run(DijetWorker& worker, int centrality) {

    // clear our output
    off_axis_result_.clear();
    cluster_sequence_.clear();
    
    // see if we need to do anything at all -
    // if there was no match in the Dijetworker, we exit out
    bool found_match = false;
    for (auto& result : worker.Dijets())
      if (result.second->found_match)
        found_match = true;
    
    if (found_match == false)
      return off_axis_result_;

    if (!LoadNextEvent(centrality)) {
      LOG(ERROR) << "OffAxisWorker could not find an event that satisfies all event cuts";
      return off_axis_result_;
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
      
      // have to remove bkg contribution from input tracks that are in the same kinematic window as the reference jet
      double max_pt_lead = ExtractDoubleFromSelector(def->lead->InitialJetDef().ConstituentSelector(), "pt >=");
      double max_pt_sub = ExtractDoubleFromSelector(def->sub->InitialJetDef().ConstituentSelector(), "pt >=");
      fastjet::Selector max_pt_sel_lead = fastjet::SelectorPtMax(max_pt_lead);
      fastjet::Selector max_pt_sel_sub = fastjet::SelectorPtMax(max_pt_sub);
      
      std::vector<fastjet::PseudoJet> lead_particles = max_pt_sel_lead(primary_particles);
      lead_particles.push_back(lead_hard_jet);
      lead_particles.push_back(sub_hard_jet);
      std::vector<fastjet::PseudoJet> sub_particles = max_pt_sel_sub(primary_particles);
      sub_particles.push_back(lead_hard_jet);
      sub_particles.push_back(sub_hard_jet);
      
      auto lead_matched_jet = RunCluster(def->lead, lead_particles, lead_hard_jet);
      auto sub_matched_jet = RunCluster(def->sub, sub_particles, sub_hard_jet);
      
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
  
    std::unique_ptr<fastjet::ClusterSequenceArea>
    cluster = make_unique<fastjet::ClusterSequenceArea>(def->MatchedJetDef().ConstituentSelector()(input),
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
  
  bool OffAxisWorker::LoadNextEvent(int centrality) {
    if (pre_selected_events_.size() == 0) {
      while (NextEvent()) {
        if (GetCentrality() == centrality)
          return true;
      }
      ReadEvent(0);
      if (GetCentrality() == centrality)
        return true;
      while (NextEvent()) {
        if (GetCentrality() == centrality)
          return true;
      }
      return false;
    }
    
    
    if (centrality < cent_min_ || centrality > cent_max_) {
      LOG(ERROR) << "submitted event outside of specified centrality range: exiting";
      return false;
    }
    int cent_bin = centrality - cent_min_;
    auto& events = pre_selected_events_.at(cent_bin);
    
    if (events.size() == 0) {
      LOG(ERROR) << "Do not have any events in specified centrality bin: exiting";
      return false;
    }
    
    unsigned& current_idx = current_event_.at(cent_bin);
    
    if (current_idx == events.size())
      current_idx = 0;
    
    ReadEvent(current_idx);
    current_idx++;
    return true;
  }
  
} // namespace dijetcore
