#ifndef DIJETCORE_WORKER_MC_DETECTOR_SIM_WORKER_PYTHIA_STAR_SIM_H
#define DIJETCORE_WORKER_MC_DETECTOR_SIM_WORKER_PYTHIA_STAR_SIM_H

#include "dijetcore/lib/memory.h"
#include "dijetcore/util/data/efficiency/run6_eff.h"

#include <random>
#include <string>
#include <vector>

#include "Pythia8/Pythia.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/Selector.hh"

#include "TF1.h"

namespace dijetcore {

class PythiaStarSim {
public:
  enum class EfficiencyMode { None, GaussianSmearing, STAR };

  PythiaStarSim();

  ~PythiaStarSim(){};

  bool Initialize(unsigned n_events = 0);

  // sets if the generator should run fastjet clustering
  // default on
  void DoClustering(bool flag) { do_clustering_ = flag; }

  // generates a new event
  void Run();

  // will generate events until max_events has been reached
  bool Next();

  // returns particles for current event
  std::vector<fastjet::PseudoJet> &GenParticles() { return truth_particles_; }
  std::vector<fastjet::PseudoJet> &DetectorParticles() {
    return det_particles_;
  }

  // returns jets for current event
  std::vector<fastjet::PseudoJet> &GenJets() { return truth_jets_; }
  std::vector<fastjet::PseudoJet> &DetectorJets() { return det_jets_; }

  void SetEfficiencyMode(EfficiencyMode mode) { eff_mode_ = mode; }
  EfficiencyMode GetEfficiencyMode() const { return eff_mode_; }
  void SetSmearingParams(double mean, double sigma);

  void SetRandomSeed(int seed) { seed_ = seed; }

  void SetPtHatRange(double min, double max) {
    pt_hat_min_ = min;
    pt_hat_max_ = max;
  };
  void AddSettingString(const std::string &str) { gen_args_.push_back(str); }

  void SetJetPtMin(double pt) { jet_pt_min_ = pt; }
  double JetPtMin() const { return jet_pt_min_; }

  void SetTrackPtMin(double pt) { track_pt_min_ = pt; }
  double TrackPtMin() const { return track_pt_min_; }

  void SetTrackPtMax(double pt) { track_pt_max_ = pt; }
  double TrackPtMax() const { return track_pt_max_; }

  void SetTrackEtaMax(double eta) { track_eta_max_ = fabs(eta); }
  double TrackEtaMax() const { return track_eta_max_; }

  // if this is zero at initialization time, it will be inferred
  // from jet radius and track eta
  void SetJetEtaMax(double eta) { jet_eta_max_ = fabs(eta); }
  double JetEtaMax() const { return jet_eta_max_; }

  void SetJetRadius(double rad) { jet_radius_ = rad; }
  double JetRadius() const { return jet_radius_; }

  void SetJetAlgorithm(fastjet::JetAlgorithm alg) { alg_ = alg; }
  fastjet::JetAlgorithm JetAlgorithm() const { return alg_; }

  Pythia8::Pythia &Pythia() { return pythia_; }

  std::string GetStatus();
  unsigned CurrentEvent() { return current_event_; }

private:
  // number of events to generate
  unsigned max_events_;
  // current event number
  unsigned current_event_;

  // Pythia generator settings
  double pt_hat_min_;
  double pt_hat_max_;
  std::vector<std::string> gen_args_;

  // analysis settings
  double jet_radius_;
  fastjet::JetAlgorithm alg_;

  double track_pt_min_;
  double track_pt_max_;
  double track_eta_max_;
  double jet_pt_min_;
  double jet_eta_max_;

  // generator
  int seed_;
  Pythia8::Pythia pythia_;
  bool initialized_;

  // efficiency settings
  EfficiencyMode eff_mode_;
  double gauss_mean_;
  double gauss_sigma_;
  Run6Eff efficiency_;
  TF1 gauss_smear_;

  bool do_clustering_;

  fastjet::JetDefinition jet_def_;
  fastjet::Selector track_selector_;
  fastjet::Selector jet_selector_;

  unique_ptr<fastjet::ClusterSequence> det_seq_;
  unique_ptr<fastjet::ClusterSequence> truth_seq_;
  std::vector<fastjet::PseudoJet> det_particles_;
  std::vector<fastjet::PseudoJet> truth_particles_;
  std::vector<fastjet::PseudoJet> det_jets_;
  std::vector<fastjet::PseudoJet> truth_jets_;

  // internally used functions
  void ConvertToPseudoJet();
  void AddDetectorTrack(fastjet::PseudoJet &part);

  // random number distribution
  std::uniform_real_distribution<double> flat_dis_;
  std::mt19937 gen_;
};

} // namespace dijetcore

#endif // DIJETCORE_WORKER_MC_DETECTOR_SIM_WORKER_PYTHIA_STAR_SIM_H