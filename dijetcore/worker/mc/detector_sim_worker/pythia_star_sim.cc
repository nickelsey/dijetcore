#include "dijetcore/worker/mc/detector_sim_worker/pythia_star_sim.h"
#include "dijetcore/lib/assert.h"
#include "dijetcore/lib/logging.h"
namespace dijetcore {

PythiaStarSim::PythiaStarSim()
    : max_events_(0), current_event_(0), pt_hat_min_(20.0), pt_hat_max_(30.0), jet_radius_(0.6),
      alg_(fastjet::antikt_algorithm), track_pt_min_(0.2), track_pt_max_(30.0),
      track_eta_max_(1.0), jet_pt_min_(10), jet_eta_max_(1.0), seed_(1343),
      initialized_(false), eff_mode_(EfficiencyMode::STAR), gauss_mean_(0.0),
      gauss_sigma_(0.4), gauss_smear_("gaussmear", "gaus", -5, 5),
      do_clustering_(true), flat_dis_(0.0, 1.0), gen_(152452) {
  gauss_smear_.SetParameters(1.0, gauss_mean_, gauss_sigma_);
}

bool PythiaStarSim::Initialize(unsigned n_events) {
  max_events_ = n_events;

  // build selectors for fastjet
  track_selector_ = fastjet::SelectorPtMin(track_pt_min_) &&
                    fastjet::SelectorPtMax(track_pt_max_) &&
                    fastjet::SelectorAbsRapMax(track_eta_max_);
  jet_selector_ = fastjet::SelectorPtMin(jet_pt_min_) &&
                  fastjet::SelectorAbsRapMax(track_eta_max_ - jet_radius_);
  jet_def_ = fastjet::JetDefinition(alg_, jet_radius_);

  // now setup pythia generator
  std::string pt_min_string_ =
      "PhaseSpace:pTHatMin = " + std::to_string(pt_hat_min_);
  std::string pt_max_string_ =
      "PhaseSpace:pTHatMax = " + std::to_string(pt_hat_max_);
  std::string qcd_string_ = "HardQCD:all = on";
  std::string e_cm_string_ = "Beams:eCM = 200";
  std::string rnd_command_ = "Random:setSeed = On";
  std::string rnd_seed_ = "Random:seed = " + std::to_string(seed_);
  pythia_.readString(pt_min_string_);
  pythia_.readString(pt_max_string_);
  pythia_.readString(qcd_string_);
  pythia_.readString(e_cm_string_);
  for (auto &str : gen_args_)
    pythia_.readString(str);
  pythia_.readString(rnd_command_);
  pythia_.readString(rnd_seed_);

  bool init_status_ = pythia_.init();

  if (init_status_ == false)
    LOG(ERROR) << "Pythia generator failed init - can't generate events";
  else
    initialized_ = true;
  return init_status_;
}

void PythiaStarSim::Run() {
  // clear results
  truth_particles_.clear();
  det_particles_.clear();
  truth_jets_.clear();
  det_jets_.clear();

  if (initialized_ == false)
    DIJETCORE_THROW(
        "Pythia Generator uninitialized - failed to generate event");

  // try to generate a good event
  int gen_count_ = 0;
  do {
    if (gen_count_ == 100)
      DIJETCORE_THROW(
          "generator failed to generate good event - check settings?");
    if (pythia_.next())
      break;
    gen_count_++;
  } while (true);

  // so we have a good event - translate input to PseudoJets
  ConvertToPseudoJet();

  // if clustering is asked for, run that now
  if (do_clustering_) {
    truth_seq_ =
        make_unique<fastjet::ClusterSequence>(truth_particles_, jet_def_);
    truth_jets_ = fastjet::sorted_by_pt(truth_seq_->inclusive_jets());

    det_seq_ = make_unique<fastjet::ClusterSequence>(det_particles_, jet_def_);
    det_jets_ =
        jet_selector_(fastjet::sorted_by_pt(det_seq_->inclusive_jets()));
  }
}

bool PythiaStarSim::Next() {
  if (max_events_ == 0 || (max_events_ > 0 && current_event_ < max_events_)) {
    Run();
    return true;
  }
  return false;
}

void PythiaStarSim::SetSmearingParams(double mean, double sigma) {
  gauss_smear_.SetParameters(1.0, mean, sigma);
}

void PythiaStarSim::ConvertToPseudoJet() {
  for (size_t i = 0; i < pythia_.event.size(); ++i) {
    if (pythia_.event[i].isFinal() && pythia_.event[i].isVisible()) {
      fastjet::PseudoJet tmp(pythia_.event[i].px(), pythia_.event[i].py(),
                             pythia_.event[i].pz(), pythia_.event[i].e());
      tmp.set_user_index(pythia_.event[i].charge());

      // add to the truth list
      truth_particles_.push_back(tmp);

      // for the detector particles, we have to smear them
      AddDetectorTrack(tmp);
    }
  }
  // apply track kinematic cuts to generator level particles
  det_particles_ = track_selector_(det_particles_);
}

void PythiaStarSim::AddDetectorTrack(fastjet::PseudoJet &part) {
  double eff;
  double prob;
  switch (eff_mode_) {
  case EfficiencyMode::None:
    det_particles_.push_back(part);
    break;

  case EfficiencyMode::GaussianSmearing:
    part *= 1.0 + gauss_smear_.GetRandom();
    det_particles_.push_back(part);
    break;

  case EfficiencyMode::STAR:
    if (part.user_index() == 0) {
      det_particles_.push_back(part);
      break;
    }

    prob = flat_dis_(gen_);
    eff = efficiency_.ppEff(part.pt(), part.eta());
    if (prob < eff) {
      part *= efficiency_.smearPt(part.pt());
      det_particles_.push_back(part);
    }
    break;

  default:
    DIJETCORE_THROW("Should never get here");
  }
}

} // namespace dijetcore
