#ifndef DIJETCORE_WORKER_DIJET_WORKER_ETA_CONE_WORKER_H
#define DIJETCORE_WORKER_DIJET_WORKER_ETA_CONE_WORKER_H

#include "dijetcore/worker/dijet_worker/dijet_worker.h"

#include <random>

namespace dijetcore {

struct EtaConeOutput {
  EtaConeOutput() {};
  EtaConeOutput(const fastjet::PseudoJet &j1, const fastjet::PseudoJet &j2)
      : lead_jet(j1), sub_jet(j2) {}

  fastjet::PseudoJet lead_jet;
  double lead_jet_rho = -1.0;
  double lead_jet_sigma = -1.0;

  fastjet::PseudoJet sub_jet;
  double sub_jet_rho = -1.0;
  double sub_jet_sigma = -1.0;
};

class EtaConeWorker {
public:
  EtaConeWorker();
  ~EtaConeWorker();

  std::unordered_map<std::string, unique_ptr<EtaConeOutput>> &
  run(DijetWorker &worker, std::vector<fastjet::PseudoJet> &input);

  std::unordered_map<std::string, unique_ptr<EtaConeOutput>> &EtaConeResult() {
    return eta_cone_result_;
  }

private:

  // finds an acceptable eta for an eta-offset cone that must be at least 2*r
  // away from the jet axis, as well as R away from const_max_eta.
  // returns a pair with the first being a boolean representing success/failure,
  // and a double representing the chosen value  
  std::pair<bool, double> findEta(const fastjet::PseudoJet &reference, double r,
                                  double const_max_eta);

  std::unordered_map<std::string, unique_ptr<EtaConeOutput>> eta_cone_result_;

  // generator for selecting a random eta
  std::mt19937 gen_;

};

} // namespace dijetcore

#endif // DIJETCORE_WORKER_DIJET_WORKER_ETA_CONE_WORKER_H
