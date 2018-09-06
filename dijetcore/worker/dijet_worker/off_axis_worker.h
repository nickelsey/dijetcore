#ifndef DIJETCORE_WORKER_DIJET_WORKER_OFF_AXIS_WORKER_H
#define DIJETCORE_WORKER_DIJET_WORKER_OFF_AXIS_WORKER_H

#include "dijetcore/worker/dijet_worker/dijet_worker.h"

#include "TStarJetPicoReader.h"

class DijetWorker;

namespace dijetcore {
  
  struct OffAxisOutput {
    OffAxisOutput();
    OffAxisOutput(const fastjet::PseudoJet& j1, const fastjet::PseudoJet& j2) : lead_jet(j1), sub_jet(j2) {}
    
    fastjet::PseudoJet lead_jet;
    double lead_jet_rho;
    double lead_jet_sigma;
    
    fastjet::PseudoJet sub_jet;
    double sub_jet_rho;
    double sub_jet_sigma;
  };
  
  class OffAxisWorker :public TStarJetPicoReader {
  public:
    OffAxisWorker();
    ~OffAxisWorker();
    
    // if you know which centralities you are interested in, then you can set them beforehand
    // this will speed up the process (-1 for open ended)
    // don't call this function multiple times - it has to run over the full tree
    void SetCentralityRange(int centLow, int centHigh);
    
    std::unordered_map<std::string, unique_ptr<OffAxisOutput>>& Run(DijetWorker& worker, int centrality);
    
    std::unordered_map<std::string, unique_ptr<OffAxisOutput>>& OffAxisResult() {return off_axis_result_;}
    
    int GetCentrality();
    
  private:
    
    bool LoadNextEvent(int centrality);
    
    std::pair<fastjet::PseudoJet, std::pair<double, double>>
    RunCluster(MatchDef* def, const std::vector<fastjet::PseudoJet>& input, const fastjet::PseudoJet& reference);
    
    std::vector<std::vector<unsigned>> pre_selected_events_;
    std::vector<unsigned> current_event_;
    unsigned cent_min_;
    unsigned cent_max_;
    
    std::vector<unique_ptr<fastjet::ClusterSequenceArea>> cluster_sequence_;
    std::unordered_map<std::string, unique_ptr<OffAxisOutput>> off_axis_result_;
    
  };
  
} // namespace dijetcore

#endif // DIJETCORE_WORKER_DIJET_WORKER_OFF_AXIS_WORKER_H
