#ifndef DIJETCORE_WORKER_DIJET_WORKER_DIJET_WORKER_H
#define DIJETCORE_WORKER_DIJET_WORKER_DIJET_WORKER_H

#include <unordered_map>
#include <array>

#include "dijetcore/lib/memory.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_matrix.h"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

namespace dijetcore {
  
  // struct used to return the results of a dijet clustering
  struct ClusterOutput {
    
    ClusterOutput() :
    dijet_def(nullptr),
    found_lead(false),
    found_sublead(false),
    found_match(false) { }
    
    ClusterOutput(DijetDefinition* def) :
    dijet_def(def),
    found_lead(false),
    found_sublead(false),
    found_match(false) { }
    
    inline bool FoundDijet() {return found_match;}
    
    DijetDefinition* dijet_def;
    
    bool found_lead;
    bool found_sublead;
    bool found_match;
    
    fastjet::PseudoJet lead_hard;
    double lead_hard_rho;
    double lead_hard_sigma;
    
    fastjet::PseudoJet sublead_hard;
    double sublead_hard_rho;
    double sublead_hard_sigma;
    
    fastjet::PseudoJet lead_match;
    double lead_match_rho;
    double lead_match_sigma;
    
    fastjet::PseudoJet sublead_match;
    double sublead_match_rho;
    double sublead_match_sigma;
  };
  
  class DijetWorker : public DijetMatrix {
  public:
    DijetWorker();
    
    // single entry construction
    DijetWorker(fastjet::JetAlgorithm jet_alg_in,
                double lead_jet_pt_in = 20.0,
                double lead_jet_R_in = 0.4,
                double lead_match_jet_R_in = 0.4,
                double sub_jet_pt_in = 10.0,
                double sub_jet_R_in = 0.4,
                double sub_match_jet_R_in = 0.4,
                double const_lead_pt_init_in = 2.0,
                double const_lead_pt_match_in = 0.2,
                double const_sub_pt_init_in = 2.0,
                double const_sub_pt_match_in = 0.2,
                double eta_in = 1.0);
    
    
    // set construction
    DijetWorker(std::set<fastjet::JetAlgorithm> jet_alg_in,
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
                std::set<double> eta_in);
    
    DijetWorker(const DijetMatrix& rhs);
    
    // clusters a set of input particles (fastjet::PseudoJets)
    // and gives as output a dictionary containing the dijet pairs
    // of all DijetDefinitions where a dijet pair was found
    std::unordered_map<std::string, unique_ptr<ClusterOutput>>&
    Run(const std::vector<fastjet::PseudoJet>& input);
    
    // second method to access the dictionary, outside
    // of calling run.
    std::unordered_map<std::string, unique_ptr<ClusterOutput>>& Dijets() {return cluster_result;}
    
  private:
    
    // clears all of the containers
    void ClearResults();
    
    // Internal step for running the clustering for each DijetDefinition
    bool RunClustering(const std::vector<fastjet::PseudoJet>& input, const DijetDefinition& def,
                       string cluster_identifier, double min_lead_pt, double min_sub_pt);
    
    // Internal step for identifying the correct cluster sequences for
    // a specific DijetDefinition
    std::array<fastjet::ClusterSequenceArea*, 4> LoadProperClusterSequence(string cluster_identifier);
    
    // Used to select leading and subleading hard jets
    // currently only supports selecting highest pT jet as leading -
    // intend to implement trigger requirement later
    fastjet::PseudoJet SelectLeadingHardJet(const std::vector<fastjet::PseudoJet>& clustered_jets);
    fastjet::PseudoJet SelectSubLeadingHardJet(const fastjet::PseudoJet& trigger_jet,
                                               const std::vector<fastjet::PseudoJet>& sublead_candidates,
                                               double dPhi_range);
    
    // function to estimate background energy density
    std::pair<double, double> EstimateBackgroundDensity(const std::vector<fastjet::PseudoJet>& input,
                                                        const JetDef& jet_def);
    
    // returns a list of background subtracted jets
    // returns rho and sigma the same was EstimateBackgroundDensity
    std::pair<double, double> SubtractBackgroundFromJets(const std::vector<fastjet::PseudoJet>& input,
                                                         const JetDef& jet_def,
                                                         const std::vector<fastjet::PseudoJet>& jets,
                                                         std::vector<fastjet::PseudoJet>& subtracted_jets);
    
    // or get the background subtractor itself
    std::pair<unique_ptr<fastjet::JetMedianBackgroundEstimator>,
              unique_ptr<fastjet::Subtractor>> GetBackgroundSubtractor(const std::vector<fastjet::PseudoJet>& input,
                                                                       const JetDef& jet_def);
    
    
    // matches a reclustered jet to the selected leading/subleading hard jet
    fastjet::PseudoJet MatchJets(const fastjet::PseudoJet& target,
                                 const MatchDef& matchdef,
                                 const std::vector<fastjet::PseudoJet>& candidates);
    
    
    // the next three functions are used to try and cut down
    // on run time by only clustering leading/subleading jet
    // separately when necessary - that is, when the input
    // constituents will be different, when the jet definition
    // is not identical, or when the area definitions don't match
    
    // checks to see if two area definitions are equivalent,
    // used by EquivalentClusterInput()
    bool EquivalentAreaDefinition(const fastjet::AreaDefinition& a1, const fastjet::AreaDefinition& a2);
    
    // used to decide if two clusterings would give the same output
    // to reduce computation time. Comparison of ConstituentSelector
    // AreaDefinition and JetDefinition
    bool EquivalentClusterInput(const JetDef& c1, const JetDef& c2);
    
    // same for background estimation
    bool EquivalentBkgEstimationInput(const JetDef& c1, const JetDef& c2);
    
    // used to keep cluster sequences in scope, in case user wants to
    // access constituents
    std::unordered_map<std::string, unique_ptr<fastjet::ClusterSequenceArea>> cluster_seq_lead_hard;
    std::unordered_map<std::string, unique_ptr<fastjet::ClusterSequenceArea>> cluster_seq_sub_hard;
    std::unordered_map<std::string, unique_ptr<fastjet::ClusterSequenceArea>> cluster_seq_lead_match;
    std::unordered_map<std::string, unique_ptr<fastjet::ClusterSequenceArea>> cluster_seq_sub_match;
    
    std::unordered_map<std::string, unique_ptr<ClusterOutput>> cluster_result;
    
    
  };
  
} // namespace dijetcore

#endif // DIJETCORE_WORKER_DIJET_WORKER_DIJET_WORKER_H
