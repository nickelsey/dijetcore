#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <random>
#include <exception>

#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/vector_conversion.h"
#include "dijetcore/util/data/trigger_lookup.h"
#include "dijetcore/worker/dijet_worker/dijet_worker.h"
#include "dijetcore/worker/dijet_worker/off_axis_worker.h"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile2D.h"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"


using std::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(offAxisInput, "", "input file for off-axis bkg effect estimation");
DIJETCORE_DEFINE_string(outputDir, "tmp", "directory for output");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(towList, "", "bad tower list");
DIJETCORE_DEFINE_string(triggers, "y7ht", "trigger selection");
DIJETCORE_DEFINE_string(constEta, "1.0", "constitutent eta cuts (comma separated)");
DIJETCORE_DEFINE_string(leadConstPt, "2.0", "leading constituent pT cut (comma separated)");
DIJETCORE_DEFINE_string(subConstPt, "2.0", "subleading constituent pT cut (comma separated)");
DIJETCORE_DEFINE_string(leadConstPtMatch, "0.2", "matched leading constituent pT cut (comma separated)");
DIJETCORE_DEFINE_string(subConstPtMatch, "0.2", "matched subleading constituent pT cut (comma separated)");
DIJETCORE_DEFINE_string(leadR, "0.4", "leading jet R (comma separated");
DIJETCORE_DEFINE_string(leadMatchR, "0.4", "leading matched jet R (comma separated");
DIJETCORE_DEFINE_string(subR, "0.4", "subleading jet R (comma separated");
DIJETCORE_DEFINE_string(subMatchR, "0.4", "subleading matched jet R (comma separated");
DIJETCORE_DEFINE_string(leadJetPt, "20.0", "leading jet pT cut (comma separated)");
DIJETCORE_DEFINE_string(subJetPt, "10.0", "subleading jet pT cut (comma separated)");
DIJETCORE_DEFINE_bool(forceConstituentPtEquality, true, "Only use DijetDefinitions where pT const is equal in leading/subleading jets");
DIJETCORE_DEFINE_bool(forceConstituentEtaEquality, true, "Only use DijetDefinitions where eta const is equal in leading/subleading jets");
DIJETCORE_DEFINE_bool(forceJetResolutionEquality, true, "Only use DijetDefinitions where leading/subleading R are equivalent");
DIJETCORE_DEFINE_bool(forceMatchJetResolutionEquality, false, "Only use DijetDefinitions where initial and matched R are equivalent");

int main(int argc, char* argv[]) {
  
  string usage = "Run 7 differential di-jet imbalance analysis routine";
  
  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);
  
  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(FLAGS_input)) {
    LOG(ERROR) << "input file does not exist: " << FLAGS_input;;
    return 1;
  }
  
  // first, build our input chain
  TChain* chain = dijetcore::NewChainFromInput(FLAGS_input);
  
  // and a chain for the off axis worker if it exists
  TChain* off_axis_chain = nullptr;
  if (boost::filesystem::exists(FLAGS_offAxisInput))
    off_axis_chain = dijetcore::NewChainFromInput(FLAGS_offAxisInput);
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);
  
  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outputDir + "/" + FLAGS_name + dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // initialize the reader
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  dijetcore::InitReaderWithDefaults(reader, chain, FLAGS_towList);
  
  // initialize the OffAxisWorker
  dijetcore::OffAxisWorker* off_axis_worker = nullptr;
  if (off_axis_chain != nullptr) {
    off_axis_worker = new dijetcore::OffAxisWorker();
    dijetcore::InitReaderWithDefaults(dynamic_cast<TStarJetPicoReader*>(off_axis_worker), off_axis_chain, FLAGS_towList);
  }
  
  // get the trigger IDs that will be used
  std::set<unsigned> triggers = dijetcore::GetTriggerIDs(FLAGS_triggers);
  
  LOG(INFO) << "taking triggers: " << FLAGS_triggers << " for analysis";
  LOG(INFO) << "trigger ids: ";
  for (auto i : triggers)
    LOG(INFO) << i << " ";
  LOG(INFO);
  
  // parse jetfinding variables
  // --------------------------
  
  // first, hard code the algorithm to be anti-kt
  std::set<fastjet::JetAlgorithm> alg{fastjet::antikt_algorithm};
  
  // constituent range
  std::set<double> const_eta              = dijetcore::ParseArgString<double>(FLAGS_constEta);
  
  // leading jet
  std::set<double> lead_const_hard_pt     = dijetcore::ParseArgString<double>(FLAGS_leadConstPt);
  std::set<double> lead_const_match_pt    = dijetcore::ParseArgString<double>(FLAGS_leadConstPtMatch);
  std::set<double> lead_R                 = dijetcore::ParseArgString<double>(FLAGS_leadR);
  std::set<double> lead_R_match           = dijetcore::ParseArgString<double>(FLAGS_leadMatchR);
  std::set<double> lead_hard_pt           = dijetcore::ParseArgString<double>(FLAGS_leadJetPt);
  
  // subleading jet
  std::set<double> sublead_const_hard_pt  = dijetcore::ParseArgString<double>(FLAGS_subConstPt);
  std::set<double> sublead_const_match_pt = dijetcore::ParseArgString<double>(FLAGS_subConstPtMatch);
  std::set<double> sublead_R              = dijetcore::ParseArgString<double>(FLAGS_subR);
  std::set<double> sublead_R_match        = dijetcore::ParseArgString<double>(FLAGS_subMatchR);
  std::set<double> sublead_hard_pt        = dijetcore::ParseArgString<double>(FLAGS_subJetPt);
  
  // here we can initialize the worker
  LOG(INFO) << "initializing worker...";
  dijetcore::DijetWorker worker(alg, lead_hard_pt, lead_R, lead_R_match, sublead_hard_pt, sublead_R,
                                sublead_R_match, lead_const_hard_pt, lead_const_match_pt,
                                sublead_const_hard_pt, sublead_const_match_pt, const_eta);
  worker.ForceConstituentPtEquality(FLAGS_forceConstituentPtEquality);
  worker.ForceConstituentEtaEquality(FLAGS_forceConstituentEtaEquality);
  worker.ForceJetResolutionEquality(FLAGS_forceJetResolutionEquality);
  worker.ForceMatchJetResolutionEquality(FLAGS_forceMatchJetResolutionEquality);
  worker.Initialize();
  
  std::set<std::string> keys = worker.Keys();
  
  // define the centrality cuts for year 7
  std::vector<unsigned> cent_boundaries = {485, 399, 269, 178, 114, 69, 39, 21, 10};
  
  for (auto key : keys)
    LOG(INFO) << key;
  
  // create an output tree for each definition
  // -----------------------------------------
  
  std::unordered_map<std::string, std::shared_ptr<TTree>> trees;
  
  // and the necessary branches
  std::unordered_map<std::string, int> run_id_dict;
  std::unordered_map<std::string, int> event_id_dict;
  std::unordered_map<std::string, double> vz_dict;
  std::unordered_map<std::string, int> refmult_dict;
  std::unordered_map<std::string, int> grefmult_dict;
  std::unordered_map<std::string, double> refmultcorr_dict;
  std::unordered_map<std::string, double> grefmultcorr_dict;
  std::unordered_map<std::string, int> cent_dict;
  std::unordered_map<std::string, double> zdcrate_dict;
  std::unordered_map<std::string, double> reactionplane_dict;
  std::unordered_map<std::string, int> nglobal_dict;
  std::unordered_map<std::string, int> npart_dict;
  std::unordered_map<std::string, TLorentzVector> lead_hard_jet_dict;
  std::unordered_map<std::string, TH1D*> lead_hard_jet_const_pt_dict;
  std::unordered_map<std::string, TH1D*> lead_hard_jet_const_dr_dict;
  std::unordered_map<std::string, int> lead_hard_jet_nconst_dict;
  std::unordered_map<std::string, double> lead_hard_rho_dict;
  std::unordered_map<std::string, double> lead_hard_sigma_dict;
  std::unordered_map<std::string, double> lead_hard_area_dict;
  std::unordered_map<std::string, TLorentzVector> lead_match_jet_dict;
  std::unordered_map<std::string, TH1D*> lead_match_jet_const_pt_dict;
  std::unordered_map<std::string, TH1D*> lead_match_jet_const_dr_dict;
  std::unordered_map<std::string, int> lead_match_jet_nconst_dict;
  std::unordered_map<std::string, double> lead_match_rho_dict;
  std::unordered_map<std::string, double> lead_match_sigma_dict;
  std::unordered_map<std::string, double> lead_match_area_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_hard_jet_dict;
  std::unordered_map<std::string, TH1D*> sublead_hard_jet_const_pt_dict;
  std::unordered_map<std::string, TH1D*> sublead_hard_jet_const_dr_dict;
  std::unordered_map<std::string, int> sublead_hard_jet_nconst_dict;
  std::unordered_map<std::string, double> sublead_hard_rho_dict;
  std::unordered_map<std::string, double> sublead_hard_sigma_dict;
  std::unordered_map<std::string, double> sublead_hard_area_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_match_jet_dict;
  std::unordered_map<std::string, TH1D*> sublead_match_jet_const_pt_dict;
  std::unordered_map<std::string, TH1D*> sublead_match_jet_const_dr_dict;
  std::unordered_map<std::string, int> sublead_match_jet_nconst_dict;
  std::unordered_map<std::string, double> sublead_match_rho_dict;
  std::unordered_map<std::string, double> sublead_match_sigma_dict;
  std::unordered_map<std::string, double> sublead_match_area_dict;
  std::unordered_map<std::string, TLorentzVector> lead_off_axis_jet_dict;
  std::unordered_map<std::string, int> lead_off_axis_nconst_dict;
  std::unordered_map<std::string, double> lead_off_axis_rho_dict;
  std::unordered_map<std::string, double> lead_off_axis_sigma_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_off_axis_jet_dict;
  std::unordered_map<std::string, int> sublead_off_axis_nconst_dict;
  std::unordered_map<std::string, double> sublead_off_axis_rho_dict;
  std::unordered_map<std::string, double> sublead_off_axis_sigma_dict;
  // fill the maps first, so that they don't decide to resize/move themselves
  // after branch creation...
  for (auto key : keys) {
    run_id_dict.insert({key, 0});
    event_id_dict.insert({key, 0});
    vz_dict.insert({key, 0});
    refmult_dict.insert({key, 0});
    grefmult_dict.insert({key, 0});
    refmultcorr_dict.insert({key, 0});
    grefmultcorr_dict.insert({key, 0});
    cent_dict.insert({key,0});
    zdcrate_dict.insert({key, 0});
    reactionplane_dict.insert({key, 0});
    nglobal_dict.insert({key, 0});
    npart_dict.insert({key, 0});
    lead_hard_jet_dict.insert({key, TLorentzVector()});
    lead_hard_jet_const_pt_dict.insert({key, new TH1D(dijetcore::MakeString(key, "jlconstpt").c_str(),
                                                      "", 90, 0, 30)});
    lead_hard_jet_const_dr_dict.insert({key, new TH1D(dijetcore::MakeString(key, "jlconstdr").c_str(),
                                                      "", 100, 0.0, 1.0)});
    lead_hard_jet_nconst_dict.insert({key, 0});
    lead_hard_rho_dict.insert({key, 0});
    lead_hard_sigma_dict.insert({key, 0});
    lead_hard_area_dict.insert({key, 0});
    lead_match_jet_dict.insert({key, TLorentzVector()});
    lead_match_jet_const_pt_dict.insert({key, new TH1D(dijetcore::MakeString(key, "jlmconstpt").c_str(),
                                                      "", 90, 0, 30)});
    lead_match_jet_const_dr_dict.insert({key, new TH1D(dijetcore::MakeString(key, "jlmconstdr").c_str(),
                                                      "", 100, 0.0, 1.0)});
    lead_match_jet_nconst_dict.insert({key, 0});
    lead_match_rho_dict.insert({key, 0});
    lead_match_sigma_dict.insert({key, 0});
    lead_match_area_dict.insert({key, 0});
    sublead_hard_jet_dict.insert({key, TLorentzVector()});
    sublead_hard_jet_const_pt_dict.insert({key, new TH1D(dijetcore::MakeString(key, "jsconstpt").c_str(),
                                                         "", 90, 0, 30)});
    sublead_hard_jet_const_dr_dict.insert({key, new TH1D(dijetcore::MakeString(key, "jsconstdr").c_str(),
                                                         "", 100, 0.0, 1.0)});
    sublead_hard_jet_nconst_dict.insert({key, 0});
    sublead_hard_rho_dict.insert({key, 0});
    sublead_hard_sigma_dict.insert({key, 0});
    sublead_hard_area_dict.insert({key, 0});
    sublead_match_jet_dict.insert({key, TLorentzVector()});
    sublead_match_jet_const_pt_dict.insert({key, new TH1D(dijetcore::MakeString(key, "jsmconstpt").c_str(),
                                                          "", 90, 0, 30)});
    sublead_match_jet_const_dr_dict.insert({key, new TH1D(dijetcore::MakeString(key, "jsmconstdr").c_str(),
                                                          "", 100, 0.0, 1.0)});
    sublead_match_jet_nconst_dict.insert({key, 0});
    sublead_match_rho_dict.insert({key, 0});
    sublead_match_sigma_dict.insert({key, 0});
    sublead_match_area_dict.insert({key, 0});
    lead_off_axis_jet_dict.insert({key, TLorentzVector()});
    lead_off_axis_rho_dict.insert({key, 0});
    lead_off_axis_sigma_dict.insert({key, 0});
    sublead_off_axis_jet_dict.insert({key, TLorentzVector()});
    sublead_off_axis_rho_dict.insert({key, 0});
    sublead_off_axis_sigma_dict.insert({key, 0});
  }
  
  for (auto key : keys) {
    std::shared_ptr<TTree> tmp = std::make_shared<TTree>(key.c_str(), key.c_str());
    // create branches for the tree
    tmp->Branch("runid", &run_id_dict[key]);
    tmp->Branch("eventid", &event_id_dict[key]);
    tmp->Branch("vz", &vz_dict[key]);
    tmp->Branch("refmult", &refmult_dict[key]);
    tmp->Branch("grefmult", &grefmult_dict[key]);
    tmp->Branch("refmultcorr", &refmultcorr_dict[key]);
    tmp->Branch("grefmultcorr", &grefmultcorr_dict[key]);
    tmp->Branch("cent", &cent_dict[key]);
    tmp->Branch("zdcrate", &zdcrate_dict[key]);
    tmp->Branch("rp", &reactionplane_dict[key]);
    tmp->Branch("nglobal", &nglobal_dict[key]);
    tmp->Branch("npart", &npart_dict[key]);
    tmp->Branch("jl", &lead_hard_jet_dict[key]);
    tmp->Branch("js", &sublead_hard_jet_dict[key]);
    tmp->Branch("jlm", &lead_match_jet_dict[key]);
    tmp->Branch("jsm", &sublead_match_jet_dict[key]);
    tmp->Branch("jloa", &lead_off_axis_jet_dict[key]);
    tmp->Branch("jsoa", &sublead_off_axis_jet_dict[key]);
    tmp->Branch("jlconst", &lead_hard_jet_nconst_dict[key]);
    tmp->Branch("jlrho", &lead_hard_rho_dict[key]);
    tmp->Branch("jlsig", &lead_hard_sigma_dict[key]);
    tmp->Branch("jlarea", &lead_hard_area_dict[key]);
    tmp->Branch("jlmconst", &lead_match_jet_nconst_dict[key]);
    tmp->Branch("jlmrho", &lead_match_rho_dict[key]);
    tmp->Branch("jlmsig", &lead_match_sigma_dict[key]);
    tmp->Branch("jlmarea", &lead_match_area_dict[key]);
    tmp->Branch("jsconst", &sublead_hard_jet_nconst_dict[key]);
    tmp->Branch("jsrho", &sublead_hard_rho_dict[key]);
    tmp->Branch("jssig", &sublead_hard_sigma_dict[key]);
    tmp->Branch("jsarea", &sublead_hard_area_dict[key]);
    tmp->Branch("jsmconst", &sublead_match_jet_nconst_dict[key]);
    tmp->Branch("jsmrho", &sublead_match_rho_dict[key]);
    tmp->Branch("jsmsig", &sublead_match_sigma_dict[key]);
    tmp->Branch("jsmarea", &sublead_match_area_dict[key]);
    tmp->Branch("jloaconst", &lead_off_axis_nconst_dict[key]);
    tmp->Branch("jloarho", &lead_off_axis_rho_dict[key]);
    tmp->Branch("jloasig", &lead_off_axis_sigma_dict[key]);
    tmp->Branch("jsoaconst", &sublead_off_axis_nconst_dict[key]);
    tmp->Branch("jsoarho", &sublead_off_axis_rho_dict[key]);
    tmp->Branch("jsoasig", &sublead_off_axis_sigma_dict[key]);
    
    trees.insert({key, tmp});
  }
  
  // histograms
  // ----------
  
  std::unordered_map<string, TH1D*> lead_jet_count_dict;
  std::unordered_map<string, TH1D*> sublead_jet_count_dict;
  
  for (auto key : keys) {
    // create a unique histogram name for each key
    string lead_name = key + "_lead_count";
    string sublead_name = key + "_sublead_count";
    TH1D* lead_tmp = new TH1D(lead_name.c_str(), "count lead jets", 800, 0.5, 800.5);
    TH1D* sublead_tmp = new TH1D(sublead_name.c_str(), "count sublead jets", 800, 0.5, 800.5);
    
    lead_jet_count_dict.insert({key, lead_tmp});
    sublead_jet_count_dict.insert({key, sublead_tmp});
  }
  
  // define a selector to reject low momentum tracks
  fastjet::Selector track_pt_min_selector = fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  
  try {
    while (reader->NextEvent()) {
      // Print out reader status every 10 seconds
      reader->PrintStatus(10);
      
      // headers for convenience
      TStarJetPicoEventHeader* header = reader->GetEvent()->GetHeader();
      
      // check if event fired a trigger we will use
      if (triggers.size() != 0) {
        bool use_event = false;
        for (auto trigger : triggers)
          if (header->HasTriggerId(trigger))
            use_event = true;
        if (!use_event) continue;
      }
      
      // get centrality
      unsigned centrality_bin = -1;
      for (int i = 0; i < cent_boundaries.size(); ++i) {
        if (header->GetProperReferenceMultiplicity() > cent_boundaries[i]) {
          centrality_bin = i;
          break;
        }
      }
      
      // get the vector container
      TStarJetVectorContainer<TStarJetVector>* container = reader->GetOutputContainer();
      std::vector<fastjet::PseudoJet> primary_particles;
      dijetcore::ConvertTStarJetVector(container, primary_particles);
      
      // select tracks above the minimum pt threshold
      primary_particles = track_pt_min_selector(primary_particles);
      
      // run the worker
      auto& worker_out = worker.Run(primary_particles);
      
      std::unordered_map<std::string, dijetcore::unique_ptr<dijetcore::OffAxisOutput>> off_axis_worker_out;
      if (off_axis_worker)
        off_axis_worker_out = std::move(off_axis_worker->Run(worker, centrality_bin));
      
      // process any found di-jet pairs
      
      for (auto& result : worker_out) {
        std::string key = result.first;
        dijetcore::ClusterOutput& out = *result.second.get();
        
        if (out.found_lead)
          lead_jet_count_dict[key]->Fill(header->GetReferenceMultiplicity());
        if (out.found_sublead)
          sublead_jet_count_dict[key]->Fill(header->GetReferenceMultiplicity());
  
        // now fill dijet results
        if (out.found_match) {
          
          // fill all branches for that key
          run_id_dict[key] = header->GetRunId();
          event_id_dict[key] = header->GetEventId();
          vz_dict[key] = header->GetPrimaryVertexZ();
          refmult_dict[key] = header->GetReferenceMultiplicity();
          grefmult_dict[key] = header->GetGReferenceMultiplicity();
          refmultcorr_dict[key] = header->GetCorrectedReferenceMultiplicity();
          grefmultcorr_dict[key] = header->GetCorrectedGReferenceMultiplicity();
          cent_dict[key] = centrality_bin;
          zdcrate_dict[key] = header->GetZdcCoincidenceRate();
          reactionplane_dict[key] = header->GetReactionPlaneAngle();
          nglobal_dict[key] = header->GetNGlobalTracks();
          npart_dict[key] = primary_particles.size();
          
          // set the four jets
          lead_hard_jet_dict[key] = TLorentzVector(out.lead_hard.px(),
                                                   out.lead_hard.py(),
                                                   out.lead_hard.pz(),
                                                   out.lead_hard.E());
          int l_const = 0;
          for (int i = 0; i < out.lead_hard.constituents().size(); ++i) {
            if (out.lead_hard.constituents()[i].pt() > 0.01) {
              l_const++;
              lead_hard_jet_const_pt_dict[key]->Fill(out.lead_hard.constituents()[i].pt());
              lead_hard_jet_const_dr_dict[key]->Fill(out.lead_hard.constituents()[i].delta_R(out.lead_hard));
            }
          }
          lead_hard_jet_nconst_dict[key] = l_const;
          lead_hard_rho_dict[key] = out.lead_hard_rho;
          lead_hard_sigma_dict[key] = out.lead_hard_sigma;
          lead_hard_area_dict[key] = out.lead_hard.area();
          lead_match_jet_dict[key] = TLorentzVector(out.lead_match.px(),
                                                    out.lead_match.py(),
                                                    out.lead_match.pz(),
                                                    out.lead_match.E());
          int lm_const = 0;
          for (int i = 0; i < out.lead_match.constituents().size(); ++i) {
            if (out.lead_match.constituents()[i].pt() > 0.01) {
                lm_const++;
              lead_match_jet_const_pt_dict[key]->Fill(out.lead_match.constituents()[i].pt());
              lead_match_jet_const_dr_dict[key]->Fill(out.lead_match.constituents()[i].delta_R(out.lead_match));
            }
          }
          lead_match_jet_nconst_dict[key] = lm_const;
          lead_match_rho_dict[key] = out.lead_match_rho;
          lead_match_sigma_dict[key] = out.lead_match_sigma;
          lead_match_area_dict[key] = out.lead_match.area();
          sublead_hard_jet_dict[key] = TLorentzVector(out.sublead_hard.px(),
                                                      out.sublead_hard.py(),
                                                      out.sublead_hard.pz(),
                                                      out.sublead_hard.E());
          int s_const = 0;
          for (int i = 0; i < out.sublead_hard.constituents().size(); ++i) {
            if (out.sublead_hard.constituents()[i].pt() > 0.01) {
              s_const++;
              sublead_hard_jet_const_pt_dict[key]->Fill(out.sublead_hard.constituents()[i].pt());
              sublead_hard_jet_const_dr_dict[key]->Fill(out.sublead_hard.constituents()[i].delta_R(out.sublead_hard));
            }
          }
          sublead_hard_jet_nconst_dict[key] = s_const;
          sublead_hard_rho_dict[key] = out.sublead_hard_rho;
          sublead_hard_sigma_dict[key] = out.sublead_hard_sigma;
          sublead_hard_area_dict[key] = out.sublead_hard.area();
          sublead_match_jet_dict[key] = TLorentzVector(out.sublead_match.px(),
                                                       out.sublead_match.py(),
                                                       out.sublead_match.pz(),
                                                       out.sublead_match.E());
          int sm_const = 0;
          for (int i = 0; i < out.sublead_match.constituents().size(); ++i) {
            if (out.sublead_match.constituents()[i].pt() > 0.01) {
              sm_const++;
              sublead_match_jet_const_pt_dict[key]->Fill(out.sublead_match.constituents()[i].pt());
              sublead_match_jet_const_dr_dict[key]->Fill(out.sublead_match.constituents()[i].delta_R(out.sublead_match));
            }
          }
          sublead_match_jet_nconst_dict[key] = sm_const;
          sublead_match_rho_dict[key] = out.sublead_match_rho;
          sublead_match_sigma_dict[key] = out.sublead_match_sigma;
          sublead_match_area_dict[key] = out.sublead_match.area();
          
          if (off_axis_worker != nullptr &&
              off_axis_worker_out.find(key) != off_axis_worker_out.end() &&
              off_axis_worker_out[key]->lead_jet.has_constituents() &&
              off_axis_worker_out[key]->sub_jet.has_constituents()) {
            
            auto& off_axis_out = *off_axis_worker_out[key].get();
            lead_off_axis_jet_dict[key] = TLorentzVector(off_axis_out.lead_jet.px(),
                                                         off_axis_out.lead_jet.py(),
                                                         off_axis_out.lead_jet.pz(),
                                                         off_axis_out.lead_jet.E());
            int loa_const = 0;
            for (int i = 0; i < off_axis_out.lead_jet.constituents().size(); ++i)
              if (off_axis_out.lead_jet.constituents()[i].pt() > 0.01)
                loa_const++;
            lead_off_axis_nconst_dict[key] = loa_const;
            lead_off_axis_rho_dict[key] = off_axis_out.lead_jet_rho;
            lead_off_axis_sigma_dict[key] = off_axis_out.lead_jet_sigma;
            sublead_off_axis_jet_dict[key] = TLorentzVector(off_axis_out.sub_jet.px(),
                                                            off_axis_out.sub_jet.py(),
                                                            off_axis_out.sub_jet.pz(),
                                                            off_axis_out.sub_jet.E());
            int soa_const = 0;
            for (int i = 0; i < off_axis_out.sub_jet.constituents().size(); ++i) {
              if (off_axis_out.sub_jet.constituents()[i].pt() > 0.01) {
                soa_const++;
              }
            }
            sublead_off_axis_nconst_dict[key] = soa_const;
            sublead_off_axis_rho_dict[key] = off_axis_out.sub_jet_rho;
            sublead_off_axis_sigma_dict[key] = off_axis_out.sub_jet_sigma;
          }
          else {
            lead_off_axis_jet_dict[key] = TLorentzVector();
            lead_off_axis_nconst_dict[key] = 0;
            lead_off_axis_rho_dict[key] = 0;
            lead_off_axis_sigma_dict[key] = 0;
            sublead_off_axis_jet_dict[key] = TLorentzVector();
            sublead_off_axis_nconst_dict[key] = 0;
            sublead_off_axis_rho_dict[key] = 0;
            sublead_off_axis_sigma_dict[key] = 0;
          }
          
          trees[key]->Fill();
        }
      }
    }
  } catch(std::exception& e) {
    LOG(ERROR) << "Caught: " << e.what() << " during analysis loop.";
  }
  
  out.Write();
  out.Close();
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}
