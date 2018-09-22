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
#include "dijetcore/util/data/centrality/centrality_run7.h"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile2D.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"


using std::string;

DIJETCORE_DEFINE_string(outputDir, "tmp", "directory for output");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_int(nEvents, 1000000, "number of pythia events to generate");
DIJETCORE_DEFINE_double(pTHatMin, 10.0, "minimum pT for hard interaction");
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
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);
  
  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outputDir + "/" + FLAGS_name + dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // create the generator
  // --------------------
  std::stringstream pt_stream;
  pt_stream << FLAGS_pTHatMin;
  string pt_low_string = "PhaseSpace:pTHatMin = " + pt_stream.str();
  string pt_high_string = "PhaseSpace:pTHatMax = 100.0";
  string qcd_string = "HardQCD:all = on";
  string e_cm_string = "Beams:eCM = 200";
  Pythia8::Pythia pythia;
  pythia.readString(qcd_string);
  pythia.readString(e_cm_string);
  pythia.readString(pt_low_string);
  pythia.readString(pt_high_string);
  pythia.readString("Random:setSeed = On");
  pythia.readString("Random:seed = 0");
  pythia.readString("HadronLevel:Decay = on");
  pythia.init();
  
  // initialize Run 7 centrality
  dijetcore::CentralityRun7 centrality;
  
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
  std::unordered_map<std::string, int> off_axis_centrality;
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
    
    off_axis_centrality.insert({key, 0});
    lead_off_axis_jet_dict.insert({key, TLorentzVector()});
    lead_off_axis_rho_dict.insert({key, 0});
    lead_off_axis_sigma_dict.insert({key, 0});
    sublead_off_axis_jet_dict.insert({key, TLorentzVector()});
    sublead_off_axis_rho_dict.insert({key, 0});
    sublead_off_axis_sigma_dict.insert({key, 0});
  }
  
  for (auto key : keys) {
    dijetcore::unique_ptr<TTree> tmp = dijetcore::make_unique<TTree>(key.c_str(), key.c_str());
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
    tmp->Branch("oacent", &off_axis_centrality[key]);
    tmp->Branch("jloaconst", &lead_off_axis_nconst_dict[key]);
    tmp->Branch("jloarho", &lead_off_axis_rho_dict[key]);
    tmp->Branch("jloasig", &lead_off_axis_sigma_dict[key]);
    tmp->Branch("jsoaconst", &sublead_off_axis_nconst_dict[key]);
    tmp->Branch("jsoarho", &sublead_off_axis_rho_dict[key]);
    tmp->Branch("jsoasig", &sublead_off_axis_sigma_dict[key]);
    
    trees[key] = std::move(tmp);
    trees[key]->SetDirectory(0);
  }
  
  // histograms
  // ----------
  
  std::unordered_map<string, dijetcore::unique_ptr<TH1D>> lead_jet_count_dict;
  std::unordered_map<string, dijetcore::unique_ptr<TH1D>> sublead_jet_count_dict;
  
  for (auto key : keys) {
    // create a unique histogram name for each key
    string lead_name = key + "_lead_count";
    string sublead_name = key + "_sublead_count";
    dijetcore::unique_ptr<TH1D> lead_tmp = dijetcore::make_unique<TH1D>(lead_name.c_str(), "count lead jets", 800, 0.5, 800.5);
    dijetcore::unique_ptr<TH1D> sublead_tmp = dijetcore::make_unique<TH1D>(sublead_name.c_str(), "count sublead jets", 800, 0.5, 800.5);
    
    lead_jet_count_dict[key] = std::move(lead_tmp);
    sublead_jet_count_dict[key] = std::move(sublead_tmp);
    lead_jet_count_dict[key]->SetDirectory(0);
    sublead_jet_count_dict[key]->SetDirectory(0);
  }
  
  // define a selector to reject low momentum tracks
  fastjet::Selector track_pt_min_selector = fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  
  // event loop
  for (int i = 0; i < FLAGS_nEvents; ++i) {
    
    // generate next event
    if (!pythia.next())
      continue;
    
    std::vector<fastjet::PseudoJet> particles;
    
    for (int j = 0; j < pythia.event.size(); ++j) {
      if (pythia.event[j].isFinal() && pythia.event[j].isVisible()) {
        fastjet::PseudoJet tmp(pythia.event[j].px(), pythia.event[j].py(), pythia.event[j].pz(), pythia.event[j].e());
        tmp.set_user_index(pythia.event[j].charge());
        particles.push_back(tmp);
      }
    }
    
    auto& worker_out = worker.Run(particles);
    
    // process any found di-jet pairs
    for (auto& result : worker_out) {
      std::string key = result.first;
      dijetcore::ClusterOutput& out = *result.second.get();
      
      // now fill dijet results
      if (out.found_match) {
        // fill all branches for that key
        npart_dict[key] = particles.size();
        
        // set the four jets
        lead_hard_jet_dict[key] = TLorentzVector(out.lead_hard.px(),
                                                 out.lead_hard.py(),
                                                 out.lead_hard.pz(),
                                                 out.lead_hard.E());
        lead_hard_jet_nconst_dict[key] = out.lead_hard.constituents().size();
        lead_hard_rho_dict[key] = out.lead_hard_rho;
        lead_hard_sigma_dict[key] = out.lead_hard_sigma;
        lead_hard_area_dict[key] = out.lead_hard.area();
        lead_match_jet_dict[key] = TLorentzVector(out.lead_match.px(),
                                                  out.lead_match.py(),
                                                  out.lead_match.pz(),
                                                  out.lead_match.E());
        lead_match_jet_nconst_dict[key] = out.lead_match.constituents().size();
        lead_match_rho_dict[key] = out.lead_match_rho;
        lead_match_sigma_dict[key] = out.lead_match_sigma;
        lead_match_area_dict[key] = out.lead_match.area();
        sublead_hard_jet_dict[key] = TLorentzVector(out.sublead_hard.px(),
                                                    out.sublead_hard.py(),
                                                    out.sublead_hard.pz(),
                                                    out.sublead_hard.E());
        sublead_hard_jet_nconst_dict[key] = out.sublead_hard.constituents().size();
        sublead_hard_rho_dict[key] = out.sublead_hard_rho;
        sublead_hard_sigma_dict[key] = out.sublead_hard_sigma;
        sublead_hard_area_dict[key] = out.sublead_hard.area();
        sublead_match_jet_dict[key] = TLorentzVector(out.sublead_match.px(),
                                                     out.sublead_match.py(),
                                                     out.sublead_match.pz(),
                                                     out.sublead_match.E());
        sublead_match_jet_nconst_dict[key] = out.sublead_match.constituents().size();
        sublead_match_rho_dict[key] = out.sublead_match_rho;
        sublead_match_sigma_dict[key] = out.sublead_match_sigma;
        sublead_match_area_dict[key] = out.sublead_match.area();
        
        trees[key]->Fill();
      }
    }
  }
  

  
  out.cd();
  for (auto& entry : trees) {
    entry.second->Write();
  }
  out.Close();
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}

