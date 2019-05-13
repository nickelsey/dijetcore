#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <random>
#include <exception>
#include <limits>

#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/worker/dijet_worker/dijet_worker.h"
#include "dijetcore/util/mc/jewel_reader.h"

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"


using std::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(outputDir, "tmp", "directory for output");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
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
DIJETCORE_DEFINE_int(nEvents, 10000, "number of events (-1 for all)");

int main(int argc, char* argv[]) {
  
  string usage = "JEWEL differential di-jet imbalance analysis routine";
  
  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);
  
  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(FLAGS_input)) {
    LOG(ERROR) << "input file does not exist: " << FLAGS_input;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);
  
  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outputDir + "/" + FLAGS_name + dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");
  
  // attempt to initialize the jewel reader
  dijetcore::JewelReader reader(FLAGS_input);
  
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
  worker.forceConstituentPtEquality(FLAGS_forceConstituentPtEquality);
  worker.forceConstituentEtaEquality(FLAGS_forceConstituentEtaEquality);
  worker.forceJetResolutionEquality(FLAGS_forceJetResolutionEquality);
  worker.forceMatchJetResolutionEquality(FLAGS_forceMatchJetResolutionEquality);
  worker.initialize();
  
  std::set<std::string> keys = worker.keys();
  
  for (auto key : keys)
    LOG(INFO) << key;
  
  // create an output tree for each definition
  // -----------------------------------------
  
  std::unordered_map<std::string, std::shared_ptr<TTree>> trees;
  
  // and the necessary branches
  std::unordered_map<std::string, double> vx_dict;
  std::unordered_map<std::string, double> vy_dict;
  std::unordered_map<std::string, double> w_dict;
  std::unordered_map<std::string, double> xsec_dict;
  std::unordered_map<std::string, double> trig_dict;

  // initiating partons
  std::unordered_map<std::string, TLorentzVector> lead_parton_dict;
  std::unordered_map<std::string, TLorentzVector> sub_parton_dict;
  std::unordered_map<std::string, int> lead_parton_id_dict;
  std::unordered_map<std::string, int> sub_parton_id_dict;

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

  // fill the maps first, so that they don't decide to resize/move themselves
  // after branch creation...
  for (auto key : keys) {
    vx_dict.insert({key, 0.0});
    vy_dict.insert({key, 0.0});
    w_dict.insert({key, 0.0});
    xsec_dict.insert({key, 0.0});
    trig_dict.insert({key, 0.0});

    lead_parton_dict.insert({key, TLorentzVector()});
    sub_parton_dict.insert({key, TLorentzVector()});
    lead_parton_id_dict.insert({key, 0});
    sub_parton_id_dict.insert({key, 0});

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
    
    
  }
  
  for (auto key : keys) {
    dijetcore::unique_ptr<TTree> tmp = dijetcore::make_unique<TTree>(key.c_str(), key.c_str());
    // create branches for the tree
    tmp->Branch("vx", &vx_dict[key]);
    tmp->Branch("vy", &vy_dict[key]);
    tmp->Branch("w", &w_dict[key]);
    tmp->Branch("xsec", &xsec_dict[key]);
    tmp->Branch("e", &trig_dict[key]);

    tmp->Branch("pl", &lead_parton_dict[key]);
    tmp->Branch("ps", &sub_parton_dict[key]);
    tmp->Branch("plid", &lead_parton_id_dict[key]);
    tmp->Branch("psid", &sub_parton_id_dict[key]);

    tmp->Branch("jl", &lead_hard_jet_dict[key]);
    tmp->Branch("js", &sublead_hard_jet_dict[key]);
    tmp->Branch("jlm", &lead_match_jet_dict[key]);
    tmp->Branch("jsm", &sublead_match_jet_dict[key]);
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
  
  int event = 0;
  size_t max_event = 0;
  if (FLAGS_nEvents < 0)
    max_event = std::numeric_limits<size_t>::max();
  else 
    max_event = FLAGS_nEvents;
  try {
    while (reader.next()) {
      if (event >= max_event)
        break;
      if (event % 500 == 0) {
        LOG(INFO) << "Event: " << event;
      }
      event++;

      std::vector<fastjet::PseudoJet> primary_particles = reader.event();
      
      // select tracks above the minimum pt threshold
      primary_particles = track_pt_min_selector(primary_particles);  
      
      // run the worker
      auto& worker_out = worker.run(primary_particles);

      // process any found di-jet pairs

      for (auto& result : worker_out) {
        std::string key = result.first;
        dijetcore::ClusterOutput& out = *result.second.get();

        // now fill dijet results
        if (out.found_match) {

          // fill all branches for that key
          vx_dict[key] = reader.vx();
          vy_dict[key] = reader.vy();
          w_dict[key] = reader.weight();
          xsec_dict[key] = reader.totalXsec();

          fastjet::PseudoJet lead_p = reader.leadingParton();
          fastjet::PseudoJet sub_p = reader.subParton();

          lead_parton_dict[key] = TLorentzVector(lead_p.px(), lead_p.py(), lead_p.pz(), lead_p.E());
          sub_parton_dict[key] = TLorentzVector(sub_p.px(), sub_p.py(), sub_p.pz(), sub_p.E());
          lead_parton_id_dict[key] = lead_p.user_index();
          sub_parton_id_dict[key] = sub_p.user_index();

          // find trigger
          double trig_pt = 0.0;
          for (auto& p : primary_particles) {
            if (p.pt() > trig_pt)
              trig_pt = p.pt();
          }
          trig_dict[key] = trig_pt;

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

          trees[key]->Fill();
        }
      }
    }
  } catch(std::exception& e) {
    LOG(ERROR) << "Caught: " << e.what() << " during analysis loop.";
  }
  
  out.cd();
  for (auto& entry : trees) {
    entry.second->Write();
  }
  out.Close();
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}
