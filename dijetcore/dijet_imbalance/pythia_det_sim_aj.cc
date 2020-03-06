// Aj of pythia in a detector simulation of the STAR detector
// for both generator and detector level jets

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/root/compile_vec_tlorentzvector.h"
#include "dijetcore/worker/dijet_worker/dijet_worker.h"
#include "dijetcore/worker/mc/detector_sim_worker/pythia_star_sim.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TROOT.h"
#include "TTree.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "output_directory", // directory for output root file to written to
    "n_events",         // number of events to generate
    "pthat_min",        // minimum pT hat for generator
    "pthat_max",        // maximum pT hat for generator
    "track_pt_min",     // minimum pT for tracks at detector level
    "track_pt_max",     // maximum pT for tracks at detector level
    "track_eta_max",    // maximum eta for tracks at detector level
    "det_sim_mode",     // detector simulation method (none, gauss, star)

    "const_eta",           // constituent  eta
    "const_eta_gen",       // constituent  eta
    "lead_const_pt",       // lead jet constituent pT cut (comma separated)
    "sublead_const_pt",    // sublead jet constituent pT cut (comma separated)
    "lead_const_pt_match", // same, for leading matched jets
    "lead_const_pt_match_gen",    // for leading matched jets generator level
    "sublead_const_pt_match",     // same, for subleading matched jets
    "sublead_const_pt_match_gen", // for subleading matched jets generator level
    "lead_r",                     // leading jet R (comma separated)
    "sublead_r",                  // leading jet R (comma separated)
    "lead_r_match",               // same for matched leading jet
    "sublead_r_match",            // same for matched subleading jet
    "lead_jet_pt",        // lead hard-core jet pt minimum (comma separated)
    "lead_jet_pt_gen",    // lead hard-core jet pt minimum (comma separated)
    "sublead_jet_pt",     // lead hard-core jet pt minimum (comma separated)
    "sublead_jet_pt_gen", // lead hard-core jet pt minimum (comma separated)
};

int main(int argc, char *argv[]) {

  string usage = "Aj for Pythia with detector sim of STAR detector";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // parse configuration file
  dijetcore::json config;
  try {
    config = dijetcore::LoadConfig(FLAGS_config, required_params);
  } catch (std::exception &e) {
    LOG(ERROR) << "error loading config: " << e.what();
    return 1;
  }

  // find output directory from configuration file
  string output_dir = config["output_directory"];

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(output_dir);
  boost::filesystem::create_directories(dir);

  // copy config file to output directory
  boost::filesystem::path input_file(FLAGS_config.c_str());
  boost::filesystem::path copy_path(dir);
  copy_path /= input_file.filename();
  boost::filesystem::copy_file(
      input_file, copy_path,
      boost::filesystem::copy_option::overwrite_if_exists);

  // create output file from the given directory, name & id
  string outfile_name = output_dir + "/" + FLAGS_name +
                        dijetcore::MakeString("_", FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // create generator
  dijetcore::PythiaStarSim generator;
  generator.DoClustering(false);
  generator.SetRandomSeed(FLAGS_id);
  generator.SetPtHatRange(config["pthat_min"].get<double>(),
                          config["pthat_max"].get<double>());
  generator.SetTrackPtMin(config["track_pt_min"].get<double>());
  generator.SetTrackPtMax(config["track_pt_max"].get<double>());
  generator.SetTrackEtaMax(config["track_eta_max"].get<double>());

  string opt_tag = "opt";
  int opt_idx = 1;
  while (config.count(dijetcore::MakeString(opt_tag, opt_idx)) > 0) {
    generator.AddSettingString(config[dijetcore::MakeString(opt_tag, opt_idx)]);
    opt_idx++;
  }

  string sim_mode = config["det_sim_mode"];

  if (sim_mode == "none") {
    generator.SetEfficiencyMode(dijetcore::PythiaStarSim::EfficiencyMode::None);
  } else if (sim_mode == "gauss") {
    generator.SetEfficiencyMode(
        dijetcore::PythiaStarSim::EfficiencyMode::GaussianSmearing);
  } else if (sim_mode == "star") {
    generator.SetEfficiencyMode(dijetcore::PythiaStarSim::EfficiencyMode::STAR);
  } else {
    LOG(ERROR) << "unrecognized simulation mode: " << sim_mode
               << " select from none, gauss or star";
  }

  // initialize the generator
  generator.Initialize();

  // parse jetfinding variables
  // --------------------------

  // first, hard code the algorithm to be anti-kt
  std::set<fastjet::JetAlgorithm> alg{fastjet::antikt_algorithm};

  // constituent range
  std::set<double> const_eta =
      dijetcore::ParseArgString<double>(config["const_eta"]);
  std::set<double> const_eta_gen =
      dijetcore::ParseArgString<double>(config["const_eta_gen"]);

  // leading jet
  std::set<double> lead_const_hard_pt =
      dijetcore::ParseArgString<double>(config["lead_const_pt"]);
  std::set<double> lead_const_match_pt =
      dijetcore::ParseArgString<double>(config["lead_const_pt_match"]);
  std::set<double> lead_const_match_pt_gen =
      dijetcore::ParseArgString<double>(config["lead_const_pt_match_gen"]);
  std::set<double> lead_R = dijetcore::ParseArgString<double>(config["lead_r"]);
  std::set<double> lead_R_match =
      dijetcore::ParseArgString<double>(config["lead_r_match"]);
  std::set<double> lead_hard_pt =
      dijetcore::ParseArgString<double>(config["lead_jet_pt"]);
  std::set<double> lead_hard_pt_gen =
      dijetcore::ParseArgString<double>(config["lead_jet_pt_gen"]);

  // subleading jet
  std::set<double> sublead_const_hard_pt =
      dijetcore::ParseArgString<double>(config["sublead_const_pt"]);
  std::set<double> sublead_const_match_pt =
      dijetcore::ParseArgString<double>(config["sublead_const_pt_match"]);
  std::set<double> sublead_const_match_pt_gen =
      dijetcore::ParseArgString<double>(config["sublead_const_pt_match_gen"]);
  std::set<double> sublead_R =
      dijetcore::ParseArgString<double>(config["sublead_r"]);
  std::set<double> sublead_R_match =
      dijetcore::ParseArgString<double>(config["sublead_r_match"]);
  std::set<double> sublead_hard_pt =
      dijetcore::ParseArgString<double>(config["sublead_jet_pt"]);
  std::set<double> sublead_hard_pt_gen =
      dijetcore::ParseArgString<double>(config["sublead_jet_pt_gen"]);

  LOG(INFO) << "grid variables ";
  LOG(INFO) << "constituent eta: " << const_eta;
  LOG(INFO) << "lead constituent hard pt cut: " << lead_const_hard_pt;
  LOG(INFO) << "lead constituent match pt cut: " << lead_const_match_pt;
  LOG(INFO) << "lead hard R: " << lead_R;
  LOG(INFO) << "lead match R: " << lead_R_match;
  LOG(INFO) << "lead hard jet pt: " << lead_hard_pt;
  LOG(INFO) << "sublead constituent hard pt cut: " << sublead_const_hard_pt;
  LOG(INFO) << "sublead constituent match pt cut: " << sublead_const_match_pt;
  LOG(INFO) << "sublead hard R: " << sublead_R;
  LOG(INFO) << "sublead match R: " << sublead_R_match;
  LOG(INFO) << "sublead hard jet pt: " << sublead_hard_pt;

  // here we can initialize the worker
  LOG(INFO) << "initializing worker...";
  dijetcore::DijetWorker worker_gen(
      alg, lead_hard_pt_gen, lead_R, lead_R_match, sublead_hard_pt_gen,
      sublead_R, sublead_R_match, lead_const_hard_pt, lead_const_match_pt_gen,
      sublead_const_hard_pt, sublead_const_match_pt_gen, const_eta_gen);
  worker_gen.forceConstituentPtEquality(true);
  worker_gen.forceConstituentEtaEquality(true);
  worker_gen.forceJetResolutionEquality(true);
  worker_gen.forceMatchJetResolutionEquality(true);
  worker_gen.initialize();

  dijetcore::DijetWorker worker_det(
      alg, lead_hard_pt, lead_R, lead_R_match, sublead_hard_pt, sublead_R,
      sublead_R_match, lead_const_hard_pt, lead_const_match_pt,
      sublead_const_hard_pt, sublead_const_match_pt, const_eta);
  worker_det.forceConstituentPtEquality(true);
  worker_det.forceConstituentEtaEquality(true);
  worker_det.forceJetResolutionEquality(true);
  worker_det.forceMatchJetResolutionEquality(true);
  worker_det.initialize();

  LOG(INFO) << "worker initialized - number of dijet definitions: "
            << worker_gen.size();

  std::set<std::string> keys = worker_gen.keys();
  for (auto key : keys)
    LOG(INFO) << key;

  if (worker_gen.keys().size() != 1 || worker_det.keys().size() != 1) {
    LOG(INFO) << "we do not support more than 1 key at a time currently";
    return 1;
  }

  // create an output tree
  // ---------------------
  std::string tree_title = dijetcore::MakeString(
      "leadpt_", config["lead_jet_pt"], "_subpt_", config["sublead_jet_pt"],
      "_leadr_", config["lead_r"], "_subleadr_", config["sublead_r"],
      "_leadconstpt_", config["lead_const_pt"], "_subleadconstpt_",
      config["sublead_const_pt"]);

  TTree *result_tree = new TTree("resultTree", tree_title.c_str());

  unsigned event_id;
  TLorentzVector lead_hc_gen;
  TLorentzVector lead_match_gen;
  TLorentzVector sub_hc_gen;
  TLorentzVector sub_match_gen;
  TLorentzVector lead_hc_det;
  TLorentzVector lead_match_det;
  TLorentzVector sub_hc_det;
  TLorentzVector sub_match_det;
  TLorentzVector trigger_gen;
  TLorentzVector trigger_det;
  TLorentzVector gen1;
  TLorentzVector gen2;
  bool gen1_gluon;
  bool gen2_gluon;

  result_tree->Branch("eventid", &event_id);
  result_tree->Branch("leadhcgen", &lead_hc_gen);
  result_tree->Branch("leadmatchgen", &lead_match_gen);
  result_tree->Branch("subhcgen", &sub_hc_gen);
  result_tree->Branch("submatchgen", &sub_match_gen);
  result_tree->Branch("leadhcdet", &lead_hc_det);
  result_tree->Branch("leadmatchdet", &lead_match_det);
  result_tree->Branch("subhcdet", &sub_hc_det);
  result_tree->Branch("submatchdet", &sub_match_det);
  result_tree->Branch("triggen", &trigger_gen);
  result_tree->Branch("trigdet", &trigger_det);
  result_tree->Branch("gen1", &gen1);
  result_tree->Branch("gen2", &gen2);
  result_tree->Branch("gen1_g", &gen1_gluon);
  result_tree->Branch("gen2_g", &gen2_gluon);

  // start event loop
  int n_events = config["n_events"].get<int>();
  for (int i = 0; i < n_events; ++i) {
    lead_hc_det = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    lead_match_det = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    sub_hc_det = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    sub_match_det = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    lead_hc_gen = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    lead_match_gen = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    sub_hc_gen = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    sub_match_gen = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    trigger_gen = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    trigger_det = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    gen1 = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    gen2 = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    gen1_gluon = false;
    gen2_gluon = false;

    generator.Run();

    std::vector<fastjet::PseudoJet> gen_part = generator.GenParticles();
    std::vector<fastjet::PseudoJet> det_part = generator.DetectorParticles();

    // run the workers
    auto &worker_out_gen = worker_gen.run(gen_part);
    auto &worker_out_det = worker_det.run(det_part);

    // fill tree entries
    event_id = i;

    // process any found di-jet pairs
    // there is only one key, so we can loop over the keys without
    // worrying about mixing
    for (auto &result : worker_out_gen) {
      std::string key = result.first;
      dijetcore::ClusterOutput &out = *result.second.get();

      if (out.found_match) {
        lead_hc_gen = TLorentzVector(out.lead_hard.px(), out.lead_hard.py(),
                                     out.lead_hard.pz(), out.lead_hard.E());
        lead_match_gen =
            TLorentzVector(out.lead_match.px(), out.lead_match.py(),
                           out.lead_match.pz(), out.lead_match.E());
        sub_hc_gen =
            TLorentzVector(out.sublead_hard.px(), out.sublead_hard.py(),
                           out.sublead_hard.pz(), out.sublead_hard.E());
        sub_match_gen =
            TLorentzVector(out.sublead_match.px(), out.sublead_match.py(),
                           out.sublead_match.pz(), out.sublead_match.E());
        // find highest pt neutral object for generator level
        auto possible_triggers = fastjet::sorted_by_pt(gen_part);
        for (auto &p : possible_triggers) {
          if (p.user_index() == 0) {
            trigger_gen = TLorentzVector(p.px(), p.py(), p.pz(), p.E());
            break;
          }
        }
      }
    }
    for (auto &result : worker_out_det) {
      std::string key = result.first;
      dijetcore::ClusterOutput &out = *result.second.get();

      if (out.found_match) {
        lead_hc_det = TLorentzVector(out.lead_hard.px(), out.lead_hard.py(),
                                     out.lead_hard.pz(), out.lead_hard.E());
        lead_match_det =
            TLorentzVector(out.lead_match.px(), out.lead_match.py(),
                           out.lead_match.pz(), out.lead_match.E());
        sub_hc_det =
            TLorentzVector(out.sublead_hard.px(), out.sublead_hard.py(),
                           out.sublead_hard.pz(), out.sublead_hard.E());
        sub_match_det =
            TLorentzVector(out.sublead_match.px(), out.sublead_match.py(),
                           out.sublead_match.pz(), out.sublead_match.E());

        // identify jet initiator (quark/gluon)
        gen1 = TLorentzVector(
            generator.Pythia().event[5].px(), generator.Pythia().event[5].py(),
            generator.Pythia().event[5].pz(), generator.Pythia().event[5].e());
        gen2 = TLorentzVector(
            generator.Pythia().event[6].px(), generator.Pythia().event[6].py(),
            generator.Pythia().event[6].pz(), generator.Pythia().event[6].e());

        gen1_gluon = generator.Pythia().event[5].id() == 21;
        gen2_gluon = generator.Pythia().event[6].id() == 21;

        // find highest pt neutral object for detector level
        auto possible_triggers = fastjet::sorted_by_pt(det_part);
        for (auto &p : possible_triggers) {
          if (p.user_index() == 0) {
            trigger_det = TLorentzVector(p.px(), p.py(), p.pz(), p.E());
            break;
          }
        }
      }
    }
    if (lead_hc_det.Pt() > 0.0 || lead_hc_gen.Pt() > 0.0)
      result_tree->Fill(); 
  }

  out.Write();
  out.Close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}