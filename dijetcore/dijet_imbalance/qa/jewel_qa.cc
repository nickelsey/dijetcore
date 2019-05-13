#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <unordered_map>

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/mc/jewel_reader.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TTree.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

// DIJETCORE_DEFINE_string(input, "", "input file");
// DIJETCORE_DEFINE_string(outputDir, "tmp", "directory for output");
// DIJETCORE_DEFINE_string(name, "job", "name for output file");
// DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
// DIJETCORE_DEFINE_int(nEvents, -1, "number of events (-1 for all)");

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "output_directory", // directory for output root file to written to
    "n_events",          // number of events to run over
};

int main(int argc, char *argv[]) {

  string usage = "JEWEL QA for trigger surface bias";

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

  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(FLAGS_input)) {
    LOG(ERROR) << "input file does not exist: " << FLAGS_input;
    ;
    return 1;
  }

  // attempt to initialize the jewel reader
  dijetcore::JewelReader reader(FLAGS_input);

  // build output directory if it doesn't exist, using boost::filesystem
  string output_dir = config["output_directory"];
  if (output_dir.empty())
    output_dir = "tmp";
  boost::filesystem::path dir(output_dir.c_str());
  boost::filesystem::create_directories(dir);

  // create output file from the given directory, name & id
  string outfile_name = output_dir + "/" + FLAGS_name +
                        dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // define the clustering jet definition
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);

  // define a selector to reject low momentum tracks
  fastjet::Selector track_pt_min_selector =
      fastjet::SelectorPtMin(2.0) && fastjet::SelectorAbsRapMax(1.0);
  fastjet::Selector jet_sel =
      fastjet::SelectorPtMin(5.0) && fastjet::SelectorAbsRapMax(0.6);

  LOG(INFO) << "building tree";
  TTree* result_tree = new TTree("T", "tree");
  TLorentzVector trigger_jet;
  TLorentzVector trigger_particle;
  TLorentzVector p1, p2;
  double w, xsec;
  double vx, vy;

  result_tree->Branch("tj", &trigger_jet);
  result_tree->Branch("tp", &trigger_particle);
  result_tree->Branch("p1", &p1);
  result_tree->Branch("p2", &p2);
  result_tree->Branch("w", &w);
  result_tree->Branch("xsec", &xsec);
  result_tree->Branch("vx", &vx);
  result_tree->Branch("vy", &vy);

  int event = 0;
  size_t max_event = 0;
  if (config["n_events"] < 0)
    max_event = std::numeric_limits<size_t>::max();
  else
    max_event = config["n_events"];
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

      // cluster
      fastjet::ClusterSequence cs(primary_particles, jet_def);

      std::vector<fastjet::PseudoJet> results =
          fastjet::sorted_by_pt(jet_sel(cs.inclusive_jets()));

      if (results.size() == 0)
        continue;

      trigger_jet = TLorentzVector(results[0].px(), results[0].py(),
                                   results[0].pz(), results[0].E());

      // find trigger
      for (auto &p : primary_particles) {
        if (p.pt() > trigger_particle.Pt())
          trigger_particle = TLorentzVector(p.px(), p.py(), p.pz(), p.E());
      }

      // fill all branches for that key
      vx = reader.vx();
      vy = reader.vy();
      w = reader.weight();
      xsec = reader.totalXsec();

      fastjet::PseudoJet lead_p = reader.leadingParton();
      fastjet::PseudoJet sub_p = reader.subParton();

      p1 = TLorentzVector(lead_p.px(), lead_p.py(), lead_p.pz(), lead_p.E());
      p2 = TLorentzVector(sub_p.px(), sub_p.py(), sub_p.pz(), sub_p.E());

      result_tree->Fill();
    }
  } catch (std::exception &e) {
    LOG(ERROR) << "Caught: " << e.what() << " during analysis loop.";
  }

  out.cd();
  result_tree->Write();
  out.Close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}
