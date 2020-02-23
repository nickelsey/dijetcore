// Aj of pythia in a detector simulation of the STAR detector
// for both generator and detector level jets

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"

#include <vector>
#include <algorithm>
#include <random>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include "RooUnfold/RooUnfold.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "output_directory",   // directory for output root file to written to
    "response_frac",      // fraction of events to use for response matrix
    "match_radius",       // matching radius for det/gen matching
    "lead_jet_pt",        // lead hard-core jet pt minimum det
    "lead_jet_pt_gen",    // lead hard-core jet pt minimum gen
    "sublead_jet_pt",     // sublead hard-core jet pt minimum det
    "sublead_jet_pt_gen", // sublead hard-core jet pt minimum gen
};

struct MatchEvent {
  double lead_hc_gen_pt;
  double sub_hc_gen_pt;
  double lead_match_gen_pt;
  double sub_match_gen_pt;
  double lead_hc_det_pt;
  double sub_hc_det_pt;
  double lead_match_det_pt;
  double sub_match_det_pt;
};

int main(int argc, char *argv[]) {

  string usage =
      "Unfolding of Pythia Aj with STAR detector sim";

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

  // load input file and setup TTreeReader
  TFile input(FLAGS_input.c_str(), "READ");
  TTreeReader input_reader("resultTree", &input);
  TTreeReaderValue<TLorentzVector> lead_hc_gen_r(input_reader, "leadhcgen");
  TTreeReaderValue<TLorentzVector> sub_hc_gen_r(input_reader, "subhcgen");
  TTreeReaderValue<TLorentzVector> lead_match_gen_r(input_reader,
                                                    "leadmatchgen");
  TTreeReaderValue<TLorentzVector> sub_match_gen_r(input_reader, "submatchgen");
  TTreeReaderValue<TLorentzVector> lead_hc_det_r(input_reader, "leadhcdet");
  TTreeReaderValue<TLorentzVector> sub_hc_det_r(input_reader, "subhcdet");
  TTreeReaderValue<TLorentzVector> lead_match_det_r(input_reader,
                                                    "leadmatchdet");
  TTreeReaderValue<TLorentzVector> sub_match_det_r(input_reader, "submatchdet");

  // get pt cuts and match radius
  double match_rad = config["match_radius"];
  double lead_gen_min_pt = config["lead_jet_pt_gen"];
  double lead_min_pt = config["lead_jet_pt"];
  double sub_gen_min_pt = config["sublead_jet_pt_gen"];
  double sub_min_pt = config["sublead_jet_pt"];
  std::vector<MatchEvent> matched_events;
  matched_events.reserve(10000);
  while (input_reader.Next()) {
    fastjet::PseudoJet lead_hc_gen(*lead_hc_gen_r);
    fastjet::PseudoJet sub_hc_gen(*sub_hc_gen_r);
    fastjet::PseudoJet lead_match_gen(*lead_match_gen_r);
    fastjet::PseudoJet sub_match_gen(*sub_match_gen_r);
    fastjet::PseudoJet lead_hc_det(*lead_hc_det_r);
    fastjet::PseudoJet sub_hc_det(*sub_hc_det_r);
    fastjet::PseudoJet lead_match_det(*lead_match_det_r);
    fastjet::PseudoJet sub_match_det(*sub_match_det_r);

    if (lead_hc_gen.pt() < lead_gen_min_pt ||
        sub_hc_gen.pt() < sub_gen_min_pt || lead_hc_det.pt() < lead_min_pt ||
        sub_hc_det.pt() < sub_min_pt)
      continue;

    bool match_lead = lead_hc_gen.delta_R(lead_hc_det) < match_rad ||
                      lead_hc_gen.delta_R(sub_hc_det) < match_rad;
    bool match_sub = sub_hc_gen.delta_R(lead_hc_det) < match_rad ||
                     sub_hc_gen.delta_R(sub_hc_det) < match_rad;
    if (!match_lead || !match_sub)
      continue;

    MatchEvent tmp;
    tmp.lead_hc_gen_pt = lead_hc_gen.pt();
    tmp.lead_match_gen_pt = lead_match_gen.pt();
    tmp.sub_hc_gen_pt = sub_hc_gen.pt();
    tmp.sub_match_gen_pt = sub_match_gen.pt();
    tmp.lead_hc_det_pt = lead_hc_det.pt();
    tmp.lead_match_det_pt = lead_match_det.pt();
    tmp.sub_hc_det_pt = sub_hc_det.pt();
    tmp.sub_match_det_pt = sub_match_det.pt();
    matched_events.push_back(tmp);
  }

  LOG(INFO) << "accepted events" << matched_events.size();

  // shuffle the events just in case
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(matched_events.begin(), matched_events.end(), g);

  // RooUnfoldResponse *response =
  //    new RooUnfoldResponse("aj3dresponse", "3D Aj response for Pythia");

  gflags::ShutDownCommandLineFlags();
  return 0;
}