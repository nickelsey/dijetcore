#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/data/centrality/centrality_run7.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/trigger_lookup.h"

#include "jetreader/reader/event_selector.h"
#include "jetreader/reader/reader.h"
#include "jetreader/reader/tower_selector.h"
#include "jetreader/reader/track_selector.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TTree.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "output_directory", // directory for output root file to written to
    "bad_tower_list",   // bad tower list for TStarJetPicoReader
    "bad_run_list",     // bad run list for TStarJetPicoReader
    "triggers",         // trigger id set (defined in
                        // dijetcore/util/data/trigger_lookup.h)
};

int main(int argc, char *argv[]) {
  string usage = "y14 qa routine";

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
  string outfile_name =
      output_dir + "/" + FLAGS_name + dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // initialize the reader
  jetreader::Reader reader(FLAGS_input);
  
  if (!std::string(config["bad_run_list"]).empty())
    reader.eventSelector()->addBadRuns(std::string(config["bad_run_list"]));
  if (!std::string(config["triggers"]).empty())
    reader.eventSelector()->addTriggerIds(config["triggers"]);
  reader.eventSelector()->setVzRange(-30, 30);
  reader.eventSelector()->setdVzMax(3.0);

  reader.trackSelector()->setDcaMax(3.0);
  reader.trackSelector()->setNHitsMin(15);
  reader.trackSelector()->setNHitsFracMin(0.52);
  reader.trackSelector()->setPtMax(30.0);
  reader.trackSelector()->setPtMin(0.2);
  reader.trackSelector()->rejectEventOnPtFailure();
  
  if (!std::string(config["bad_tower_list"]).empty())
    reader.towerSelector()->addBadTowers(std::string(config["bad_tower_list"]));
  reader.towerSelector()->setEtMax(30.0);
  reader.towerSelector()->setEtMin(0.2);
  reader.towerSelector()->rejectEventOnEtFailure();

  reader.init();

  // turn on centrality for low/mid luminosity
  reader.centrality().loadCentralityDef(jetreader::CentDefId::Run14LowMid);

  // output tree
  TTree *output = new TTree("ResultTree", "aj for auau");

  try {
    while (reader.next()) {

      // set event quantities for tree
      unsigned runid = reader.picoDst()->event()->runId();
      unsigned eventid = reader.picoDst()->event()->eventId();
      unsigned ref = reader.picoDst()->event()->refMult();
      unsigned gref = reader.picoDst()->event()->grefMult();
      double vz = reader.picoDst()->event()->primaryVertex().Z();
      unsigned cent = reader.centrality16();

      std::vector<fastjet::PseudoJet> primary_particles = reader.pseudojets();

      LOG(INFO) << "cent: " << cent;
    }
  } catch (std::exception &e) {
    LOG(ERROR) << "Caught: " << e.what() << " during analysis loop.";
  }

  output->Write();
  out.Close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}