// check if the tower matching algorithm in our trees
// and the StPicoDst algorithm give similar results
// input must be from the same data

#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/data/reader_util.h"

#include <set>
#include <string>

#include "StPicoEvent/StPicoDstReader.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoReader.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetPicoTrackCuts.h"

DIJETCORE_DEFINE_string(stPicoFiles, "", "file or file list of StPicoDsts");
DIJETCORE_DEFINE_string(jetPicoFiles, "", "file or file list of TStarJetPicos");
DIJETCORE_DEFINE_string(outputDir, "tmp", "output directory");
DIJETCORE_DEFINE_string(
    name, "compare_pico",
    "output filename and identifier for batch job monitoring");
DIJETCORE_DEFINE_int(id, -1,
                     "job id for batch jobs and output naming. Leave negative "
                     "to ignore (when running single jobs");
DIJETCORE_DEFINE_int(nEvents, 50, "number of events to compare");
DIJETCORE_DEFINE_string(towList, "resources/bad_tower_lists/empty_list.txt",
                        "bad tower list");

using std::pair;
using std::string;
using dijetcore::dijetcore_map;
using dijetcore::PairHash;

using key_type = pair<unsigned, unsigned>;
using val_type = size_t;
using event_map = dijetcore_map<key_type, val_type, PairHash>;

// enumerate events for each dataset, returns a map where each key is 
// an event {runid, eventid} pair, and the value is the entry position in the TChain
event_map EnumerateEvents(TChain* chain) {
  event_map map;
  TTreeReader reader(chain);

  string runid_name = "fEventHeader.fRunId";
  string eventid_name = "fEventHeader.fEventId";
  if (string(chain->GetName()) == string("PicoDst")) {
    runid_name = "Event.mRunId";
    eventid_name = "Event.mEventId";
  }

  TTreeReaderValue<int> runid(reader, runid_name.c_str());
  TTreeReaderValue<int> eventid(reader, eventid_name.c_str());

  while (reader.Next()) {
    size_t entry = reader.GetCurrentEntry();
    key_type evt{*runid, *eventid};
    map[evt] = entry;
  }
  return map;
}

std::set<key_type> ExtractKeys(event_map& map) {
  std::set<key_type> ret;
  for (auto& entry : map) {
    ret.insert(entry.first);
  }
  return ret;
}


int main(int argc, char* argv[]) {
  // setup command line flags
  dijetcore::SetUsageMessage(
      "Compare tower/track matching in TStarJetPico and StPicoDst files");
  dijetcore::ParseCommandLineFlags(&argc, argv);

  // initialize logging - search for logging related command line flags
  dijetcore::InitLogging(argv[0]);

  // attempt to build directory for output files, if it does not already
  if (!dijetcore::CreateDirectory(FLAGS_outputDir)) {
    LOG(FATAL) << "Could not create output directory: " << FLAGS_outputDir;
  };

  // build our input chain of ROOT tree(s)
  TChain* chain_jetpico = dijetcore::NewChainFromInput(FLAGS_jetPicoFiles);
  TChain* chain_stpico =
      dijetcore::NewChainFromInput(FLAGS_stPicoFiles, "PicoDst");


  // get the event maps for both chains
  event_map jetpico_events = EnumerateEvents(chain_jetpico);
  event_map stpico_events = EnumerateEvents(chain_stpico);

  // get the overlap of the two event sets
  std::set<key_type> jetpico_keys = ExtractKeys(jetpico_events);
  std::set<key_type> stpico_keys = ExtractKeys(stpico_events);
  std::set<key_type> common_keys;
  for (auto& key : jetpico_keys) {
    if (stpico_keys.find(key) != stpico_keys.end()) {
      common_keys.insert(key);
    }
  }

  LOG(INFO) << "Number of matched events: " << common_keys.size();

  // initialize the JetTree reader
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  if (chain_jetpico != nullptr) reader->SetInputChain(chain_jetpico);
  reader->GetTowerCuts()->AddBadTowers(FLAGS_towList.c_str());
  // set hadronic correction
  reader->SetApplyFractionHadronicCorrection(true);
  reader->SetFractionHadronicCorrection(0.999);
  reader->SetRejectTowerElectrons(kFALSE);
  reader->GetTrackCuts()->SetDCACut(
      3.0);  // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(
      15);  // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(
      0.52);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(
      1000);  // essentially infinity - cut in eventcuts
  reader->GetTowerCuts()->SetMaxEtCut(
      1000);  // essentially infinity - cut in eventcuts
  reader->GetEventCuts()->SetMaxEventPtCut(1000);  // Set Maximum track Pt
  reader->GetEventCuts()->SetMaxEventEtCut(1000);  // Set Maximum tower Et
  reader->GetEventCuts()->SetVertexZCut(30);  // vertex z range (z = beam axis)
  reader->GetEventCuts()->SetVertexZDiffCut(3);  // cut on Vz - VPD Vz
  reader->Init(FLAGS_nEvents);

  // create file name from job name and (conditionally) the job id
  std::string filename = FLAGS_name;
  if (FLAGS_id < 0)
    filename += ".root";
  else
    filename = dijetcore::MakeString(filename, "_", FLAGS_id, ".root");
  std::string outputFileName =
      dijetcore::ConcatenatePath(FLAGS_outputDir, filename);

  // create output file - set to create or overwrite current file
  TFile outputFile(outputFileName.c_str(), "RECREATE");

  google::ShutdownGoogleLogging();
  google::ShutDownCommandLineFlags();
  return 0;
}