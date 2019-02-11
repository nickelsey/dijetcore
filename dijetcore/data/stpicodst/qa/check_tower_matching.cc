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

#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoReader.h"
#include "TStarJetPicoTower.h"
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
DIJETCORE_DEFINE_int(nEvents, 100, "number of events to compare");
DIJETCORE_DEFINE_string(towList, "resources/bad_tower_lists/empty_list.txt",
                        "bad tower list");

using dijetcore::dijetcore_map;
using dijetcore::PairHash;
using std::pair;
using std::string;

using key_type = pair<unsigned, unsigned>;
using val_type = size_t;
using event_map = dijetcore_map<key_type, val_type, PairHash>;

struct Tower {
  int id = 0;
  int id2 = 0;
  int id3 = 0;
  double eta = 0;
  double phi = 0;
  double e;
};

struct Track {
  int id = 0;
  double eta = 0;
  double phi = 0;
  double pt = 0;
  Tower tow;
};

bool MatchTracks(Track& t1, Track& t2) {
  if (fabs(t1.eta - t2.eta) < 0.01 && fabs(t1.phi - t2.phi) < 0.01 &&
      fabs(t1.pt - t2.pt) < 0.01)
    return true;
  return false;
}

// enumerate events for each dataset, returns a map where each key is
// an event {runid, eventid} pair, and the value is the entry position in the
// TChain
event_map EnumerateEventsTStarJetPico(TChain* chain) {
  event_map map;
  TTreeReader reader(chain);

  string runid_name = "fEventHeader.fRunId";
  string eventid_name = "fEventHeader.fEventId";

  TTreeReaderValue<int> runid(reader, runid_name.c_str());
  TTreeReaderValue<int> eventid(reader, eventid_name.c_str());

  while (reader.Next()) {
    size_t entry = reader.GetCurrentEntry();
    key_type evt{*runid, *eventid};
    map[evt] = entry;
  }
  return map;
}

event_map EnumerateEventsStPicoDst(TChain* chain) {
  event_map map;
  TTreeReader reader(chain);

  string event_name = "Event";

  TTreeReaderArray<StPicoEvent> event(reader, event_name.c_str());

  while (reader.Next()) {
    size_t entry = reader.GetCurrentEntry();
    StPicoEvent& evt = event[0];

    key_type event_key{evt.runId(), evt.eventId()};
    map[event_key] = entry;
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
  event_map jetpico_events = EnumerateEventsTStarJetPico(chain_jetpico);

  event_map stpico_events = EnumerateEventsStPicoDst(chain_stpico);

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
  reader->Init(-1);

  // build StPicoDstReader
  StPicoDstReader* pico_reader = new StPicoDstReader(FLAGS_stPicoFiles.c_str());
  pico_reader->Init();

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

  int event_nr = 0;
  int tracks = 0;
  int matches = 0;
  for (auto& key : common_keys) {
  
    if (event_nr >= FLAGS_nEvents) break;

    // read in the matched events
    size_t jet_entry = jetpico_events[key];
    size_t pico_entry = stpico_events[key];

    reader->ReadEvent(jet_entry);
    pico_reader->tree()->GetEntry(pico_entry);

    // get all matched tracks from TStarJetPicos
    std::vector<Track> jetpico_tracks;
    jetpico_tracks.reserve(
        reader->GetEvent()->GetHeader()->GetNOfPrimaryTracks());

    for (Int_t ntower = 0;
         ntower < reader->GetEvent()->GetHeader()->GetNOfTowers(); ntower++) {
      TStarJetPicoTower* ptower = reader->GetEvent()->GetTower(ntower);
      for (Int_t ntrack = 0; ntrack < ptower->GetNAssocTracks(); ntrack++) {
        Int_t idx = ptower->GetMatchedTrackIndex(ntrack);
        TStarJetPicoPrimaryTrack* ptrack =
            reader->GetEvent()->GetPrimaryTrack(idx);
        Track jet_track;
        jet_track.id = idx;
        jet_track.eta = ptrack->GetEta();
        jet_track.phi = ptrack->GetPhi();
        jet_track.pt = ptrack->GetPt();
        jet_track.tow.id = ptower->GetId();
        jet_track.tow.eta = ptower->GetEta();
        jet_track.tow.phi = ptower->GetPhi();
        jet_track.tow.e = ptower->GetEnergy();

        jetpico_tracks.push_back(jet_track);
      }
    }

    // get all matched tracks from StPicos
    std::vector<Track> stpico_tracks;
    stpico_tracks.reserve(pico_reader->picoDst()->numberOfBEmcPidTraits());
    
    for (int i = 0; i < pico_reader->picoDst()->numberOfBEmcPidTraits(); ++i) {
      
      StPicoBEmcPidTraits* pid = pico_reader->picoDst()->bemcPidTraits(i);
      StPicoTrack* track = pico_reader->picoDst()->track(pid->trackIndex());
      StPicoBTowHit* tower = pico_reader->picoDst()->btowHit(pid->btowId());

      if (!track->isPrimary())
        continue;

      Track pico_track;
      pico_track.id = pid->trackIndex();
      pico_track.eta = track->pMom().Eta();
      pico_track.phi = track->pMom().Phi();
      pico_track.pt = track->pPt();
      pico_track.tow.id = pid->btowId();
      pico_track.tow.id2 = pid->btowId2();
      pico_track.tow.id3 = pid->btowId3();
      pico_track.tow.eta = pid->btowEtaDist();
      pico_track.tow.phi = pid->btowPhiDist();
      pico_track.tow.e = pid->btowE();

      stpico_tracks.push_back(pico_track);
    }

    // now find matches - we'll start with stpicos because they have tighter cuts
    for (auto& sttrack : stpico_tracks) {
      for (auto& jettrack : jetpico_tracks) {
        if (MatchTracks(sttrack, jettrack)) {
          tracks++;
          if (sttrack.tow.id == jettrack.tow.id ||
              sttrack.tow.id2 == jettrack.tow.id ||
              sttrack.tow.id3 == jettrack.tow.id) {
            matches++; 
          }
        }
      }
    }

    event_nr++;
  }

  LOG(INFO) << "percentage of matched towers: " << (double) matches / tracks;

  google::ShutdownGoogleLogging();
  google::ShutDownCommandLineFlags();
  return 0;
}
