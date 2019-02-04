// General QA of Run 14 data - generates QA plots for 
// run-level, event-level and track-level variables
// of interest.

#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/worker/data/run14_qa_worker/run14_qa_worker.h"
#include "dijetcore/util/data/trigger_lookup.h"

#include <string>
#include <set>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"

DIJETCORE_DEFINE_string(input, "", "input root TFile (or txt file with TFile list) with JetTree(s)");
DIJETCORE_DEFINE_string(outputDir, "tmp", "directory to write output TFile in");
DIJETCORE_DEFINE_string(name, "run14_qa", "name for output files");
DIJETCORE_DEFINE_int(id, -1, "job id: used for batch jobs: keep negative to ignore for output file naming");
DIJETCORE_DEFINE_bool(doTrackQA, true, "flag to turn on/off track QA histograms");
DIJETCORE_DEFINE_bool(doTowerQA, true, "flag to turn on/off tower QA histograms");
DIJETCORE_DEFINE_string(histogramPrefix, "", "can be used to set a prefix for all histogram names");
DIJETCORE_DEFINE_string(runIDFile, "", "file containing TTree 'runid' that contains list of run 14 runids");
DIJETCORE_DEFINE_string(triggers, "y14vpdmb30", "trigger selection for accepted events - defaults to VPD MB 30, a minbias trigger");
DIJETCORE_DEFINE_string(towList, "resources/bad_tower_lists/empty_list.txt", "bad tower list");

int main(int argc, char* argv[]) {
  // setup command line flags
  dijetcore::SetUsageMessage("Run 14 QA routine for run level, event level and track/tower level variables.\n \
                              Run QA: luminosity, Number of events, average vertex position, average number of tracks\n \
                              Event QA: multiplicity, vertex, highest energy track/tower \n \
                              Track QA: momentum, dca, nhit, nhit/nhitposs, etc \n \
                              Tower QA: energy, transverse energy, rate of activity above threshold, ADC");
  dijetcore::ParseCommandLineFlags(&argc, argv);
  
  // initialize logging - search for logging related command line flags 
  dijetcore::InitLogging(argv[0]);

  // attempt to build directory for output files, if it does not already
  if (!dijetcore::CreateDirectory(FLAGS_outputDir)) {
    LOG(FATAL) << "Could not create output directory: " << FLAGS_outputDir;
  };

  // build our input chain of ROOT tree(s)
  TChain* chain = dijetcore::NewChainFromInput(FLAGS_input);

  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  dijetcore::InitReaderWithDefaults(reader, chain, "submit/empty_list.txt", "");
  reader->GetTrackCuts()->SetDCACut(3.0);                // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(15);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.52);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  reader->GetTowerCuts()->SetMaxEtCut(1000);             // essentially infinity - cut in eventcuts
  reader->GetTowerCuts()->AddBadTowers(FLAGS_towList.c_str()); // bad tower list
  reader->GetEventCuts()->SetMaxEventPtCut(1000);        // Set Maximum track Pt
  reader->GetEventCuts()->SetMaxEventEtCut(1000);        // Set Maximum tower Et
  reader->GetEventCuts()->SetVertexZCut(30);             // vertex z range (z = beam axis)
  reader->GetEventCuts()->SetVertexZDiffCut(3);          // cut on Vz - VPD Vz
    
  // create file name from job name and (conditionally) the job id 
  std::string filename = FLAGS_name;
  if (FLAGS_id < 0)
    filename += ".root";
  else 
    filename = dijetcore::MakeString(filename, "_", FLAGS_id, ".root");
  std::string outputFileName = dijetcore::ConcatenatePath(FLAGS_outputDir, filename);

  // create output file - set to create or overwrite current file
  TFile outputFile(outputFileName.c_str(), "RECREATE");

  // initialize worker
  dijetcore::Run14QAWorker worker;
  worker.DoTowerQA(FLAGS_doTowerQA);
  worker.DoTrackQA(FLAGS_doTrackQA);

  // check if we are going to do runID dependent QA, if so, load runIDs
  if (!FLAGS_runIDFile.empty()) {
    TFile runIDFile(FLAGS_runIDFile.c_str(), "READ");
    TTreeReader runid_reader("runid", &runIDFile);
    TTreeReaderValue<unsigned> runid(runid_reader, "runid");
    
    // create set to naturally eliminate duplicates
    std::set<unsigned> runid_set;
    while (runid_reader.Next())
      runid_set.insert(*runid);
    
    // add runids to worker
    worker.DoRunQA(runid_set);
  }

  // initialize worker
  worker.Init(FLAGS_histogramPrefix);

  // get selected triggers
  std::set<unsigned> triggers = dijetcore::GetTriggerIDs(FLAGS_triggers);

  // start the event loop
  // --------------------
  try {
    while (reader->NextEvent()) {
      // Print out reader status every 10 seconds
      reader->PrintStatus(10);
      TStarJetPicoEventHeader* header = reader->GetEvent()->GetHeader();

      // check if event fired a trigger we will use
      if (triggers.size() != 0) {
        bool use_event = false;
        for (auto trigger : triggers)
          if (header->HasTriggerId(trigger))
            use_event = true;
        if (!use_event) continue;
      }

      // run
      worker.Run(*reader);

    }
  } catch(std::exception& e) {
    std::cerr << "Caught: " << e.what() << " during analysis loop." << std::endl;
  }

  // write to file
  worker.WriteTo(outputFile);
  outputFile.Close();

  return 0;
}