// routine for generating a TTree of all unique RunIDs
// in a set of TStarJetPicoEvents

#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/data/reader_util.h"

#include <string>
#include <set>

#include "TFile.h"
#include "TTree.h"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"

DIJETCORE_DEFINE_string(tree, "", "input root TFile (or txt file with TFile list) with JetTree(s)");
DIJETCORE_DEFINE_string(outputDir, "tmp", "directory to write output TFile in");
DIJETCORE_DEFINE_string(jobName, "run14_qa", "name for output files");
DIJETCORE_DEFINE_int(id, -1, "job id: used for batch jobs: keep negative to ignore for output file naming");

int main(int argc, char* argv[]) {
  // setup command line flags
  dijetcore::SetUsageMessage("Generate a tree of all unique runids");
  dijetcore::ParseCommandLineFlags(&argc, argv);
  
  // initialize logging - search for logging related command line flags 
  dijetcore::InitLogging(argv[0]);

  // attempt to build directory for output files, if it does not already
  if (!dijetcore::CreateDirectory(FLAGS_outputDir)) {
    LOG(FATAL) << "Could not create output directory: " << FLAGS_outputDir;
  };

  // build our input chain of ROOT tree(s)
  TChain* chain = dijetcore::NewChainFromInput(FLAGS_tree);

  // initialize the reader(s)
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  dijetcore::InitReaderWithDefaults(reader, chain, "submit/empty_list.txt", "");
  reader->GetTrackCuts()->SetDCACut(3.0);                // distance of closest approach to primary vtx
  reader->GetTrackCuts()->SetMinNFitPointsCut(15);       // minimum fit points in track reco
  reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.52);  // minimum ratio of fit points used over possible
  reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
  reader->GetTowerCuts()->SetMaxEtCut(1000);             // essentially infinity - cut in eventcuts
  reader->GetEventCuts()->SetMaxEventPtCut(1000);        // Set Maximum track Pt
  reader->GetEventCuts()->SetMaxEventEtCut(1000);        // Set Maximum tower Et
  reader->GetEventCuts()->SetVertexZCut(30);             // vertex z range (z = beam axis)
  reader->GetEventCuts()->SetVertexZDiffCut(3);          // cut on Vz - VPD Vz
    
  // create file name from job name and (conditionally) the job id 
  std::string filename = FLAGS_jobName;
  if (FLAGS_id < 0)
    filename += ".root";
  else 
    filename = dijetcore::MakeString(filename, "_", FLAGS_id, ".root");
  std::string outputFileName = dijetcore::ConcatenatePath(FLAGS_outputDir, filename);

  // create output file - set to create or overwrite current file
  TFile outputFile(outputFileName.c_str(), "RECREATE");

  // create set to contain all runids
  std::set<unsigned> all_runids;

  try {
    while(reader->NextEvent()) {
      all_runids.insert(reader->GetEvent()->GetHeader()->GetRunId());
    }
  } catch(std::exception& e) {
    std::cerr << "Caught: " << e.what() << " during analysis loop." << std::endl;
  }

  // create output TTree
  TTree* tree = new TTree("runid", "runid");
  unsigned runid_val = 0;
  tree->Branch("runid", &runid_val);
  for (auto& val : all_runids) {
    runid_val = val;
    tree->Fill();
  }

  outputFile.cd();
  tree->Write();
  outputFile.Close();

  return 0;
}