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
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

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

  LOG(INFO) << "Got here";
  // build our input chain of ROOT tree(s)
  TChain* chain = dijetcore::NewChainFromInput(FLAGS_tree);
  LOG(INFO) << "Got here";
  // initialize the reader
  TTreeReader reader(chain);
  TTreeReaderValue<int> runid_value(reader, "fEventHeader.fRunId");
  LOG(INFO) << "Got here reader: " << reader.GetEntries();
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
  LOG(INFO) << "Got here loop";
  try {
    while (reader.Next()) {
      LOG(INFO) << "in loop";
      all_runids.insert(*runid_value);
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