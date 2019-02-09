// routine for generating a comparison of refmult and grefmult
// for run 14 productions

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
DIJETCORE_DEFINE_string(jobName, "run14_refmult_qa", "name for output files");
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

  // initialize the reader
  TTreeReader reader(chain);
  TTreeReaderValue<int> runid_value(reader, "fEventHeader.fRunId");
  TTreeReaderValue<int> eventid_value(reader, "fEventHeader.fEventId");
  TTreeReaderValue<int> refmult_value(reader, "fEventHeader.fRefMult");
  TTreeReaderValue<int> grefmult_value(reader, "fEventHeader.fGRefMult");

  // create file name from job name and (conditionally) the job id 
  std::string filename = FLAGS_jobName;
  if (FLAGS_id < 0)
    filename += ".root";
  else 
    filename = dijetcore::MakeString(filename, "_", FLAGS_id, ".root");
  std::string outputFileName = dijetcore::ConcatenatePath(FLAGS_outputDir, filename);

  // create output file - set to create or overwrite current file
  TFile outputFile(outputFileName.c_str(), "RECREATE");

  // create output TTree
  TTree* tree = new TTree("refmult", "refmult tree");
  unsigned runid = 0;
  unsigned eventid = 0;
  unsigned refmult = 0;
  unsigned grefmult = 0;
  tree->Branch("runid", &runid);
  tree->Branch("eventid", &eventid);
  tree->Branch("refmult", &refmult);
  tree->Branch("grefmult", &grefmult);
  
  try {
    while (reader.Next()) {
      runid = *runid_value;
      eventid = *eventid_value;
      refmult = *refmult_value;
      grefmult = *grefmult_value;
      tree->Fill();

      if (*grefmult_value > 1000) {
        LOG(INFO) << "grefmult above 1000: " << *grefmult_value;
        LOG(INFO) << "runid: " << *runid_value << " eventid: " << *eventid_value;
      
    }
  } catch(std::exception& e) {
    std::cerr << "Caught: " << e.what() << " during analysis loop." << std::endl;
  }

  outputFile.cd();
  tree->Write();
  outputFile.Close();

  return 0;
}