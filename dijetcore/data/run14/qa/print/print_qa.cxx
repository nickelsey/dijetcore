#include "dijetcore/util/root/root_print_routines.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/flags.h"

#include <string>
#include <set>

#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

DIJETCORE_DEFINE_string(input, "", "root file with results from run_qa");
DIJETCORE_DEFINE_string(outputDir, "results", "directory for output");
DIJETCORE_DEFINE_string(histPrefixes, "", "histogram prefixes, separated by commas");

using std::string;

int main(int argc, char* argv[]) {
  
  string usage = "Run 7 differential di-jet imbalance print routine";
  
  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);
  
  // turn off print messages
  gErrorIgnoreLevel = kInfo+1;

  // load input file
  if (FLAGS_input.empty()) {
    LOG(ERROR) << "no input file specified";
    return 1;
  }

  TFile input_file(FLAGS_input.c_str(), "READ");

  if (!input_file.IsOpen()) {
    LOG(ERROR) << "Could not open input file: " << FLAGS_input << ", exiting";
    return 1;
  }

  // load histograms
  // first, find all histogram prefixes
  std::set<string> prefixes;
  dijetcore::SplitString(FLAGS_histPrefixes, prefixes, ",");

  gflags::ShutDownCommandLineFlags();
  return 0;
}