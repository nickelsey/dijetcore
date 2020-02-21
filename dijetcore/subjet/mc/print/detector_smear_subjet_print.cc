// Print routines for detector smearing subjet analysis
// comparing resolution as a function of other jet parameters

#include "dijetcore/lib/containers.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/root/grid_print_utils.h"
#include "dijetcore/util/root/histogram_utils.h"
#include "dijetcore/util/root/root_print_utils.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <string>
#include <vector>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "",
                        "input root file from detector_smear_subjet");
DIJETCORE_DEFINE_string(outputdir, "tmp", "output directory to print results");

std::vector<std::string>
get_bin_labels(std::vector<std::pair<double, double>> &v, std::string var, std::string units) {
  std::vector<std::string> ret;
  for (auto &e : v) {
    ret.push_back(dijetcore::MakeString(e.first, " < ", var, " < ", e.second, " ", units));
  }
  return ret;
}

void split_and_print_2d(TH2D *hist,
                        std::vector<std::pair<double, double>> &bins,
                        std::vector<std::string> &bin_names,
                        std::string outputdir, std::string name,
                        dijetcore::histogramOpts hopts,
                        dijetcore::canvasOpts copts,
                        std::string x_title, std::string y_title) {
  std::vector<TH1D *> hists;
  for (auto &bin : bins) {
    std::string hist_name =
        dijetcore::MakeString(hist->GetName(), "_", bin.first, "_", bin.second);
    hist->GetXaxis()->SetRangeUser(bin.first, bin.second);
    TH1D* tmp = hist->ProjectionY(hist_name.c_str());
    tmp->GetXaxis()->SetRangeUser(-0.5, 0.5);
    tmp->Scale(1.0 / tmp->Integral());
    hists.push_back(tmp);
  }
  hist->GetXaxis()->SetRange(0, -1);


  dijetcore::Overlay1D(hists, bin_names, hopts, copts, outputdir, name, "", x_title, y_title, "");
}

void print_profiles(TH2D* split1, TH2D* split2, std::string outputdir, std::string name, dijetcore::histogramOpts hopts, dijetcore::canvasOpts copts, std::string x_title, std::string y_title) {
  TProfile* p1 = split1->ProfileX("tmp1", 1, -1, "");
  TProfile* p2 = split2->ProfileX("tmp2", 1, -1, "");
  p1->GetYaxis()->SetRangeUser(-0.2, 0.2);

  dijetcore::Overlay1D(p1, p2, "split 1", "split 2", hopts, copts, outputdir, name, "", x_title, y_title);
}

int main(int argc, char *argv[]) {

  string usage = "Print analysis of variation in subjet quantities due to "
                 "detector response.";

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

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outputdir);
  boost::filesystem::create_directories(dir);

  // check to make sure the input file exists
  if (!boost::filesystem::exists(FLAGS_input)) {
    std::cerr << "input file does not exist: " << FLAGS_input << std::endl;
    return 1;
  }

  TFile input(FLAGS_input.c_str(), "READ");

  // load all histograms
  TH2D *ptzgsmear1 = (TH2D *)input.Get("ptzgsmear1");
  TH2D *ptzgsmear2 = (TH2D *)input.Get("ptzgsmear2");
  TH2D *ptrgsmear1 = (TH2D *)input.Get("ptrgsmear1");
  TH2D *ptrgsmear2 = (TH2D *)input.Get("ptrgsmear2");
  TH2D *pttfsmear1 = (TH2D *)input.Get("pttfsmear1");
  TH2D *pttfsmear2 = (TH2D *)input.Get("pttfsmear2");
  TH2D *zgzgsmear1 = (TH2D *)input.Get("zgzgsmear1");
  TH2D *zgzgsmear2 = (TH2D *)input.Get("zgzgsmear2");
  TH2D *zgrgsmear1 = (TH2D *)input.Get("zgrgsmear1");
  TH2D *zgrgsmear2 = (TH2D *)input.Get("zgrgsmear2");
  TH2D *zgtfsmear1 = (TH2D *)input.Get("zgtfsmear1");
  TH2D *zgtfsmear2 = (TH2D *)input.Get("zgtfsmear2");
  TH2D *rgzgsmear1 = (TH2D *)input.Get("rgzgsmear1");
  TH2D *rgzgsmear2 = (TH2D *)input.Get("rgzgsmear2");
  TH2D *rgrgsmear1 = (TH2D *)input.Get("rgrgsmear1");
  TH2D *rgrgsmear2 = (TH2D *)input.Get("rgrgsmear2");
  TH2D *rgtfsmear1 = (TH2D *)input.Get("rgtfsmear1");
  TH2D *rgtfsmear2 = (TH2D *)input.Get("rgtfsmear2");

  // set x axis splits for each variable
  std::vector<std::pair<double, double>> pt_splits{
      {15, 20}, {20, 25}, {25, 30}, {30, 40}};
  std::vector<std::string> pt_split_strings =
      get_bin_labels(pt_splits, "p_{T}", "GeV/c");
  std::vector<std::pair<double, double>> zg_splits{
      {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}};
  std::vector<std::string> zg_split_strings = get_bin_labels(zg_splits, "z_{g}", "");
  std::vector<std::pair<double, double>> rg_splits{
      {0.0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}};
  std::vector<std::string> rg_split_strings = get_bin_labels(rg_splits, "R_{g}", "");

  // create histogram and canvas opts
  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;


  // print the individual distributions
  split_and_print_2d(ptzgsmear1, pt_splits, pt_split_strings, FLAGS_outputdir, "pt_zg_split_1", hopts, copts, "#Delta z_{g}/z_{g}", "fraction");
  split_and_print_2d(ptzgsmear2, pt_splits, pt_split_strings, FLAGS_outputdir, "pt_zg_split_2", hopts, copts, "#Delta z_{g}/z_{g}", "fraction");
  split_and_print_2d(ptrgsmear1, pt_splits, pt_split_strings, FLAGS_outputdir, "pt_rg_split_1", hopts, copts, "#Delta R_{g}/R_{g}", "fraction");
  split_and_print_2d(ptrgsmear2, pt_splits, pt_split_strings, FLAGS_outputdir, "pt_rg_split_2", hopts, copts, "#Delta R_{g}/R_{g}", "fraction");
  split_and_print_2d(pttfsmear1, pt_splits, pt_split_strings, FLAGS_outputdir, "pt_tf_split_1", hopts, copts, "#Delta t_{f}/t_{f}", "fraction");
  split_and_print_2d(pttfsmear2, pt_splits, pt_split_strings, FLAGS_outputdir, "pt_tf_split_2", hopts, copts, "#Delta t_{f}/t_{f}", "fraction");

  split_and_print_2d(zgzgsmear1, zg_splits, zg_split_strings, FLAGS_outputdir, "zg_zg_split_1", hopts, copts, "#Delta z_{g}/z_{g}", "fraction");
  split_and_print_2d(zgzgsmear2, zg_splits, zg_split_strings, FLAGS_outputdir, "zg_zg_split_2", hopts, copts, "#Delta z_{g}/z_{g}", "fraction");
  split_and_print_2d(zgrgsmear1, zg_splits, zg_split_strings, FLAGS_outputdir, "zg_rg_split_1", hopts, copts, "#Delta R_{g}/R_{g}", "fraction");
  split_and_print_2d(zgrgsmear2, zg_splits, zg_split_strings, FLAGS_outputdir, "zg_rg_split_2", hopts, copts, "#Delta R_{g}/R_{g}", "fraction");
  split_and_print_2d(zgtfsmear1, zg_splits, zg_split_strings, FLAGS_outputdir, "zg_tf_split_1", hopts, copts, "#Delta t_{f}/t_{f}", "fraction");
  split_and_print_2d(zgtfsmear2, zg_splits, zg_split_strings, FLAGS_outputdir, "zg_tf_split_2", hopts, copts, "#Delta t_{f}/t_{f}", "fraction");

  split_and_print_2d(rgzgsmear1, rg_splits, rg_split_strings, FLAGS_outputdir, "rg_zg_split_1", hopts, copts, "#Delta z_{g}/z_{g}", "fraction");
  split_and_print_2d(rgzgsmear2, rg_splits, rg_split_strings, FLAGS_outputdir, "rg_zg_split_2", hopts, copts, "#Delta z_{g}/z_{g}", "fraction");
  split_and_print_2d(rgrgsmear1, rg_splits, rg_split_strings, FLAGS_outputdir, "rg_rg_split_1", hopts, copts, "#Delta R_{g}/R_{g}", "fraction");
  split_and_print_2d(rgrgsmear2, rg_splits, rg_split_strings, FLAGS_outputdir, "rg_rg_split_2", hopts, copts, "#Delta R_{g}/R_{g}", "fraction");
  split_and_print_2d(rgtfsmear1, rg_splits, rg_split_strings, FLAGS_outputdir, "rg_tf_split_1", hopts, copts, "#Delta t_{f}/t_{f}", "fraction");
  split_and_print_2d(rgtfsmear2, rg_splits, rg_split_strings, FLAGS_outputdir, "rg_tf_split_2", hopts, copts, "#Delta t_{f}/t_{f}", "fraction");

  // print the profiles
  print_profiles(ptzgsmear1, ptzgsmear2, FLAGS_outputdir, "ptzg", hopts, copts, "p_{T}", "#Delta z_{g}/z_{g}");
  print_profiles(ptrgsmear1, ptrgsmear2, FLAGS_outputdir, "ptrg", hopts, copts, "p_{T}", "#Delta R_{g}/R_{g}");
  print_profiles(pttfsmear1, pttfsmear2, FLAGS_outputdir, "pttf", hopts, copts, "p_{T}", "#Delta t_{f}/t_{f}");

  print_profiles(zgzgsmear1, zgzgsmear2, FLAGS_outputdir, "zgzg", hopts, copts, "z_{g}", "#Delta z_{g}/z_{g}");
  print_profiles(zgrgsmear1, zgrgsmear2, FLAGS_outputdir, "zgrg", hopts, copts, "z_{g}", "#Delta R_{g}/R_{g}");
  print_profiles(zgtfsmear1, zgtfsmear2, FLAGS_outputdir, "zgtf", hopts, copts, "z_{g}", "#Delta t_{f}/t_{f}");

  print_profiles(rgzgsmear1, rgzgsmear2, FLAGS_outputdir, "rgzg", hopts, copts, "R_{g}", "#Delta z_{g}/z_{g}");
  print_profiles(rgrgsmear1, rgrgsmear2, FLAGS_outputdir, "rgrg", hopts, copts, "R_{g}", "#Delta R_{g}/R_{g}");
  print_profiles(rgtfsmear1, rgtfsmear2, FLAGS_outputdir, "rgtf", hopts, copts, "R_{g}", "#Delta t_{f}/t_{f}");

  gflags::ShutDownCommandLineFlags();
  return 0;
}
