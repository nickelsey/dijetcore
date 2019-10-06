// plots comparison of rho/sigma for various constituent pT cuts
// requires 16-8 gev leading/subleading pT, anti-kt, R=0.4 di-jets,
// and compares 0.2 GeV, 1.0 GeV, and 2.0 GeV
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/data/centrality/centrality_run7.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/trigger_lookup.h"
#include "dijetcore/util/data/vector_conversion.h"
#include "dijetcore/util/root/root_print_utils.h"

#include <string>

#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

DIJETCORE_DEFINE_string(input, "", "input file with result trees");
DIJETCORE_DEFINE_string(outputDir, "tmp", "output location");

int main(int argc, char *argv[]) {
  std::string usage = "prints comparisons of rho and sigma for di-jets with "
                      "different constituent pT cuts";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outputDir);
  boost::filesystem::create_directories(dir);

  // set drawing preferences for histograms and graphs
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  TFile input(FLAGS_input.c_str(), "READ");

  std::string dijet_key_2_0 =
      "LEAD_INIT_R_0.4_alg_2_pt_16_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_"
      "0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_8_const_eta_1_const_"
      "pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2";
  std::string dijet_key_1_0 =
      "LEAD_INIT_R_0.4_alg_2_pt_16_const_eta_1_const_pt_1_MATCH_R_0.4_alg_2_pt_"
      "0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_8_const_eta_1_const_"
      "pt_1_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2";

  TTreeReader tree_2_0(dijet_key_2_0.c_str(), &input);
  TTreeReaderValue<int> run_2(tree_2_0, "runid");
  TTreeReaderValue<int> event_2(tree_2_0, "eventid");
  TTreeReaderValue<int> refmult_2(tree_2_0, "refmult");
  TTreeReaderValue<int> cent_2(tree_2_0, "cent");
  TTreeReaderValue<TLorentzVector> jl_2(tree_2_0, "jl");
  TTreeReaderValue<TLorentzVector> js_2(tree_2_0, "js");
  TTreeReaderValue<TLorentzVector> jlm_2(tree_2_0, "jlm");
  TTreeReaderValue<TLorentzVector> jsm_2(tree_2_0, "jsm");
  TTreeReaderValue<double> rho_2(tree_2_0, "jlmrho");
  TTreeReaderValue<double> rho_2_hard(tree_2_0, "jlrho");
  TTreeReaderValue<double> sig_2(tree_2_0, "jlmsig");
  TTreeReaderValue<double> sig_2_hard(tree_2_0, "jlsig");

  TTreeReader tree_1_0(dijet_key_1_0.c_str(), &input);
  TTreeReaderValue<int> run_1(tree_1_0, "runid");
  TTreeReaderValue<int> event_1(tree_1_0, "eventid");
  TTreeReaderValue<int> refmult_1(tree_1_0, "refmult");
  TTreeReaderValue<int> cent_1(tree_1_0, "cent");
  TTreeReaderValue<TLorentzVector> jl_1(tree_1_0, "jl");
  TTreeReaderValue<TLorentzVector> js_1(tree_1_0, "js");
  TTreeReaderValue<TLorentzVector> jlm_1(tree_1_0, "jlm");
  TTreeReaderValue<TLorentzVector> jsm_1(tree_1_0, "jsm");
  TTreeReaderValue<double> rho_1(tree_1_0, "jlmrho");
  TTreeReaderValue<double> rho_1_hard(tree_1_0, "jlrho");
  TTreeReaderValue<double> sig_1(tree_1_0, "jlmsig");
  TTreeReaderValue<double> sig_1_hard(tree_1_0, "jlsig");

  // define the histograms
  TH1D *rho_2_0 = new TH1D("rho20", ";#rho;event fraction", 100, 0, 100);
  TH1D *rho_1_0 = new TH1D("rho10", ";#rho;event fraction", 100, 0, 100);
  TH1D *rho_0_2 = new TH1D("rho02", ";#rho;event fraction", 100, 0, 100);
  TH1D *sig_2_0 = new TH1D("sig20", ";#sigma;event fraction", 100, 0, 20);
  TH1D *sig_1_0 = new TH1D("sig10", ";#sigma;event fraction", 100, 0, 20);
  TH1D *sig_0_2 = new TH1D("sig02", ";#sigma;event fraction", 100, 0, 20);

  while (tree_2_0.Next()) {
    if (*cent_2 < 0 || *cent_2 >= 3)
      continue;
    rho_0_2->Fill(*rho_2);
    rho_2_0->Fill(*rho_2_hard);
    sig_0_2->Fill(*sig_2);
    sig_2_0->Fill(*sig_2_hard);
  }

  while (tree_1_0.Next()) {
    if (*cent_1 < 0 || *cent_1 >= 3)
      continue;
    rho_1_0->Fill(*rho_1_hard);
    sig_1_0->Fill(*sig_1_hard);
  }

  rho_2_0->Scale(1.0 / rho_2_0->Integral());
  rho_1_0->Scale(1.0 / rho_1_0->Integral());
  rho_0_2->Scale(1.0 / rho_0_2->Integral());
  sig_2_0->Scale(1.0 / sig_2_0->Integral());
  sig_1_0->Scale(1.0 / sig_1_0->Integral());
  sig_0_2->Scale(1.0 / sig_0_2->Integral());

  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;
  copts.log_y = true;

  std::string rho_x_axis = "#rho [GeV/A]";
  std::string sigma_x_axis = "#sigma [GeV/A]";
  std::string y_axis = "event fraction";
  std::vector<TH1D *> rho_hists{rho_2_0, rho_1_0, rho_0_2};
  std::vector<TH1D *> sigma_hists{sig_2_0, sig_1_0, sig_0_2};
  std::vector<std::string> hist_labels{"p_{T}^{const}>2.0 GeV/c",
                                       "p_{T}^{const}>1.0 GeV/c",
                                       "p_{T}^{const}>0.2 GeV/c"};
  dijetcore::Overlay1D(rho_hists, hist_labels, hopts, copts, FLAGS_outputDir, "pt_dep_rho", "", rho_x_axis, y_axis, "", false);
  dijetcore::Overlay1D(sigma_hists, hist_labels, hopts, copts, FLAGS_outputDir, "pt_dep_sigma", "", sigma_x_axis, y_axis, "", false);
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}
