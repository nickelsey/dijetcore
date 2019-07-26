// generate efficiency curves from embedding data: two sets of histograms are
// given, mc tracks and reconstructed tracks in bins of luminosity, centrality,
// pt, eta and phi this division of reco / mc gives effective efficiency

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/root/root_print_utils.h"

#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TStyle.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

struct Options {
  string name = "y14_effic"; /* name for output root file */
  string input = "";      /* root file, should contain both y6 & y12 data (with
                             different histogram      prefixes) */
  string out_dir = "tmp"; /* directory to save output in */
  int lumi_bins = 3;      /* number of bins in luminosity */
  int cent_bins = 16;     /* number of bins in centrality */
  double max_pt = 5.0;    /* maximum pt cutoff */
};

DIJETCORE_DEFINE_string(name, "y14_effic", "name for output root file");
DIJETCORE_DEFINE_string(input, "", "input root file");
DIJETCORE_DEFINE_string(outdir, "tmp", "output directory for results");
DIJETCORE_DEFINE_int(lumiBins, 3, "number of luminosity bins");
DIJETCORE_DEFINE_int(centBins, 16, "number of centrality bins");
DIJETCORE_DEFINE_double(maxPt, 5.0, "maximum track pT");

int main(int argc, char *argv[]) {

  // setup command line flags
  dijetcore::SetUsageMessage("Generate efficiency curves for run14");
  dijetcore::ParseCommandLineFlags(&argc, argv);

  // initialize logging - search for logging related command line flags
  dijetcore::InitLogging(argv[0]);

  if (FLAGS_outdir.empty())
    FLAGS_outdir = "tmp";
  boost::filesystem::path dir(FLAGS_outdir.c_str());
  boost::filesystem::create_directories(dir);

  std::vector<string> refcent_string{"0-5%",   "5-10%",  "10-15%", "15-20%",
                                     "20-25%", "25-30%", "30-35%", "35-40%",
                                     "40-45%", "45-50%", "50-55%", "55-60%",
                                     "60-65%", "65-70%", "70-75%", "75-80%"};
  std::vector<string> lumi_string{"0-33 khz", "33-66 khz", "66-100 khz"};

  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);

  // check to make sure the input file exists
  if (!boost::filesystem::exists(FLAGS_input)) {
    std::cerr << "input file does not exist: " << FLAGS_input << std::endl;
    ;
    return 1;
  }

  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outdir + "/" + FLAGS_name + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // load input file and read in histograms
  TFile in(FLAGS_input.c_str(), "READ");
  std::vector<std::vector<TH3D *>> in_mc;
  std::vector<std::vector<TH3D *>> in_match;

  for (int lumi = 0; lumi < FLAGS_lumiBins; ++lumi) {
    in_mc.push_back(std::vector<TH3D *>());
    in_match.push_back(std::vector<TH3D *>());
    for (int cent = 0; cent < FLAGS_centBins; ++cent) {

      string mc_name =
          "mc_lumi_" + std::to_string(lumi) + "_cent_" + std::to_string(cent);
      string match_name = "match_lumi_" + std::to_string(lumi) + "_cent_" +
                          std::to_string(cent);

      in_mc[lumi].push_back((TH3D *)in.Get(mc_name.c_str()));
      in_match[lumi].push_back((TH3D *)in.Get(match_name.c_str()));

      in_mc[lumi][cent]->GetXaxis()->SetRange(1, 50);
      in_match[lumi][cent]->GetXaxis()->SetRange(1, 50);
    }
  }

  // get 2D histograms - integrate out phi dependence for now
  std::vector<std::vector<TH2D *>> mc_2d;
  std::vector<std::vector<TH2D *>> match_2d;
  std::vector<std::vector<TH2D *>> eff_curves;
  std::vector<std::vector<TH2D *>> eff_curves_ratio;
  std::vector<std::vector<TH1D *>> eff_curves_1d;
  std::vector<std::vector<TH1D *>> eff_curves_1d_ratio;
  TCanvas c;
  for (int i = 0; i < FLAGS_lumiBins; ++i) {

    mc_2d.push_back(std::vector<TH2D *>());
    match_2d.push_back(std::vector<TH2D *>());
    eff_curves.push_back(std::vector<TH2D *>());
    eff_curves_ratio.push_back(std::vector<TH2D *>());
    eff_curves_1d.push_back(std::vector<TH1D *>());
    eff_curves_1d_ratio.push_back(std::vector<TH1D *>());

    for (int j = 0; j < FLAGS_centBins; ++j) {
      mc_2d[i].push_back((TH2D *)in_mc[i][j]->Project3D("YX"));
      match_2d[i].push_back((TH2D *)in_match[i][j]->Project3D("YX"));

      string eff_name =
          "efficiency_lumi_" + std::to_string(i) + "_cent_" + std::to_string(j);

      eff_curves[i].push_back((TH2D *)match_2d[i][j]->Clone(eff_name.c_str()));
      eff_curves[i][j]->Divide(mc_2d[i][j]);
      eff_curves_ratio[i].push_back((TH2D *)eff_curves[i][j]->Clone());

      string title = lumi_string[i] + " " + refcent_string[j] + " central";

      eff_curves[i][j]->SetTitle(title.c_str());
      eff_curves[i][j]->GetZaxis()->SetTitle("efficiency");
      eff_curves[i][j]->GetXaxis()->SetRangeUser(0.0, FLAGS_maxPt);
      eff_curves[i][j]->GetZaxis()->SetRangeUser(0.0, 1.05);
      eff_curves[i][j]->Draw("surf1");
      string full_eff_name = FLAGS_outdir + "/" + eff_name + ".pdf";
      c.SaveAs(full_eff_name.c_str());

      eff_curves_1d[i].push_back(match_2d[i][j]->ProjectionX());
      eff_curves_1d[i][j]->Divide(mc_2d[i][j]->ProjectionX());
      eff_curves_1d_ratio[i].push_back((TH1D *)eff_curves_1d[i][j]->Clone());
      if (i > 0) {
        eff_curves_ratio[i][j]->Divide(eff_curves_ratio[0][j]);
        eff_curves_1d_ratio[i][j]->Divide(eff_curves_1d_ratio[0][j]);

        full_eff_name = FLAGS_outdir + "/" + eff_name + "_ratio.pdf";
        eff_curves_ratio[i][j]->Draw("COLZ");
        c.SaveAs(full_eff_name.c_str());
      }
    }
  }

  dijetcore::histogramOpts hOpts;
  dijetcore::canvasOpts cOpts;
  cOpts.leg_upper_bound = 0.55;
  cOpts.leg_left_bound = 0.9;
  cOpts.leg_lower_bound = 0.20;
  cOpts.leg_right_bound = 0.5;

  for (int i = 0; i < eff_curves_1d.size(); ++i) {

    std::string canvas_name;
    std::string file_name;
    if (i == 0) {
      canvas_name = "low luminosity";
      file_name = "low_lumi_eff";
    } else if (i == 1) {
      canvas_name = "mid luminosity";
      file_name = "mid_lumi_eff";
    } else if (i == 2) {
      canvas_name = "high luminosity";
      file_name = "high_lumi_eff";
    } else {
      canvas_name = "Centrality";
      file_name = "centrality_eff";
    }
    eff_curves_1d[i][0]->GetYaxis()->SetRangeUser(0.0, 1.05);
    Overlay1D(eff_curves_1d[i], refcent_string, hOpts, cOpts, FLAGS_outdir,
              file_name, "", "p_{T}", "efficiency", canvas_name);

    if (i != 0) {
      if (i == 1) {
        canvas_name = "mid / low";
        file_name = "mid_low_ratio";
      } else if (i == 2) {
        canvas_name = "high / low";
        file_name = "high_low_ratio.pdf";
      }
      Overlay1D(eff_curves_1d_ratio[i], refcent_string, hOpts, cOpts,
                FLAGS_outdir, file_name, "", "p_{T}", "efficiency",
                canvas_name);
    }
  }

  // get errors
  std::vector<std::vector<TH2D *>> errors;
  for (int i = 0; i < FLAGS_lumiBins; ++i) {
    errors.push_back(std::vector<TH2D *>());
    for (int j = 0; j < FLAGS_centBins; ++j) {
      string err_name =
          "error_lumi_" + std::to_string(i) + "_cent_" + std::to_string(j);
      errors[i].push_back((TH2D *)eff_curves[i][j]->Clone(err_name.c_str()));

      for (int k = 0; k < (errors[i][j]->GetNbinsX() + 2) *
                              (errors[i][j]->GetNbinsY() + 1);
           ++k) {
        errors[i][j]->SetBinContent(k, errors[i][j]->GetBinError(k));
      }

      string title =
          lumi_string[i] + " " + refcent_string[j] + " central: error";

      err_name = FLAGS_outdir + "/" + err_name + ".pdf";
      errors[i][j]->SetTitle(title.c_str());
      errors[i][j]->GetZaxis()->SetRangeUser(0.0, 0.1);
      errors[i][j]->Draw("colz");
      c.SaveAs(err_name.c_str());
    }
  }

  // other histograms for QA
  TH3D *ptmatched = (TH3D *)in.Get("mcptvsmatchptvseta");
  ptmatched->GetXaxis()->SetRangeUser(0.0, 5.0);
  ptmatched->GetYaxis()->SetRangeUser(0.0, 6.0);
  ((TH2D *)ptmatched->Project3D("YX"))->Draw("colz");
  c.SetLogz();
  string matched_name = FLAGS_outdir + "/" + "match_pt.pdf";
  c.SaveAs(matched_name.c_str());

  TH2D *mcvsmatch = (TH2D *)in.Get("mcvsmatched");
  mcvsmatch->Draw("colz");
  string mcvsmatch_name = FLAGS_outdir + "/" + "mcvsmatch.pdf";
  c.SaveAs(mcvsmatch_name.c_str());

  TH2D *refzdc = (TH2D *)in.Get("refzdc");
  string refzdc_name = FLAGS_outdir + "/" + "refzdc.pdf";
  refzdc->Draw("colz");
  c.SaveAs(refzdc_name.c_str());
  c.SetLogz(false);

  TH1D *fit = (TH1D *)in.Get("fitpoints");
  string fit_name = FLAGS_outdir + "/" + "fitpoints.pdf";
  fit->Draw();
  c.SaveAs(fit_name.c_str());

  TH2D *dcapt = (TH2D *)in.Get("dcapt");
  TH1D *dca = dcapt->ProjectionX();
  string dca_name = FLAGS_outdir + "/" + "dca.pdf";
  dca->GetXaxis()->SetRangeUser(0.0, 3.0);
  dca->Draw();
  c.SetLogy();
  c.SaveAs(dca_name.c_str());
  c.SetLogy(false);
  // write to file
  out.cd();

  for (auto vec : eff_curves) {
    for (auto hist : vec) {
      hist->Write();
    }
  }

  out.Close();
  return 0;
}
