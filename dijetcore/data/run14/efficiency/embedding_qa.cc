// simple QA from the output of StEfficiencyAssessor
// (StEfficiencyAssessor is not part of the dijetcore repo, but can be found on
// my github ;)

#include "dijetcore/lib/assert.h"
#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/root/histogram_utils.h"
#include "dijetcore/util/root/root_print_utils.h"

#include <set>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TTree.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoReader.h"

DIJETCORE_DEFINE_string(
    input, "",
    "input root TFile (or txt file with TFile list) with embedding qa");
DIJETCORE_DEFINE_string(outputDir, "tmp", "output directory");
DIJETCORE_DEFINE_string(name, "embed_qa", "output filename");

template <class T> using vec2D = std::vector<std::vector<T>>;
template <class T>
using histContainer = std::vector<std::pair<vec2D<T>, std::string>>;

std::vector<TH1D *> ptByCent(TH3D *h) {
  std::vector<TH1D *> ret;
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    std::string name = dijetcore::MakeString(h->GetName(), i, "pt");
    ret.push_back((TH1D *)h->ProjectionY(name.c_str(), i, i));
    ret.back()->SetName(name.c_str());
    ret.back()->Scale(1.0 / ret.back()->Integral());
  }
  return ret;
}

std::vector<TH2D *> splitByCent(TH3D *h) {
  std::vector<TH2D *> ret;
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    std::string name = dijetcore::MakeString(h->GetName(), i);
    h->GetXaxis()->SetRange(i, i);
    ret.push_back((TH2D *)h->Project3D("zy"));
    ret.back()->SetName(name.c_str());
  }
  return ret;
}

vec2D<TH1D *> splitByCentAndPt(TH3D *h, std::string tag = "") {
  vec2D<TH1D *> ret;
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    ret.push_back(std::vector<TH1D *>());
    for (int j = 1; j <= h->GetYaxis()->GetNbins(); ++j) {
      std::string name = dijetcore::MakeString(h->GetName(), tag, i, "_", j);
      ret.back().push_back((TH1D *)h->ProjectionZ(name.c_str(), i, i, j, j));
    }
  }
  return ret;
}

vec2D<TH1D *> splitByCentAndPtNormalized(TH3D *h) {
  vec2D<TH1D *> ret = splitByCentAndPt(h);
  for (auto &cont : ret)
    for (auto &h : cont)
      if (h->Integral() > 0)
        h->Scale(1.0 / h->Integral());
  return ret;
}

template <class T>
void printComparisonsByBin(
    histContainer<T> &hists,
    vec2D<std::pair<std::string, std::string>> &grid_names,
    std::vector<unsigned> &cent_bins, std::vector<unsigned> &pt_bins,
    dijetcore::histogramOpts hopts, dijetcore::canvasOpts copts,
    std::string output_dir, std::string output_prefix, std::string x_label,
    std::string y_label, double y_min, double y_max) {
  for (int i = 0; i < cent_bins.size(); ++i) {
    for (int j = 0; j < pt_bins.size(); ++j) {
      int cent_bin = cent_bins[i];
      int pt_bin = pt_bins[j];
      std::vector<TH1D *> print_hists;
      std::vector<std::string> hist_names;
      for (auto &grid : hists) {
        print_hists.push_back(grid.first[cent_bin][pt_bin]);
        hist_names.push_back(grid.second);
      }

      std::string bin_name = dijetcore::MakeString(
          output_prefix, "_", grid_names[cent_bin][pt_bin].first, "_",
          grid_names[cent_bin][pt_bin].second);
      print_hists[0]->GetYaxis()->SetRangeUser(y_min, y_max);
      dijetcore::Overlay1D(print_hists, hist_names, hopts, copts, output_dir,
                           bin_name, "", x_label, y_label);
    }
  }
}

int main(int argc, char *argv[]) {
  // setup command line flags
  dijetcore::SetUsageMessage(
      "Generate plots of track/embedding QA for STAR embedding.");
  dijetcore::ParseCommandLineFlags(&argc, argv);

  // initialize logging - search for logging related command line flags
  dijetcore::InitLogging(argv[0]);

  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);

  // define centrality and pT bins
  std::vector<std::string> cent_bin_names{"0_5",    "5_10",  "10_20",
                                          "20_30",  "30_40", "40_50",
                                          "50_60%", "60_70", "70_80"};
  std::vector<std::string> pt_bin_names{"0_pt_1", "1_pt_2", "2_pt_3", "3_pt_4",
                                        "4_pt_5"};

  // select pt bins to print
  std::vector<unsigned> cent_bins{0, 1, 5};
  std::vector<unsigned> pt_bins{0, 3};

  // load input file
  TFile input(FLAGS_input.c_str(), "READ");

  // load histograms
  TH1D *vz = (TH1D *)input.Get("vz");
  TH1D *refmult = (TH1D *)input.Get("refmult");
  TH1D *grefmult = (TH1D *)input.Get("grefmult");
  TH1D *centrality = (TH1D *)input.Get("centrality");

  TH3D *mcrecotracks = (TH3D *)input.Get("mcrecotracks");
  TH2D *mctracks = (TH2D *)input.Get("mctracks");
  TH2D *recotracks = (TH2D *)input.Get("recotracks");

  TH3D *mceta = (TH3D *)input.Get("mceta");
  TH3D *mcphi = (TH3D *)input.Get("mcphi");

  TH3D *reconhit = (TH3D *)input.Get("reconhit");
  TH3D *recodca = (TH3D *)input.Get("recodca");
  TH3D *reconhitposs = (TH3D *)input.Get("reconhitposs");
  TH3D *recoeta = (TH3D *)input.Get("recoeta");
  TH3D *recophi = (TH3D *)input.Get("recophi");
  TH3D *recofitfrac = (TH3D *)input.Get("recofitfrac");
  TH3D *recodcascale = (TH3D *)input.Get("recodcascale");

  TH3D *reconhitcut = (TH3D *)input.Get("reconhitcut");
  TH3D *recodcacut = (TH3D *)input.Get("recodcacut");
  TH3D *reconhitposscut = (TH3D *)input.Get("reconhitposscut");
  TH3D *recoetacut = (TH3D *)input.Get("recoetacut");
  TH3D *recophicut = (TH3D *)input.Get("recophicut");
  TH3D *recofitfraccut = (TH3D *)input.Get("recocutfitfrac");

  TH3D *datanhit = (TH3D *)input.Get("datanhit");
  TH3D *datadca = (TH3D *)input.Get("datadca");
  TH3D *datanhitposs = (TH3D *)input.Get("datanhitposs");
  TH3D *dataeta = (TH3D *)input.Get("dataeta");
  TH3D *dataphi = (TH3D *)input.Get("dataphi");
  TH3D *datafitfrac = (TH3D *)input.Get("datafitfrac");
  TH3D *datadcascale = (TH3D *)input.Get("datadcascale");

  TH3D *recodcaext = (TH3D *)input.Get("recocutdcaext");
  TH3D *datadcaext = (TH3D *)input.Get("datadcaext");

  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  // create our histogram and canvas options
  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;
  dijetcore::canvasOpts coptsNoLeg;
  coptsNoLeg.do_legend = false;
  dijetcore::canvasOpts coptslogz;
  coptslogz.log_z = true;
  dijetcore::canvasOpts coptslogy;
  coptslogy.log_y = true;
  dijetcore::canvasOpts cOptsBottomLeg;
  cOptsBottomLeg.leg_upper_bound = 0.28;
  cOptsBottomLeg.leg_lower_bound = 0.20;
  cOptsBottomLeg.leg_right_bound = 0.89;
  cOptsBottomLeg.leg_left_bound = 0.7;
  dijetcore::canvasOpts cOptsBottomLeftLeg;
  cOptsBottomLeftLeg.leg_upper_bound = 0.4;
  cOptsBottomLeftLeg.leg_lower_bound = 0.18;
  cOptsBottomLeftLeg.leg_right_bound = 0.18;
  cOptsBottomLeftLeg.leg_left_bound = 0.4;
  dijetcore::canvasOpts cOptsBottomLeftLegLogy;
  cOptsBottomLeftLegLogy.log_y = true;
  cOptsBottomLeftLegLogy.leg_upper_bound = 0.4;
  cOptsBottomLeftLegLogy.leg_lower_bound = 0.18;
  cOptsBottomLeftLegLogy.leg_right_bound = 0.18;
  cOptsBottomLeftLegLogy.leg_left_bound = 0.4;
  dijetcore::canvasOpts cOptsTopLeftLeg;
  cOptsTopLeftLeg.leg_right_bound = 0.18;
  cOptsTopLeftLeg.leg_left_bound = 0.4;

  vec2D<std::pair<std::string, std::string>> grid_names;
  grid_names.resize(cent_bin_names.size());
  for (int i = 0; i < cent_bin_names.size(); ++i) {
    for (int j = 0; j < pt_bin_names.size(); ++j) {
      grid_names[i].push_back({cent_bin_names[i], pt_bin_names[j]});
    }
  }

  dijetcore::PrettyPrint1D(vz, hopts, copts, "embedding sample",
                           FLAGS_outputDir, "vz", "", "v_{z}", "counts");
  dijetcore::PrettyPrint1D(centrality, hopts, cOptsBottomLeg,
                           "embedding sample", FLAGS_outputDir, "centrality",
                           "", "centrality bin", "counts");
  dijetcore::PrettyPrint1D(refmult, hopts, coptslogy, "embedding sample",
                           FLAGS_outputDir, "refmult", "", "refmult", "counts");
  dijetcore::PrettyPrint1D(grefmult, hopts, coptslogy, "embedding sample",
                           FLAGS_outputDir, "grefmult", "", "grefmult",
                           "counts");

  // first do the reweighted DCA
  vec2D<TH1D *> dca_ext_reco = splitByCentAndPtNormalized(recodcaext);
  vec2D<TH1D *> dca_ext_data = splitByCentAndPtNormalized(datadcaext);

  vec2D<TH1D *> dca_ext_reco_resum;
  vec2D<TH1D *> dca_ext_data_resum;
  int reweight_pt_bins = 5;
  for (int i = 0; i < dca_ext_reco.size(); ++i) {
    dca_ext_reco_resum.push_back(
        std::vector<TH1D *>(reweight_pt_bins, nullptr));
    dca_ext_data_resum.push_back(
        std::vector<TH1D *>(reweight_pt_bins, nullptr));
    for (int j = 0; j < dca_ext_reco[i].size(); ++j) {
      unsigned bin = j / (dca_ext_reco[i].size() / reweight_pt_bins);
      if (dca_ext_reco_resum[i][bin] == nullptr) {
        std::string name_reco = dijetcore::MakeString("reco_resum_", bin);
        dca_ext_reco_resum[i][bin] =
            (TH1D *)dca_ext_reco[i][j]->Clone(name_reco.c_str());
        std::string name_data = dijetcore::MakeString("data_resum_", bin);
        dca_ext_data_resum[i][bin] =
            (TH1D *)dca_ext_data[i][j]->Clone(name_data.c_str());
      } else {
        dca_ext_reco_resum[i][bin]->Add(dca_ext_reco[i][j]);
        dca_ext_data_resum[i][bin]->Add(dca_ext_data[i][j]);
      }
    }
  }

  for (auto &c : dca_ext_reco_resum)
    for (auto &h : c) {
      h->Scale(1.0 / h->Integral());
    }
  for (auto &c : dca_ext_data_resum)
    for (auto &h : c) {
      h->Scale(1.0 / h->Integral());
    }

  histContainer<TH1D *> dca_ext_container{{dca_ext_reco_resum, "reco w/ cuts "},
                                          {dca_ext_data_resum, "muDst"}};

  printComparisonsByBin(dca_ext_container, grid_names, cent_bins, pt_bins,
                        hopts, coptslogy, FLAGS_outputDir, "dca_ext",
                        "DCA [cm]", "fraction", 0.0001, 0.5);

  // print out bins from the 3D histograms in slices of pT and centrality

  // printout pT
  std::vector<TH1D *> mc_pt = ptByCent(mceta);
  std::vector<TH1D *> reco_pt = ptByCent(recoeta);
  std::vector<TH1D *> reco_pt_cuts = ptByCent(recoetacut);
  std::vector<TH1D *> data_pt = ptByCent(dataeta);

  std::vector<TH1D *> cent_0_5_pt{mc_pt[0], reco_pt[0], reco_pt_cuts[0]};
  std::vector<TH1D *> cent_40_50_pt{mc_pt[5], reco_pt[5], reco_pt_cuts[5]};
  std::vector<std::string> data_names{"MC", "reco", "reco w/ cuts"};
  cent_0_5_pt.front()->GetYaxis()->SetRangeUser(0.01, 0.08);
  cent_40_50_pt.front()->GetYaxis()->SetRangeUser(0.01, 0.08);

  dijetcore::Overlay1D(cent_0_5_pt, data_names, hopts, copts, FLAGS_outputDir,
                       "pt_0_5", "", "p_{T} [GeV/c]", "fraction");
  dijetcore::Overlay1D(cent_40_50_pt, data_names, hopts, copts, FLAGS_outputDir,
                       "pt_40_50", "", "p_{T} [GeV/c]", "fraction");

  // get DCA as a function of pT
  vec2D<TH1D *> dca_reco = splitByCentAndPtNormalized(recodca);
  vec2D<TH1D *> dca_reco_cut = splitByCentAndPtNormalized(recodcacut);
  vec2D<TH1D *> dca_data = splitByCentAndPtNormalized(datadca);

  std::vector<unsigned> dca_pt_bins{1, 2, 3, 4};

  vec2D<std::pair<std::string, std::string>> dca_grid_names;
  dca_grid_names.resize(cent_bin_names.size());
  for (int i = 0; i < cent_bin_names.size(); ++i) {
    for (int j = 0; j < recodca->GetYaxis()->GetNbins(); ++j) {
      std::string pt_bin =
          dijetcore::MakeString(0.25 * j, "_pt_", 0.25 * (j + 1));
      dca_grid_names[i].push_back({cent_bin_names[i], pt_bin});
    }
  }

  histContainer<TH1D *> special_dca_container{
      {dca_reco, "reco"}, {dca_reco_cut, "reco w/ cuts"}, {dca_data, "muDst"}};

  printComparisonsByBin(special_dca_container, dca_grid_names, cent_bins,
                        dca_pt_bins, hopts, coptslogy, FLAGS_outputDir,
                        "lowdca", "DCA [cm]", "fraction", 0.0001, 0.5);

  // rebin before more general printouts of eta/phi/nhit/dca/etc
  mceta->RebinY(4);
  mceta->RebinZ(2);
  mcphi->RebinY(4);
  mcphi->RebinZ(2);
  recodca->RebinY(4);
  recodca->RebinZ(2);
  recoeta->RebinY(4);
  recoeta->RebinZ(2);
  recophi->RebinY(4);
  recophi->RebinZ(2);
  recodcascale->RebinY(4);
  recodcascale->RebinZ(2);
  recodcacut->RebinY(4);
  recodcacut->RebinZ(2);
  recoetacut->RebinY(4);
  recoetacut->RebinZ(2);
  recophicut->RebinY(4);
  recophicut->RebinZ(2);
  datadca->RebinY(4);
  datadca->RebinZ(2);
  dataeta->RebinY(4);
  dataeta->RebinZ(2);
  dataphi->RebinY(4);
  dataphi->RebinZ(2);
  datadcascale->RebinY(4);
  datadcascale->RebinZ(2);

  vec2D<TH1D *> mc_eta_split = splitByCentAndPtNormalized(mceta);
  vec2D<TH1D *> reco_eta_split = splitByCentAndPtNormalized(recoeta);
  vec2D<TH1D *> reco_cut_eta_split = splitByCentAndPtNormalized(recoetacut);
  vec2D<TH1D *> data_eta_split = splitByCentAndPtNormalized(dataeta);

  histContainer<TH1D *> eta_container{{mc_eta_split, "MC"},
                                      {reco_eta_split, "reco"},
                                      {reco_cut_eta_split, "reco w/ cuts"},
                                      {data_eta_split, "muDst"}};

  printComparisonsByBin(eta_container, grid_names, cent_bins, pt_bins, hopts,
                        copts, FLAGS_outputDir, "eta", "#eta", "fraction", 0.02,
                        0.06);

  // phi
  vec2D<TH1D *> mc_phi_split = splitByCentAndPtNormalized(mcphi);
  vec2D<TH1D *> reco_phi_split = splitByCentAndPtNormalized(recophi);
  vec2D<TH1D *> reco_cut_phi_split = splitByCentAndPtNormalized(recophicut);
  vec2D<TH1D *> data_phi_split = splitByCentAndPtNormalized(dataphi);

  histContainer<TH1D *> phi_container{{mc_phi_split, "MC"},
                                      {reco_phi_split, "reco"},
                                      {reco_cut_phi_split, "reco w/ cuts"},
                                      {data_phi_split, "muDst"}};

  printComparisonsByBin(phi_container, grid_names, cent_bins, pt_bins, hopts,
                        copts, FLAGS_outputDir, "phi", "#phi", "fraction", 0.02,
                        0.06);

  // nhit
  vec2D<TH1D *> reco_nhit_split = splitByCentAndPtNormalized(reconhit);
  vec2D<TH1D *> reco_cut_nhit_split = splitByCentAndPtNormalized(reconhitcut);
  vec2D<TH1D *> data_nhit_split = splitByCentAndPtNormalized(datanhit);

  histContainer<TH1D *> nhit_container{{reco_nhit_split, "reco"},
                                       {reco_cut_nhit_split, "reco w/ cuts"},
                                       {data_nhit_split, "muDst"}};

  printComparisonsByBin(nhit_container, grid_names, cent_bins, pt_bins, hopts,
                        cOptsTopLeftLeg, FLAGS_outputDir, "nhit", "N_{hits}",
                        "fraction", 0.00, 0.12);

  // nhitposs
  vec2D<TH1D *> reco_nhitposs_split = splitByCentAndPtNormalized(reconhitposs);
  vec2D<TH1D *> reco_cut_nhitposs_split =
      splitByCentAndPtNormalized(reconhitposscut);
  vec2D<TH1D *> data_nhitposs_split = splitByCentAndPtNormalized(datanhitposs);

  histContainer<TH1D *> nhitposs_container{
      {reco_nhitposs_split, "reco"},
      {reco_cut_nhitposs_split, "reco w/ cuts"},
      {data_nhitposs_split, "muDst"}};

  printComparisonsByBin(nhitposs_container, grid_names, cent_bins, pt_bins,
                        hopts, cOptsTopLeftLeg, FLAGS_outputDir, "nhitposs",
                        "N_{hits}^{possible}", "fraction", 0.00, 0.12);

  // fitfrac
  vec2D<TH1D *> reco_fitfrac_split = splitByCentAndPtNormalized(recofitfrac);
  vec2D<TH1D *> reco_cut_fitfrac_split =
      splitByCentAndPtNormalized(recofitfraccut);
  vec2D<TH1D *> data_fitfrac_split = splitByCentAndPtNormalized(datafitfrac);

  histContainer<TH1D *> fitfrac_container{
      {reco_fitfrac_split, "reco"},
      {reco_cut_fitfrac_split, "reco w/ cuts"},
      {data_fitfrac_split, "muDst"}};

  printComparisonsByBin(fitfrac_container, grid_names, cent_bins, pt_bins,
                        hopts, cOptsTopLeftLeg, FLAGS_outputDir, "fitfrac",
                        "N_{hit}/N_{hit}^{possible}", "fraction", 0.00, 0.15);

  // dca
  vec2D<TH1D *> reco_dca_split = splitByCentAndPtNormalized(recodca);
  vec2D<TH1D *> reco_cut_dca_split = splitByCentAndPtNormalized(recodcacut);
  vec2D<TH1D *> data_dca_split = splitByCentAndPtNormalized(datadca);

  histContainer<TH1D *> dca_container{{reco_dca_split, "reco"},
                                      {reco_cut_dca_split, "reco w/ cuts"},
                                      {data_dca_split, "muDst"}};

  printComparisonsByBin(dca_container, grid_names, cent_bins, pt_bins, hopts,
                        coptslogy, FLAGS_outputDir, "dca", "DCA [cm]",
                        "fraction", 0.0001, 0.5);

  google::ShutdownGoogleLogging();
  google::ShutDownCommandLineFlags();
  return 0;
}
