#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/memory.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/root/root_print_routines.h"

#include <set>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TStyle.h"

DIJETCORE_DEFINE_string(input, "", "root file with results from run_qa");
DIJETCORE_DEFINE_string(outputDir, "results", "directory for output");
DIJETCORE_DEFINE_string(
    histPrefixes, "",
    "histogram prefixes separated by delim (set by -delim= flag)");
DIJETCORE_DEFINE_string(legendLabels, "P18ih",
                        "labels for histogram legends (in same order as "
                        "histPrefixes), separated by same delimiter");
DIJETCORE_DEFINE_string(
    delim, ",",
    "delimiter for histPrefixes and legendLabels, must be single character");

using dijetcore::dijetcore_map;
using dijetcore::make_unique;
using dijetcore::MakeString;
using dijetcore::unique_ptr;
using std::string;

int main(int argc, char *argv[]) {
  string usage = "Run 14 QA print routines";

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
  gErrorIgnoreLevel = kInfo + 1;

  // load input file
  if (FLAGS_input.empty()) {
    LOG(ERROR) << "no input file specified";
    return 1;
  }

  // attempt to build directory for output files, if it does not already exist
  if (!dijetcore::CreateDirectory(FLAGS_outputDir)) {
    LOG(ERROR) << "Could not create output directory: " << FLAGS_outputDir;
    return 1;
  };

  TFile input_file(FLAGS_input.c_str(), "READ");

  if (!input_file.IsOpen()) {
    LOG(ERROR) << "Could not open input file: " << FLAGS_input << ", exiting";
    return 1;
  }

  // parse the histogram prefixes and legend labels
  std::vector<string> prefixes;
  std::vector<string> legend_labels;
  if (FLAGS_delim.size() > 1) {
    LOG(ERROR) << "Can not use a delimiter of size greater that 1";
    return 1;
  }
  char delim = char(FLAGS_delim[0]);
  dijetcore::SplitString<std::vector<string>>(FLAGS_histPrefixes, prefixes,
                                              delim);
  dijetcore::SplitString<std::vector<string>>(FLAGS_legendLabels, legend_labels,
                                              delim);

  if (prefixes.size() != legend_labels.size()) {
    LOG(ERROR)
        << "size of prefixes does not equal size of legend_labels; exiting";
    return 1;
  }

  // now we will zip them into a map with the prefix as key
  dijetcore_map<string, string> labels;
  for (int i = 0; i < prefixes.size(); ++i)
    labels[prefixes[i]] = legend_labels[i];

  // histograms will be stored in a map using prefix as the key
  dijetcore_map<string, TH3F *> run_id_refgref;
  dijetcore_map<string, TH3F *> run_id_nprim_nglob;
  dijetcore_map<string, TH2F *> run_id_zdc;
  dijetcore_map<string, TH2F *> run_id_bbc;
  dijetcore_map<string, TH2F *> run_id_vzvpdvz;
  dijetcore_map<string, TH2F *> run_id_vz;
  dijetcore_map<string, TH3F *> run_id_vx_vy;
  dijetcore_map<string, TH2F *> run_id_track_pt;
  dijetcore_map<string, TH2F *> run_id_track_px;
  dijetcore_map<string, TH2F *> run_id_track_py;
  dijetcore_map<string, TH2F *> run_id_track_pz;
  dijetcore_map<string, TH2F *> run_id_track_eta;
  dijetcore_map<string, TH2F *> run_id_track_phi;
  dijetcore_map<string, TH2F *> run_id_track_dca;
  dijetcore_map<string, TH2F *> run_id_track_nhit;
  dijetcore_map<string, TH2F *> run_id_track_nhitposs;
  dijetcore_map<string, TH2F *> run_id_track_nhitfrac;
  dijetcore_map<string, THnSparseS *> run_id_tower_e;
  dijetcore_map<string, THnSparseS *> run_id_tower_et;
  dijetcore_map<string, THnSparseS *> run_id_tower_adc;
  dijetcore_map<string, TH2F *> zdc_vz;
  dijetcore_map<string, TH2F *> vz_vx;
  dijetcore_map<string, TH2F *> vz_vy;
  dijetcore_map<string, TH2F *> zdc_refmult;
  dijetcore_map<string, TH2F *> zdc_grefmult;
  dijetcore_map<string, TH2F *> bbc_refmult;
  dijetcore_map<string, TH1F *> n_vertices;
  dijetcore_map<string, TH2F *> px_py;
  dijetcore_map<string, TH2F *> pz_px;
  dijetcore_map<string, TH2F *> pz_py;
  dijetcore_map<string, TH2F *> zdc_px;
  dijetcore_map<string, TH2F *> zdc_py;
  dijetcore_map<string, TH2F *> zdc_pz;
  dijetcore_map<string, TH2F *> zdc_pt;
  dijetcore_map<string, TH2F *> zdc_dca;
  dijetcore_map<string, TH2F *> zdc_nhit;
  dijetcore_map<string, TH2F *> zdc_nhitposs;
  dijetcore_map<string, TH2F *> zdc_nhitfrac;
  dijetcore_map<string, TH2F *> zdc_eta;
  dijetcore_map<string, TH2F *> zdc_phi;
  dijetcore_map<string, TH2F *> eta_phi;
  dijetcore_map<string, TH2F *> e_et;
  dijetcore_map<string, TH2F *> zdc_e;
  dijetcore_map<string, TH2F *> zdc_et;
  dijetcore_map<string, TH2F *> zdc_adc;
  dijetcore_map<string, TH2F *> tow_eta_phi;

  // first, check to make sure all prefixes are valid, if not, remove them
  prefixes.erase(
      std::remove_if(prefixes.begin(), prefixes.end(),
                     [&input_file](auto &a) {
                       return input_file.Get(MakeString(a, "zdcvz").c_str())
                                  ? false
                                  : true;
                     }),
      prefixes.end());

  // load all histograms, and we will keep track of if all labels have runid
  // info, track info and tower info
  bool do_runid = true;
  bool do_track = true;
  bool do_tower = true;
  for (auto &prefix : prefixes) {
    if (input_file.Get(MakeString(prefix, "runidrefgref").c_str()) == nullptr)
      do_runid = false;
    if (input_file.Get(MakeString(prefix, "pxpy").c_str()) == nullptr)
      do_track = false;
    if (input_file.Get(MakeString(prefix, "eet").c_str()) == nullptr)
      do_tower = false;

    run_id_refgref[prefix] =
        (TH3F *)input_file.Get(MakeString(prefix, "runidrefgref").c_str());
    run_id_nprim_nglob[prefix] =
        (TH3F *)input_file.Get(MakeString(prefix, "runidnprimnglob").c_str());
    run_id_zdc[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidzdc").c_str());
    run_id_bbc[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidbbc").c_str());
    run_id_vzvpdvz[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runiddvz").c_str());
    run_id_vz[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidvz").c_str());
    run_id_vx_vy[prefix] =
        (TH3F *)input_file.Get(MakeString(prefix, "runidvxvy").c_str());
    run_id_track_pt[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidpt").c_str());
    run_id_track_px[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidpx").c_str());
    run_id_track_py[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidpy").c_str());
    run_id_track_pz[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidpz").c_str());
    run_id_track_eta[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runideta").c_str());
    run_id_track_phi[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidphi").c_str());
    run_id_track_dca[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runiddca").c_str());
    run_id_track_nhit[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidnhit").c_str());
    run_id_track_nhitposs[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidnhitposs").c_str());
    run_id_track_nhitfrac[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "runidnhitfrac").c_str());
    run_id_tower_e[prefix] =
        (THnSparseS *)input_file.Get(MakeString(prefix, "runide").c_str());
    run_id_tower_et[prefix] =
        (THnSparseS *)input_file.Get(MakeString(prefix, "runidet").c_str());
    run_id_tower_adc[prefix] =
        (THnSparseS *)input_file.Get(MakeString(prefix, "runidadc").c_str());
    zdc_vz[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcvz").c_str());
    vz_vx[prefix] = (TH2F *)input_file.Get(MakeString(prefix, "vzvx").c_str());
    vz_vy[prefix] = (TH2F *)input_file.Get(MakeString(prefix, "vzvy").c_str());
    zdc_refmult[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcref").c_str());
    zdc_grefmult[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcgref").c_str());
    bbc_refmult[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "bbcref").c_str());
    n_vertices[prefix] =
        (TH1F *)input_file.Get(MakeString(prefix, "nvertex").c_str());
    px_py[prefix] = (TH2F *)input_file.Get(MakeString(prefix, "pxpy").c_str());
    pz_px[prefix] = (TH2F *)input_file.Get(MakeString(prefix, "pzpx").c_str());
    pz_py[prefix] = (TH2F *)input_file.Get(MakeString(prefix, "pzpy").c_str());
    zdc_px[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcpx").c_str());
    zdc_py[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcpy").c_str());
    zdc_pz[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcpz").c_str());
    zdc_pt[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcpt").c_str());
    zdc_dca[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcdca").c_str());
    zdc_nhit[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcnhit").c_str());
    zdc_nhitposs[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcnhitposs").c_str());
    zdc_nhitfrac[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcnhitfrac").c_str());
    zdc_eta[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdceta").c_str());
    zdc_phi[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcphi").c_str());
    eta_phi[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "tracketaphi").c_str());
    e_et[prefix] = (TH2F *)input_file.Get(MakeString(prefix, "eet").c_str());
    zdc_e[prefix] = (TH2F *)input_file.Get(MakeString(prefix, "zdce").c_str());
    zdc_et[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcet").c_str());
    zdc_adc[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "zdcadc").c_str());
    tow_eta_phi[prefix] =
        (TH2F *)input_file.Get(MakeString(prefix, "towetaphi").c_str());
  }

  // create our histogram and canvas options
  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;
  dijetcore::canvasOpts coptslogz;
  coptslogz.log_z = true;
  dijetcore::canvasOpts coptslogy;
  coptslogy.log_y = true;
  dijetcore::canvasOpts cOptsBottomLeg;
  cOptsBottomLeg.leg_upper_bound = 0.4;
  cOptsBottomLeg.leg_lower_bound = 0.18;
  cOptsBottomLeg.leg_right_bound = 0.9;
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

  // ***************************************************
  // event level QA
  std::vector<unique_ptr<TH1D>> refmult;
  std::vector<unique_ptr<TH1D>> grefmult;
  std::vector<unique_ptr<TH1D>> zdc;
  std::vector<unique_ptr<TH1D>> vx;
  std::vector<unique_ptr<TH1D>> vy;
  std::vector<unique_ptr<TH1D>> vz;
  std::vector<unique_ptr<TH1D>> dvz;
  for (auto &prefix : prefixes) {
    refmult.push_back(
        make_unique<TH1D>(*run_id_refgref[prefix]->ProjectionY()));
    // refmult.back()->Scale(1.0 / refmult.back()->);
    grefmult.push_back(
        make_unique<TH1D>(*run_id_refgref[prefix]->ProjectionZ()));
    zdc.push_back(make_unique<TH1D>(*run_id_zdc[prefix]->ProjectionY()));
    vx.push_back(make_unique<TH1D>(*run_id_vx_vy[prefix]->ProjectionY()));
    vy.push_back(make_unique<TH1D>(*run_id_vx_vy[prefix]->ProjectionZ()));
    vz.push_back(make_unique<TH1D>(*run_id_vz[prefix]->ProjectionY()));
    dvz.push_back(make_unique<TH1D>(*run_id_vzvpdvz[prefix]->ProjectionY()));
  }

  google::ShutDownCommandLineFlags();
  return 0;
}