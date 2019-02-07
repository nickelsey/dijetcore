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

DIJETCORE_DEFINE_string(input1, "",
                        "root file with results from run_qa for first dataset");
DIJETCORE_DEFINE_string(
    input2, "", "root file with results from run_qa for second dataset");
DIJETCORE_DEFINE_string(outputDir, "results", "directory for output");
DIJETCORE_DEFINE_string(histPrefix1, "p18ih", "histogram prefix for input1");
DIJETCORE_DEFINE_string(histPrefix2, "p17id", "histogram prefix for input2");
DIJETCORE_DEFINE_string(legendLabel1, "P18ih",
                        "label for histogram legends for input1");
DIJETCORE_DEFINE_string(legendLabel2, "P17id",
                        "label for histogram legends for input1");

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
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  // turn off print messages
  gErrorIgnoreLevel = kInfo + 1;

  // load input file
  if (FLAGS_input1.empty()) {
    LOG(ERROR) << "no input file specified for input1";
    return 1;
  }

  if (FLAGS_input2.empty()) {
    LOG(ERROR) << "no input file specified for input2";
    return 1;
  }

  // attempt to build directory for output files, if it does not already exist
  if (!dijetcore::CreateDirectory(FLAGS_outputDir)) {
    LOG(ERROR) << "Could not create output directory: " << FLAGS_outputDir;
    return 1;
  };

  unique_ptr<TFile> input_file1 = make_unique<TFile>(FLAGS_input1.c_str(), "READ");
  unique_ptr<TFile> input_file2 = make_unique<TFile>(FLAGS_input2.c_str(), "READ");

  dijetcore_map<string, unique_ptr<TFile>> files;
  files[FLAGS_histPrefix1] = std::move(input_file1);
  files[FLAGS_histPrefix2] = std::move(input_file2);

  for (auto &file : files) {
    if (!file.second->IsOpen()) {
      LOG(ERROR) << "Could not open input file: " << file.second->GetName()
                 << ", exiting";
    }
  }

  // parse the histogram prefixes and legend labels
  std::vector<string> prefixes{FLAGS_histPrefix1, FLAGS_histPrefix2};
  std::vector<string> legend_labels{FLAGS_legendLabel1, FLAGS_legendLabel2};

  // now we will zip them into a map with the prefix as key
  dijetcore_map<string, string> labels;
  for (int i = 0; i < prefixes.size(); ++i)
    labels[prefixes[i]] = legend_labels[i];

  // check to what data is present in the files (runid histograms,
  // track and tower info)
  bool do_runid_qa = true;
  bool do_track_qa = true;
  bool do_tower_qa = true;
  for (auto &prefix : prefixes) {
    // check to make sure the prefix is valid...
    if (files[prefix]->Get(MakeString(prefix, "zdcvz").c_str()) == nullptr) {
      LOG(ERROR) << "prefix: " << prefix << " not located in "
                 << files[prefix]->GetName();
      return 1;
    }

    // check if the input has runid, track and tower info
    if (files[prefix]->Get(MakeString(prefix, "runidvz").c_str()) == nullptr)
      do_runid_qa = false;
    if (files[prefix]->Get(MakeString(prefix, "pxpy").c_str()) == nullptr)
      do_track_qa = false;
    if (files[prefix]->Get(MakeString(prefix, "eet").c_str()) == nullptr)
      do_tower_qa = false;
  }

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

  // load the histograms
  for (auto &prefix : prefixes) {
    // load all default histograms
    zdc_vz[prefix] =
        (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcvz").c_str());
    vz_vx[prefix] =
        (TH2F *)files[prefix]->Get(MakeString(prefix, "vzvx").c_str());
    vz_vy[prefix] =
        (TH2F *)files[prefix]->Get(MakeString(prefix, "vzvy").c_str());
    zdc_refmult[prefix] =
        (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcref").c_str());
    zdc_grefmult[prefix] =
        (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcgref").c_str());
    bbc_refmult[prefix] =
        (TH2F *)files[prefix]->Get(MakeString(prefix, "bbcref").c_str());
    n_vertices[prefix] =
        (TH1F *)files[prefix]->Get(MakeString(prefix, "nvertex").c_str());

    if (do_runid_qa) {
      run_id_refgref[prefix] =
          (TH3F *)files[prefix]->Get(MakeString(prefix, "runidrefgref").c_str());
      run_id_nprim_nglob[prefix] = (TH3F *)files[prefix]->Get(
          MakeString(prefix, "runidnprimnglob").c_str());
      run_id_zdc[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "runidzdc").c_str());
      run_id_bbc[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "runidbbc").c_str());
      run_id_vzvpdvz[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "runiddvz").c_str());
      run_id_vz[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "runidvz").c_str());
      run_id_vx_vy[prefix] =
          (TH3F *)files[prefix]->Get(MakeString(prefix, "runidvxvy").c_str());
    }

    if (do_track_qa) {
      px_py[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "pxpy").c_str());
      pz_px[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "pzpx").c_str());
      pz_py[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "pzpy").c_str());
      zdc_px[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcpx").c_str());
      zdc_py[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcpy").c_str());
      zdc_pz[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcpz").c_str());
      zdc_pt[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcpt").c_str());
      zdc_dca[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcdca").c_str());
      zdc_nhit[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcnhit").c_str());
      zdc_nhitposs[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcnhitposs").c_str());
      zdc_nhitfrac[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcnhitfrac").c_str());
      zdc_eta[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdceta").c_str());
      zdc_phi[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcphi").c_str());
      eta_phi[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "tracketaphi").c_str());

      if (do_runid_qa) {
        run_id_track_pt[prefix] =
            (TH2F *)files[prefix]->Get(MakeString(prefix, "runidpt").c_str());
        run_id_track_px[prefix] =
            (TH2F *)files[prefix]->Get(MakeString(prefix, "runidpx").c_str());
        run_id_track_py[prefix] =
            (TH2F *)files[prefix]->Get(MakeString(prefix, "runidpy").c_str());
        run_id_track_pz[prefix] =
            (TH2F *)files[prefix]->Get(MakeString(prefix, "runidpz").c_str());
        run_id_track_eta[prefix] =
            (TH2F *)files[prefix]->Get(MakeString(prefix, "runideta").c_str());
        run_id_track_phi[prefix] =
            (TH2F *)files[prefix]->Get(MakeString(prefix, "runidphi").c_str());
        run_id_track_dca[prefix] =
            (TH2F *)files[prefix]->Get(MakeString(prefix, "runiddca").c_str());
        run_id_track_nhit[prefix] =
            (TH2F *)files[prefix]->Get(MakeString(prefix, "runidnhit").c_str());
        run_id_track_nhitposs[prefix] = (TH2F *)files[prefix]->Get(
            MakeString(prefix, "runidnhitposs").c_str());
        run_id_track_nhitfrac[prefix] = (TH2F *)files[prefix]->Get(
            MakeString(prefix, "runidnhitfrac").c_str());
      }
    }

    if (do_tower_qa) {
      e_et[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "eet").c_str());
      zdc_e[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdce").c_str());
      zdc_et[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcet").c_str());
      zdc_adc[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcadc").c_str());
      tow_eta_phi[prefix] =
          (TH2F *)files[prefix]->Get(MakeString(prefix, "towetaphi").c_str());

      if (do_runid_qa) {
        run_id_tower_e[prefix] = (THnSparseS *)files[prefix]->Get(
            MakeString(prefix, "runide").c_str());
        run_id_tower_et[prefix] = (THnSparseS *)files[prefix]->Get(
            MakeString(prefix, "runidet").c_str());
        run_id_tower_adc[prefix] = (THnSparseS *)files[prefix]->Get(
            MakeString(prefix, "runidadc").c_str());
      }
    }
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
  std::vector<TH1D *> refmult;
  std::vector<TH1D *> grefmult;
  std::vector<TH1D *> zdc;
  std::vector<TH1D *> vx;
  std::vector<TH1D *> vy;
  std::vector<TH1D *> vz;
  std::vector<TH1D *> dvz;
  for (auto &prefix : prefixes) {
    refmult.push_back(run_id_refgref[prefix]->ProjectionY());
    refmult.back()->Scale(1.0 / refmult.back()->Integral());
    grefmult.push_back(run_id_refgref[prefix]->ProjectionZ());
    grefmult.back()->Scale(1.0 / grefmult.back()->Integral());
    zdc.push_back(run_id_zdc[prefix]->ProjectionY());
    zdc.back()->Scale(1.0 / zdc.back()->Integral());
    vx.push_back(run_id_vx_vy[prefix]->ProjectionY());
    vx.back()->Scale(1.0 / vx.back()->Integral());
    vy.push_back(run_id_vx_vy[prefix]->ProjectionZ());
    vy.back()->Scale(1.0 / vy.back()->Integral());
    vz.push_back(run_id_vz[prefix]->ProjectionY());
    vz.back()->Scale(1.0 / vz.back()->Integral());
    dvz.push_back(run_id_vzvpdvz[prefix]->ProjectionY());
    dvz.back()->Scale(1.0 / dvz.back()->Integral());
  }

  // print comparisons and ratios between the two datasets

  google::ShutDownCommandLineFlags();
  return 0;
}

template<class H>
void PrintRatios(std::vector<H*> histograms, std::vector<string>& legend_entries,
                dijetcore::histogramOpts hopts, dijetcore::canvasOpts copts,
                string& loc, string& name, string& x_axis, string& y_axis, 
                string leg_title = "") {

};

// void PrintWithRatio(H* h1,
//                       H* h2,
//                       std::string h1_title,
//                       std::string h2_title,
//                       histogramOpts hopts,
//                       canvasOpts copts,
//                       std::string output_loc,
//                       std::string output_name,
//                       std::string canvas_title,
//                       std::string x_axis_label,
//                       std::string y_axis_label,
//                       std::string legend_title = "")