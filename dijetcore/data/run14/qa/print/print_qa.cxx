#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/filesystem.h"
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
DIJETCORE_DEFINE_string(histPrefix, "", "histogram prefix");
DIJETCORE_DEFINE_string(legendLabel, "P18ih", "label for histogram legends");

using dijetcore::MakeString;
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
    LOG(FATAL) << "Could not create output directory: " << FLAGS_outputDir;
  };

  TFile input_file(FLAGS_input.c_str(), "READ");

  if (!input_file.IsOpen()) {
    LOG(ERROR) << "Could not open input file: " << FLAGS_input << ", exiting";
    return 1;
  }

  // load histograms, account for any histogram name prefix
  string prefix = FLAGS_histPrefix;

  TH3F *run_id_refgref =
      (TH3F *)input_file.Get(MakeString(prefix, "runidrefgref").c_str());
  TH3F *run_id_nprim_nglob =
      (TH3F *)input_file.Get(MakeString(prefix, "runidnprimnglob").c_str());
  TH2F *run_id_zdc =
      (TH2F *)input_file.Get(MakeString(prefix, "runidzdc").c_str());
  TH2F *run_id_bbc =
      (TH2F *)input_file.Get(MakeString(prefix, "runidbbc").c_str());
  TH2F *run_id_vzvpdvz =
      (TH2F *)input_file.Get(MakeString(prefix, "runiddvz").c_str());
  TH2F *run_id_vz =
      (TH2F *)input_file.Get(MakeString(prefix, "runidvz").c_str());
  TH3F *run_id_vx_vy =
      (TH3F *)input_file.Get(MakeString(prefix, "runidvxvy").c_str());
  TH2F *run_id_track_pt =
      (TH2F *)input_file.Get(MakeString(prefix, "runidpt").c_str());
  TH2F *run_id_track_px =
      (TH2F *)input_file.Get(MakeString(prefix, "runidpx").c_str());
  TH2F *run_id_track_py =
      (TH2F *)input_file.Get(MakeString(prefix, "runidpy").c_str());
  TH2F *run_id_track_pz =
      (TH2F *)input_file.Get(MakeString(prefix, "runidpz").c_str());
  TH2F *run_id_track_eta =
      (TH2F *)input_file.Get(MakeString(prefix, "runideta").c_str());
  TH2F *run_id_track_phi =
      (TH2F *)input_file.Get(MakeString(prefix, "runidphi").c_str());
  TH2F *run_id_track_dca =
      (TH2F *)input_file.Get(MakeString(prefix, "runiddca").c_str());
  TH2F *run_id_track_nhit =
      (TH2F *)input_file.Get(MakeString(prefix, "runidnhit").c_str());
  TH2F *run_id_track_nhitposs =
      (TH2F *)input_file.Get(MakeString(prefix, "runidnhitposs").c_str());
  TH2F *run_id_track_nhitfrac =
      (TH2F *)input_file.Get(MakeString(prefix, "runidnhitfrac").c_str());
  THnSparseS *run_id_tower_e =
      (THnSparseS *)input_file.Get(MakeString(prefix, "runide").c_str());
  THnSparseS *run_id_tower_et =
      (THnSparseS *)input_file.Get(MakeString(prefix, "runidet").c_str());
  THnSparseS *run_id_tower_adc =
      (THnSparseS *)input_file.Get(MakeString(prefix, "runidadc").c_str());
  TH2F *zdc_vz = (TH2F *)input_file.Get(MakeString(prefix, "zdcvz").c_str());
  TH2F *vz_vx = (TH2F *)input_file.Get(MakeString(prefix, "vzvx").c_str());
  TH2F *vz_vy = (TH2F *)input_file.Get(MakeString(prefix, "vzvy").c_str());
  TH2F *zdc_refmult =
      (TH2F *)input_file.Get(MakeString(prefix, "zdcref").c_str());
  TH2F *bbc_refmult =
      (TH2F *)input_file.Get(MakeString(prefix, "bbcref").c_str());
  TH1F *n_vertices =
      (TH1F *)input_file.Get(MakeString(prefix, "nvertex").c_str());
  TH2F *px_py = (TH2F *)input_file.Get(MakeString(prefix, "pxpy").c_str());
  TH2F *pz_px = (TH2F *)input_file.Get(MakeString(prefix, "pzpx").c_str());
  TH2F *pz_py = (TH2F *)input_file.Get(MakeString(prefix, "pzpy").c_str());
  TH2F *zdc_px = (TH2F *)input_file.Get(MakeString(prefix, "zdcpx").c_str());
  TH2F *zdc_py = (TH2F *)input_file.Get(MakeString(prefix, "zdcpy").c_str());
  TH2F *zdc_pz = (TH2F *)input_file.Get(MakeString(prefix, "zdcpz").c_str());
  TH2F *zdc_pt = (TH2F *)input_file.Get(MakeString(prefix, "zdcpt").c_str());
  TH2F *zdc_dca = (TH2F *)input_file.Get(MakeString(prefix, "zdcdca").c_str());
  TH2F *zdc_nhit =
      (TH2F *)input_file.Get(MakeString(prefix, "zdcnhit").c_str());
  TH2F *zdc_nhitposs =
      (TH2F *)input_file.Get(MakeString(prefix, "zdcnhitposs").c_str());
  TH2F *zdc_nhitfrac =
      (TH2F *)input_file.Get(MakeString(prefix, "zdcnhitfrac").c_str());
  TH2F *zdc_eta = (TH2F *)input_file.Get(MakeString(prefix, "zdceta").c_str());
  TH2F *zdc_phi = (TH2F *)input_file.Get(MakeString(prefix, "zdcphi").c_str());
  TH2F *eta_phi =
      (TH2F *)input_file.Get(MakeString(prefix, "tracketaphi").c_str());
  TH2F *e_et = (TH2F *)input_file.Get(MakeString(prefix, "eet").c_str());
  TH2F *zdc_e = (TH2F *)input_file.Get(MakeString(prefix, "zdce").c_str());
  TH2F *zdc_et = (TH2F *)input_file.Get(MakeString(prefix, "zdcet").c_str());
  TH2F *zdc_adc = (TH2F *)input_file.Get(MakeString(prefix, "zdcadc").c_str());
  TH2F *tow_eta_phi =
      (TH2F *)input_file.Get(MakeString(prefix, "towetaphi").c_str());

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

  // start by printing the general, event level QA

  TProfile* zdc_vz_ave = zdc_vz->ProfileX();
  TProfile* vz_vx_ave = vz_vx->ProfileX();
  TProfile* vz_vy_ave = vz_vy->ProfileX();
  TProfile* zdc_ref_ave = zdc_refmult->ProfileX();

  // unique_ptr<TH2F> zdc_vz_;
  // unique_ptr<TH2F> vz_vx_;
  // unique_ptr<TH2F> vz_vy_;
  // unique_ptr<TH2F> zdc_refmult_;
  // unique_ptr<TH2F> bbc_refmult_;
  // unique_ptr<TH1F> n_vertices_;


  gflags::ShutDownCommandLineFlags();
  return 0;
}
