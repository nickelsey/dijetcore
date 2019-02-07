#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
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
  gStyle->SetOptTitle(0);
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
  TH2F *zdc_grefmult =
      (TH2F *)input_file.Get(MakeString(prefix, "zdcgref").c_str());
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
  copts.do_legend = false;
  dijetcore::canvasOpts coptslogz;
  coptslogz.log_z = true;
  coptslogz.do_legend = false;
  dijetcore::canvasOpts coptslogy;
  coptslogy.log_y = true;
  coptslogy.do_legend = false;
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
  LOG(INFO) << "Printing results for event QA";
  TProfile *zdc_vz_ave = zdc_vz->ProfileX();
  TProfile *vz_vx_ave = vz_vx->ProfileX();
  TProfile *vz_vy_ave = vz_vy->ProfileX();
  TProfile *zdc_ref_ave = zdc_refmult->ProfileX();
  TProfile *zdc_gref_ave = zdc_grefmult->ProfileX();
  TProfile *bbc_ref_ave = bbc_refmult->ProfileX();

  dijetcore::PrettyPrint1D(zdc_vz_ave, hopts, copts, FLAGS_legendLabel,
                           FLAGS_outputDir, "zdc_avg_vz", "", "zdcX [kHz]",
                           "<V_{z}> [cm]");
  dijetcore::PrettyPrint1D(vz_vx_ave, hopts, copts, FLAGS_legendLabel,
                           FLAGS_outputDir, "vz_ave_vx", "", "V_{z} [cm]",
                           "<V_{x}> [cm]");
  dijetcore::PrettyPrint1D(vz_vy_ave, hopts, copts, FLAGS_legendLabel,
                           FLAGS_outputDir, "vz_ave_vy", "", "V_{z} [cm]",
                           "<V_{y}> [cm]");
  dijetcore::PrettyPrint1D(zdc_ref_ave, hopts, copts, FLAGS_legendLabel,
                           FLAGS_outputDir, "zdc_ave_ref", "", "zdcX [kHz]",
                           "<refmult>");
  dijetcore::PrettyPrint1D(zdc_gref_ave, hopts, copts, FLAGS_legendLabel,
                           FLAGS_outputDir, "zdc_ave_gref", "", "zdcX [kHz]",
                           "<grefmult>");
  dijetcore::PrettyPrint1D(bbc_ref_ave, hopts, copts, FLAGS_legendLabel,
                           FLAGS_outputDir, "bbc_ave_ref", "", "bbcX [kHz]",
                           "<refmult>");

  // do runid QA if its present in the TFile
  if (run_id_refgref != nullptr) {
    LOG(INFO) << "Printing results for RunID QA";

    TProfile *runid_ref_ave =
        ((TH2D *)run_id_refgref->Project3D("YX"))->ProfileX();
    TProfile *runid_gref_ave =
        ((TH2D *)run_id_refgref->Project3D("ZX"))->ProfileX();
    TProfile *runid_nprim_ave =
        ((TH2D *)run_id_nprim_nglob->Project3D("YX"))->ProfileX();
    TProfile *runid_nglob_ave =
        ((TH2D *)run_id_nprim_nglob->Project3D("ZX"))->ProfileX();
    TProfile *runid_zdc_ave = run_id_zdc->ProfileX();
    TProfile *runid_bbc_ave = run_id_bbc->ProfileX();
    TProfile *runid_dvz_ave = run_id_vzvpdvz->ProfileX();
    TProfile *runid_vz_ave = run_id_vz->ProfileX();
    TProfile *runid_vx_ave =
        ((TH2D *)run_id_vx_vy->Project3D("YX"))->ProfileX();
    TProfile *runid_vy_ave =
        ((TH2D *)run_id_vx_vy->Project3D("ZX"))->ProfileX();

    dijetcore::PrettyPrint1D(runid_ref_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_ref", "", "Run ID",
                             "<refmult>");
    dijetcore::PrettyPrint1D(runid_gref_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_gref", "", "Run ID",
                             "<refmult>");
    dijetcore::PrettyPrint1D(runid_nprim_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_nprim", "", "Run ID",
                             "<N_{primary}>");
    dijetcore::PrettyPrint1D(runid_nglob_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_nglob", "", "Run ID",
                             "<N_{global}>");
    dijetcore::PrettyPrint1D(runid_zdc_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_zdc", "", "Run ID",
                             "<zdcX> [kHz]");
    dijetcore::PrettyPrint1D(runid_bbc_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_bbc", "", "Run ID",
                             "<bbcX> [kHz]");
    dijetcore::PrettyPrint1D(runid_dvz_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_dvz", "", "Run ID",
                             "<V_{z}-V_{z}^{vpd}> [cm]");
    dijetcore::PrettyPrint1D(runid_vz_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_vz", "", "Run ID",
                             "<V_{z}> [cm]");
    dijetcore::PrettyPrint1D(runid_vx_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_vx", "", "Run ID",
                             "<V_{x}> [cm]");
    dijetcore::PrettyPrint1D(runid_vy_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "runid_ave_vy", "", "Run ID",
                             "<V_{y}> [cm]");
  }

  // do track QA if its present in the TFile
  if (px_py != nullptr) {
    LOG(INFO) << "Printing results for track QA";

    dijetcore::Print2DSimple(px_py, hopts, coptslogz, FLAGS_outputDir,
                             "track_pxpy", "", "p_{X} [GeV/c]",
                             "p_{y} [GeV/c]");
    dijetcore::Print2DSimple(pz_px, hopts, coptslogz, FLAGS_outputDir,
                             "track_pzpx", "", "p_{Z} [GeV/c]",
                             "p_{X} [GeV/c]");
    dijetcore::Print2DSimple(pz_py, hopts, coptslogz, FLAGS_outputDir,
                             "track_pzpy", "", "p_{Z} [GeV/c]",
                             "p_{Y} [GeV/c]");
    dijetcore::Print2DSimple(eta_phi, hopts, copts, FLAGS_outputDir,
                             "track_etaphi", "", "#eta", "#phi");

    TProfile *zdc_px_ave = zdc_px->ProfileX();
    TProfile *zdc_py_ave = zdc_py->ProfileX();
    TProfile *zdc_pz_ave = zdc_pz->ProfileX();
    TProfile *zdc_pt_ave = zdc_pt->ProfileX();
    TProfile *zdc_dca_ave = zdc_dca->ProfileX();
    TProfile *zdc_nhit_ave = zdc_nhit->ProfileX();
    TProfile *zdc_nhitposs_ave = zdc_nhitposs->ProfileX();
    TProfile *zdc_nhitfrac_ave = zdc_nhitfrac->ProfileX();
    TProfile *zdc_eta_ave = zdc_eta->ProfileX();
    TProfile *zdc_phi_ave = zdc_phi->ProfileX();

    dijetcore::PrettyPrint1D(zdc_px_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_px", "", "zdcX [kHz]",
                             "<p_{X}> [GeV/c]");
    dijetcore::PrettyPrint1D(zdc_py_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_py", "", "zdcX [kHz]",
                             "<p_{Y}> [GeV/c]");
    dijetcore::PrettyPrint1D(zdc_pz_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_pz", "", "zdcX [kHz]",
                             "<p_{Z}> [GeV/c]");
    dijetcore::PrettyPrint1D(zdc_pt_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_pt", "", "zdcX [kHz]",
                             "<p_{T}> [GeV/c]");
    dijetcore::PrettyPrint1D(zdc_nhit_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_nhit", "", "zdcX [kHz]",
                             "<N_{hits fit}>");
    dijetcore::PrettyPrint1D(zdc_nhitposs_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_nhitposs", "",
                             "zdcX [kHz]", "<N_{hits possible}>");
    dijetcore::PrettyPrint1D(zdc_nhitfrac_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_nhitfrac", "",
                             "zdcX [kHz]", "<N_{hits fit}/N_{hits possible}>");
    dijetcore::PrettyPrint1D(zdc_eta_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_eta", "", "zdcX [kHz]",
                             "<#eta>");
    dijetcore::PrettyPrint1D(zdc_phi_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_phi", "", "zdcX [kHz]",
                             "<#phi>");

    if (run_id_track_pt != nullptr) {
      LOG(INFO) << "Printing results for run ID dependent track QA";

      TProfile *runid_pt_ave = run_id_track_pt->ProfileX();
      TProfile *runid_px_ave = run_id_track_px->ProfileX();
      TProfile *runid_py_ave = run_id_track_py->ProfileX();
      TProfile *runid_pz_ave = run_id_track_pz->ProfileX();
      TProfile *runid_eta_ave = run_id_track_eta->ProfileX();
      TProfile *runid_phi_ave = run_id_track_phi->ProfileX();
      TProfile *runid_dca_ave = run_id_track_dca->ProfileX();
      TProfile *runid_nhit_ave = run_id_track_nhit->ProfileX();
      TProfile *runid_nhitposs_ave = run_id_track_nhitposs->ProfileX();
      TProfile *runid_nhitfrac_ave = run_id_track_nhitfrac->ProfileX();

      dijetcore::PrettyPrint1D(runid_pt_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_pt", "", "Run ID",
                               "<p_{T}> [GeV/c]");
      dijetcore::PrettyPrint1D(runid_px_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_px", "", "Run ID",
                               "<p_{X}> [GeV/c]");
      dijetcore::PrettyPrint1D(runid_py_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_py", "", "Run ID",
                               "<p_{Y}> [GeV/c]");
      dijetcore::PrettyPrint1D(runid_pz_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_pz", "", "Run ID",
                               "<p_{Z}> [GeV/c]");
      dijetcore::PrettyPrint1D(runid_eta_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_eta", "", "Run ID",
                               "<#eta>");
      dijetcore::PrettyPrint1D(runid_phi_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_phi", "", "Run ID",
                               "<#phi>");
      dijetcore::PrettyPrint1D(runid_dca_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_dca", "", "Run ID",
                               "<DCA> [cm]");
      dijetcore::PrettyPrint1D(runid_nhit_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_nhit", "", "Run ID",
                               "<N_{hits fit}>");
      dijetcore::PrettyPrint1D(
          runid_nhitposs_ave, hopts, copts, FLAGS_legendLabel, FLAGS_outputDir,
          "runid_ave_nhitposs", "", "Run ID", "<N_{hits possible}>");
      dijetcore::PrettyPrint1D(runid_nhitfrac_ave, hopts, copts,
                               FLAGS_legendLabel, FLAGS_outputDir,
                               "runid_ave_nhitfrac", "", "Run ID",
                               "<N_{hits fit}/N_{hits possible}>");
    }
  }

  // do tower QA if its present in the TFile
  if (e_et != nullptr) {
    LOG(INFO) << "Printing results for track QA";

    dijetcore::Print2DSimple(tow_eta_phi, hopts, copts, FLAGS_outputDir,
                             "tower_etaphi", "", "#eta", "#phi");

    TH1D *e = e_et->ProjectionX();
    TH1D *et = e_et->ProjectionY();
    TProfile *zdc_e_ave = zdc_e->ProfileX();
    TProfile *zdc_et_ave = zdc_et->ProfileX();
    TProfile *zdc_adc_ave = zdc_adc->ProfileX();

    dijetcore::PrettyPrint1D(e, hopts, coptslogy, FLAGS_legendLabel,
                             FLAGS_outputDir, "tower_e_spectrum", "", "E [GeV]",
                             "counts");
    dijetcore::PrettyPrint1D(et, hopts, coptslogy, FLAGS_legendLabel,
                             FLAGS_outputDir, "tower_et_spectrum", "",
                             "E_{T} [GeV]", "counts");
    dijetcore::PrettyPrint1D(zdc_e_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_tow_e", "", "zdcX [kHz]",
                             "<E> [GeV]");
    dijetcore::PrettyPrint1D(zdc_et_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_tow_et", "",
                             "zdcX [kHz]", "<E_{T}> [GeV]");
    dijetcore::PrettyPrint1D(zdc_adc_ave, hopts, copts, FLAGS_legendLabel,
                             FLAGS_outputDir, "zdc_ave_tow_adc", "",
                             "zdcX [kHz]", "<ADC>");

    if (run_id_tower_e != nullptr) {
      LOG(INFO) << "Printing results for run ID dependent tower QA";
      TH2D *runid_e = run_id_tower_e->Projection(2, 0);
      TH2D *runid_et = run_id_tower_et->Projection(2, 0);
      TH2D *runid_adc = run_id_tower_adc->Projection(2, 0);
      TH2D *towerid_e = run_id_tower_e->Projection(2, 1);
      TH2D *towerid_et = run_id_tower_et->Projection(2, 1);
      TH2D *towerid_adc = run_id_tower_adc->Projection(2, 1);

      TH1D *towerid_activity = run_id_tower_e->Projection(1);

      dijetcore::Print2DSimple(runid_e, hopts, coptslogz, FLAGS_outputDir,
                               "runid_e", "", "Run ID", "E [GeV]");
      dijetcore::Print2DSimple(runid_et, hopts, coptslogz, FLAGS_outputDir,
                               "runid_et", "", "Run ID", "E_{T} [GeV]");
      dijetcore::Print2DSimple(runid_adc, hopts, coptslogz, FLAGS_outputDir,
                               "runid_adc", "", "Run ID", "ADC");
      dijetcore::Print2DSimple(towerid_e, hopts, coptslogz, FLAGS_outputDir,
                               "towerid_e", "", "Tower ID", "E [GeV]");
      dijetcore::Print2DSimple(towerid_et, hopts, coptslogz, FLAGS_outputDir,
                               "towerid_et", "", "Tower ID", "E_{T} [GeV]");
      dijetcore::Print2DSimple(towerid_adc, hopts, coptslogz, FLAGS_outputDir,
                               "towerid_adc", "", "TowerID", "ADC");
      dijetcore::PrettyPrint1D(towerid_activity, hopts, copts,
                               FLAGS_legendLabel, FLAGS_outputDir,
                               "towerid_activity", "", "Tower ID", "Counts");

      TProfile *runid_e_ave = runid_e->ProfileX();
      TProfile *runid_et_ave = runid_et->ProfileX();
      TProfile *runid_adc_ave = runid_adc->ProfileX();
      TProfile *towerid_e_ave = towerid_e->ProfileX();
      TProfile *towerid_et_ave = towerid_et->ProfileX();
      TProfile *towerid_adc_ave = towerid_adc->ProfileX();

      dijetcore::PrettyPrint1D(runid_e_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_e", "", "Run ID",
                               "<E> [GeV]");
      dijetcore::PrettyPrint1D(runid_et_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_et", "", "Run ID",
                               "<E_{T}> [GeV]");
      dijetcore::PrettyPrint1D(runid_adc_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "runid_ave_adc", "", "Run ID",
                               "<ADC> [GeV]");
      dijetcore::PrettyPrint1D(towerid_e_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "towerid_ave_e", "", "Tower ID",
                               "<E> [GeV]");
      dijetcore::PrettyPrint1D(towerid_et_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "towerid_ave_et", "",
                               "Tower ID", "<E_{T}> [GeV]");
      dijetcore::PrettyPrint1D(towerid_adc_ave, hopts, copts, FLAGS_legendLabel,
                               FLAGS_outputDir, "towerid_ave_adc", "",
                               "Tower ID", "<ADC> [GeV]");
    }
  }

  gflags::ShutDownCommandLineFlags();
  return 0;
}

// void PrettyPrint1D(H* h,
//                      histogramOpts hopts,
//                      canvasOpts copts,
//                      std::string hist_title,
//                      std::string output_loc,
//                      std::string output_name,
//                      std::string canvas_title,
//                      std::string x_axis_label,
//                      std::string y_axis_label,
//                      std::string legend_title = "")