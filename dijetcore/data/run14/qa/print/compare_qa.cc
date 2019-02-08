#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/memory.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/root/histogram_utils.h"
#include "dijetcore/util/root/root_print_utils.h"

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

using dijetcore::Diff;
using dijetcore::dijetcore_map;
using dijetcore::make_unique;
using dijetcore::MakeString;
using dijetcore::Norm;
using dijetcore::RelativeDiff;
using dijetcore::unique_ptr;
using std::string;

template <class H>
TH1D *RelativeDiff(std::vector<H *> vec) {
  return RelativeDiff(vec[0], vec[1]);
}

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

  unique_ptr<TFile> input_file1 =
      make_unique<TFile>(FLAGS_input1.c_str(), "READ");
  unique_ptr<TFile> input_file2 =
      make_unique<TFile>(FLAGS_input2.c_str(), "READ");

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
  dijetcore_map<string, TH2F *> zdc_dvz;
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
    zdc_dvz[prefix] =
        (TH2F *)files[prefix]->Get(MakeString(prefix, "zdcdvz").c_str());
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
      run_id_refgref[prefix] = (TH3F *)files[prefix]->Get(
          MakeString(prefix, "runidrefgref").c_str());
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
  dijetcore::canvasOpts coptsNoLeg;
  coptsNoLeg.do_legend = false;
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

  std::vector<TH1D *> vx;
  std::vector<TH1D *> vy;
  std::vector<TH1D *> vz;
  std::vector<TH1D *> dvz;
  std::vector<TH1D *> zdc;
  std::vector<TH1F *> nvertices;
  std::vector<TH1D *> refmult;
  std::vector<TH1D *> grefmult;
  for (auto &prefix : prefixes) {
    vx.push_back(Norm(vz_vx[prefix]->ProjectionY()));
    vy.push_back(Norm(vz_vy[prefix]->ProjectionY()));
    vz.push_back(Norm(vz_vx[prefix]->ProjectionX()));
    dvz.push_back(Norm(zdc_dvz[prefix]->ProjectionY()));
    zdc.push_back(Norm(zdc_refmult[prefix]->ProjectionX()));
    nvertices.push_back(Norm((TH1F *)n_vertices[prefix]->Clone(
        MakeString(prefix, "nvert_clone").c_str())));
    refmult.push_back(Norm(zdc_refmult[prefix]->ProjectionY()));
    grefmult.push_back(Norm(zdc_grefmult[prefix]->ProjectionY()));
  }

  // print comparisons and ratios between the two datasets
  dijetcore::PrintWithRatio(vx, legend_labels, hopts, coptslogy,
                            FLAGS_outputDir, "vx", "V_{x} [cm]", "fraction");
  dijetcore::PrintWithRatio(vy, legend_labels, hopts, coptslogy,
                            FLAGS_outputDir, "vy", "V_{y} [cm]", "fraction");
  dijetcore::PrintWithRatio(vz, legend_labels, hopts, copts, FLAGS_outputDir,
                            "vz", "V_{z} [cm]", "fraction");
  dijetcore::PrintWithRatio(zdc, legend_labels, hopts, copts, FLAGS_outputDir,
                            "zdc", "zdcX [kHz]", "fraction");
  dijetcore::PrintWithRatio(nvertices, legend_labels, hopts, copts,
                            FLAGS_outputDir, "nvertices", "N_{vertices}",
                            "fraction");
  dijetcore::PrintWithRatio(refmult, legend_labels, hopts, coptslogy,
                            FLAGS_outputDir, "refmult", "refmult", "fraction");
  dijetcore::PrintWithRatio(grefmult, legend_labels, hopts, coptslogy,
                            FLAGS_outputDir, "grefmult", "grefmult",
                            "fraction");
  dijetcore::PrintWithRatio(dvz, legend_labels, hopts, copts, FLAGS_outputDir,
                            "dvz", "#DeltaV_{z}^{TPC, VPD} [cm]", "fraction");

  // start runid QA if info is present for both datasets
  if (do_runid_qa) {
    std::vector<TProfile *> runid_ref;
    std::vector<TProfile *> runid_gref;
    std::vector<TProfile *> runid_nprim;
    std::vector<TProfile *> runid_nglob;
    std::vector<TProfile *> runid_zdc;
    std::vector<TProfile *> runid_dvz;
    std::vector<TProfile *> runid_vz;
    std::vector<TProfile *> runid_vx;
    std::vector<TProfile *> runid_vy;

    for (auto &prefix : prefixes) {
      runid_ref.push_back(
          ((TH2D *)run_id_refgref[prefix]->Project3D("YX"))->ProfileX());
      runid_gref.push_back(
          ((TH2D *)run_id_refgref[prefix]->Project3D("ZX"))->ProfileX());
      runid_nprim.push_back(
          ((TH2D *)run_id_nprim_nglob[prefix]->Project3D("YX"))->ProfileX());
      runid_nglob.push_back(
          ((TH2D *)run_id_nprim_nglob[prefix]->Project3D("ZX"))->ProfileX());
      runid_zdc.push_back(run_id_zdc[prefix]->ProfileX());
      runid_dvz.push_back(run_id_vzvpdvz[prefix]->ProfileX());
      runid_vz.push_back(run_id_vz[prefix]->ProfileX());
      runid_vx.push_back(
          ((TH2D *)run_id_vx_vy[prefix]->Project3D("YX"))->ProfileX());
      runid_vy.push_back(
          ((TH2D *)run_id_vx_vy[prefix]->Project3D("ZX"))->ProfileX());
    }

    TH1D *runid_ref_diff = RelativeDiff(runid_ref);
    TH1D *runid_gref_diff = RelativeDiff(runid_gref);
    TH1D *runid_npart_diff = RelativeDiff(runid_nprim);
    TH1D *runid_nglob_diff = RelativeDiff(runid_nglob);
    TH1D *runid_zdc_diff = RelativeDiff(runid_zdc);
    TH1D *runid_dvz_diff = RelativeDiff(runid_dvz);
    TH1D *runid_vz_diff = RelativeDiff(runid_vz);
    TH1D *runid_vx_diff = RelativeDiff(runid_vx);
    TH1D *runid_vy_diff = RelativeDiff(runid_vy);

    dijetcore::PrintWithRatio(runid_ref, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_ref", "Run ID",
                              "<refmult>");
    dijetcore::PrettyPrint1D(runid_ref_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_ref", "", "Run ID",
                             "#Delta <refmult> / <refmult>");

    dijetcore::PrintWithRatio(runid_gref, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_gref", "Run ID",
                              "<grefmult>");
    dijetcore::PrettyPrint1D(runid_gref_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_gref", "", "Run ID",
                             "#Delta <grefmult> / <grefmult>");

    dijetcore::PrintWithRatio(runid_nprim, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_npart", "Run ID",
                              "<N_{Primary}>");
    dijetcore::PrettyPrint1D(runid_npart_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_nprim", "", "Run ID",
                             "#Delta <N_{Primary}> / <N_{Primary}>");

    dijetcore::PrintWithRatio(runid_nglob, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_nglob", "Run ID",
                              "<N_{Global}>");
    dijetcore::PrettyPrint1D(runid_nglob_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_nglob", "", "Run ID",
                             "#Delta <N_{Global}> / <N_{Global}>");

    dijetcore::PrintWithRatio(runid_zdc, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_zdc", "Run ID",
                              "<zdcX> [kHz]");
    dijetcore::PrettyPrint1D(runid_zdc_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_zdc", "", "Run ID",
                             "#Delta <zdc> / <zdc>");

    dijetcore::PrintWithRatio(runid_dvz, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_dvz", "Run ID",
                              "<#Delta V_{z}> [cm]");
    dijetcore::PrettyPrint1D(runid_dvz_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_dvz", "", "Run ID",
                             "#Delta <#Delta V_{z}> / <#Delta V_{z}>");

    dijetcore::PrintWithRatio(runid_vz, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_vz", "Run ID",
                              "<V_{z}> [cm]");
    dijetcore::PrettyPrint1D(runid_vz_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_vz", "", "Run ID",
                             "#Delta <V_{z}> / <V_{z}>");

    dijetcore::PrintWithRatio(runid_vx, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_vx", "Run ID",
                              "<V_{x}> [cm]");
    dijetcore::PrettyPrint1D(runid_vx_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_vx", "", "Run ID",
                             "#Delta <V_{x}> / <V_{x}>");

    dijetcore::PrintWithRatio(runid_vy, legend_labels, hopts, copts,
                              FLAGS_outputDir, "runid_vy", "Run ID",
                              "<V_{y}> [cm]");
    dijetcore::PrettyPrint1D(runid_vy_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "runid_rel_vy", "", "Run ID",
                             "#Delta <V_{y}> / <V_{y}>");
  }

  // start track QA if info is present in both datasets
  if (do_track_qa) {
    std::vector<TH1D *> track_px;
    std::vector<TH1D *> track_py;
    std::vector<TH1D *> track_pz;
    std::vector<TH1D *> track_pt;
    std::vector<TH1D *> track_dca;
    std::vector<TH1D *> track_nhit;
    std::vector<TH1D *> track_nhitposs;
    std::vector<TH1D *> track_nhitfrac;
    std::vector<TH1D *> track_eta;
    std::vector<TH1D *> track_phi;
    std::vector<TProfile *> zdc_px_ave;
    std::vector<TProfile *> zdc_py_ave;
    std::vector<TProfile *> zdc_pz_ave;
    std::vector<TProfile *> zdc_pt_ave;
    std::vector<TProfile *> zdc_dca_ave;
    std::vector<TProfile *> zdc_nhit_ave;
    std::vector<TProfile *> zdc_nhitposs_ave;
    std::vector<TProfile *> zdc_nhitfrac_ave;
    std::vector<TProfile *> zdc_eta_ave;
    std::vector<TProfile *> zdc_phi_ave;

    for (auto &prefix : prefixes) {
      track_px.push_back(Norm(px_py[prefix]->ProjectionX()));
      track_py.push_back(Norm(px_py[prefix]->ProjectionY()));
      track_pz.push_back(Norm(pz_px[prefix]->ProjectionX()));
      track_pt.push_back(Norm(zdc_pt[prefix]->ProjectionY()));
      track_dca.push_back(Norm(zdc_dca[prefix]->ProjectionY()));
      track_nhit.push_back(Norm(zdc_nhit[prefix]->ProjectionY()));
      track_nhitposs.push_back(Norm(zdc_nhitposs[prefix]->ProjectionY()));
      track_nhitfrac.push_back(Norm(zdc_nhitfrac[prefix]->ProjectionY()));
      track_eta.push_back(Norm(zdc_eta[prefix]->ProjectionY()));
      track_phi.push_back(Norm(zdc_phi[prefix]->ProjectionY()));
      zdc_px_ave.push_back(zdc_px[prefix]->ProfileX());
      zdc_py_ave.push_back(zdc_py[prefix]->ProfileX());
      zdc_pz_ave.push_back(zdc_pz[prefix]->ProfileX());
      zdc_pt_ave.push_back(zdc_pt[prefix]->ProfileX());
      zdc_dca_ave.push_back(zdc_dca[prefix]->ProfileX());
      zdc_nhit_ave.push_back(zdc_nhit[prefix]->ProfileX());
      zdc_nhitposs_ave.push_back(zdc_nhitposs[prefix]->ProfileX());
      zdc_nhitfrac_ave.push_back(zdc_nhitfrac[prefix]->ProfileX());
      zdc_eta_ave.push_back(zdc_eta[prefix]->ProfileX());
      zdc_phi_ave.push_back(zdc_phi[prefix]->ProfileX());
    }

    dijetcore::PrintWithRatio(track_px, legend_labels, hopts, coptslogy,
                              FLAGS_outputDir, "track_px", "p_{X} [GeV/c]",
                              "fraction");
    dijetcore::PrintWithRatio(track_py, legend_labels, hopts, coptslogy,
                              FLAGS_outputDir, "track_py", "p_{Y} [GeV/c]",
                              "fraction");
    dijetcore::PrintWithRatio(track_pz, legend_labels, hopts, coptslogy,
                              FLAGS_outputDir, "track_pz", "p_{z} [GeV/c]",
                              "fraction");
    dijetcore::PrintWithRatio(track_pt, legend_labels, hopts, coptslogy,
                              FLAGS_outputDir, "track_pt", "p_{T} [GeV/c]",
                              "fraction");
    dijetcore::PrintWithRatio(track_dca, legend_labels, hopts, copts,
                              FLAGS_outputDir, "track_dca", "DCA [cm]",
                              "fraction");
    dijetcore::PrintWithRatio(track_nhit, legend_labels, hopts, copts,
                              FLAGS_outputDir, "track_nhit", "N_{fitted hits}",
                              "fraction");
    dijetcore::PrintWithRatio(track_nhitposs, legend_labels, hopts, copts,
                              FLAGS_outputDir, "track_nhitposs",
                              "N_{possible hits}", "fraction");
    dijetcore::PrintWithRatio(track_nhitfrac, legend_labels, hopts, copts,
                              FLAGS_outputDir, "track_nhitfrac",
                              "N_{fitted hits}/N_{possible}", "fraction");
    dijetcore::PrintWithRatio(track_eta, legend_labels, hopts, copts,
                              FLAGS_outputDir, "track_eta", "#eta", "fraction");
    dijetcore::PrintWithRatio(track_phi, legend_labels, hopts, copts,
                              FLAGS_outputDir, "track_phi", "#phi", "fraction");

    TH1D *zdc_px_diff = RelativeDiff(zdc_px_ave);
    TH1D *zdc_py_diff = RelativeDiff(zdc_py_ave);
    TH1D *zdc_pz_diff = RelativeDiff(zdc_pz_ave);
    TH1D *zdc_pt_diff = RelativeDiff(zdc_pt_ave);
    TH1D *zdc_dca_diff = RelativeDiff(zdc_dca_ave);
    TH1D *zdc_nhit_diff = RelativeDiff(zdc_nhit_ave);
    TH1D *zdc_nhitposs_diff = RelativeDiff(zdc_nhitposs_ave);
    TH1D *zdc_nhitfrac_diff = RelativeDiff(zdc_nhitfrac_ave);
    TH1D *zdc_eta_diff = RelativeDiff(zdc_eta_ave);
    TH1D *zdc_phi_diff = RelativeDiff(zdc_phi_ave);

    dijetcore::PrintWithRatio(zdc_px_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_px", "zdcX [kHz]",
                              "<p_{X}> [GeV/c]");
    dijetcore::PrettyPrint1D(zdc_px_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_px", "", "zdcX [kHz]",
                             "#Delta <p_{X}> / <p_{X}>");

    dijetcore::PrintWithRatio(zdc_py_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_py", "zdcX [kHz]",
                              "<p_{Y}> [GeV/c]");
    dijetcore::PrettyPrint1D(zdc_py_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_py", "", "zdcY [kHz]",
                             "#Delta <p_{Y}> / <p_{Y}>");

    dijetcore::PrintWithRatio(zdc_pz_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_pz", "zdcX [kHz]",
                              "<p_{Z}> [GeV/c]");
    dijetcore::PrettyPrint1D(zdc_pz_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_pz", "", "zdcX [kHz]",
                             "#Delta <p_{Z}> / <p_{Z}>");

    dijetcore::PrintWithRatio(zdc_pt_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_pt", "zdcX [kHz]",
                              "<p_{T}> [GeV/c]");
    dijetcore::PrettyPrint1D(zdc_pt_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_pt", "", "zdcX [kHz]",
                             "#Delta <p_{T}> / <p_{T}>");

    dijetcore::PrintWithRatio(zdc_dca_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_dca", "zdcX [kHz]",
                              "<DCA> [cm]");
    dijetcore::PrettyPrint1D(zdc_dca_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_dca", "", "zdcX [kHz]",
                             "#Delta <DCA> / <DCA>");

    dijetcore::PrintWithRatio(zdc_nhit_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_nhit", "zdcX [kHz]",
                              "<N_{hits fit}>");
    dijetcore::PrettyPrint1D(zdc_nhit_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_nhit", "", "zdcX [kHz]",
                             "#Delta <N_{hits}> / <N_{hits}>");

    dijetcore::PrintWithRatio(zdc_nhitposs_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_nhitposs", "zdcX [kHz]",
                              "<N_{hits possible}>");
    dijetcore::PrettyPrint1D(zdc_nhitposs_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_nhitposs", "",
                             "zdcX [kHz]",
                             "#Delta <N_{hits poss}> / <N_{hits poss}>");

    dijetcore::PrintWithRatio(zdc_nhitfrac_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_nhitfrac", "zdcX [kHz]",
                              "<N_{hit frac}>");
    dijetcore::PrettyPrint1D(zdc_nhitfrac_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_nhitfrac", "",
                             "zdcX [kHz]",
                             "#Delta <N_{hits frac}> / <N_{hits frac}>");

    dijetcore::PrintWithRatio(zdc_eta_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_eta", "zdcX [kHz]",
                              "<#eta>");
    dijetcore::PrettyPrint1D(zdc_eta_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_eta", "", "zdcX [kHz]",
                             "#Delta <#eta> / <#eta>");

    dijetcore::PrintWithRatio(zdc_phi_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_ave_phi", "zdcX [kHz]",
                              "<#phi>");
    dijetcore::PrettyPrint1D(zdc_phi_diff, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_phi", "", "zdcX [kHz]",
                             "#Delta <#phi> / <#phi>");

    // do track qa by run if the information is present
    if (do_runid_qa) {
      std::vector<TProfile *> runid_pt;
      std::vector<TProfile *> runid_px;
      std::vector<TProfile *> runid_py;
      std::vector<TProfile *> runid_pz;
      std::vector<TProfile *> runid_dca;
      std::vector<TProfile *> runid_nhit;
      std::vector<TProfile *> runid_nhitposs;
      std::vector<TProfile *> runid_nhitfrac;
      std::vector<TProfile *> runid_eta;
      std::vector<TProfile *> runid_phi;

      for (auto &prefix : prefixes) {
        runid_pt.push_back(run_id_track_pt[prefix]->ProfileX());
        runid_px.push_back(run_id_track_px[prefix]->ProfileX());
        runid_py.push_back(run_id_track_py[prefix]->ProfileX());
        runid_pz.push_back(run_id_track_pz[prefix]->ProfileX());
        runid_dca.push_back(run_id_track_dca[prefix]->ProfileX());
        runid_nhit.push_back(run_id_track_nhit[prefix]->ProfileX());
        runid_nhitposs.push_back(run_id_track_nhitposs[prefix]->ProfileX());
        runid_nhitfrac.push_back(run_id_track_nhitfrac[prefix]->ProfileX());
        runid_eta.push_back(run_id_track_eta[prefix]->ProfileX());
        runid_phi.push_back(run_id_track_phi[prefix]->ProfileX());
      }

      TH1D *runid_pt_diff = RelativeDiff(runid_pt);
      TH1D *runid_px_diff = RelativeDiff(runid_px);
      TH1D *runid_py_diff = RelativeDiff(runid_py);
      TH1D *runid_pz_diff = RelativeDiff(runid_pz);
      TH1D *runid_dca_diff = RelativeDiff(runid_dca);
      TH1D *runid_nhit_diff = RelativeDiff(runid_nhit);
      TH1D *runid_nhitposs_diff = RelativeDiff(runid_nhitposs);
      TH1D *runid_nhitfrac_diff = RelativeDiff(runid_nhitfrac);
      TH1D *runid_eta_diff = RelativeDiff(runid_eta);
      TH1D *runid_phi_diff = RelativeDiff(runid_phi);

      dijetcore::PrintWithRatio(runid_pt, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_pt", "Run ID",
                                "<p_{T}> [GeV/c]");
      dijetcore::PrettyPrint1D(runid_pt_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_pt", "", "Run ID",
                               "#Delta <p_{T}> / <p_{T}>");

      dijetcore::PrintWithRatio(runid_px, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_px", "Run ID",
                                "<p_{X}> [GeV/c]");
      dijetcore::PrettyPrint1D(runid_px_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_px", "", "Run ID",
                               "#Delta <p_{X}> / <p_{X}>");

      dijetcore::PrintWithRatio(runid_py, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_py", "Run ID",
                                "<p_{Y}> [GeV/c]");
      dijetcore::PrettyPrint1D(runid_py_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_py", "", "Run ID",
                               "#Delta <p_{Y}> / <p_{Y}>");

      dijetcore::PrintWithRatio(runid_pz, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_pz", "Run ID",
                                "<p_{Z}> [GeV/c]");
      dijetcore::PrettyPrint1D(runid_pz_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_pz", "", "Run ID",
                               "#Delta <p_{Z}> / <p_{Z}>");

      dijetcore::PrintWithRatio(runid_dca, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_dca", "Run ID",
                                "<DCA> [cm]");
      dijetcore::PrettyPrint1D(runid_dca_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_dca", "", "Run ID",
                               "#Delta <DCA> / <DCA>");

      dijetcore::PrintWithRatio(runid_nhit, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_nhit", "Run ID",
                                "<N_{hit}>");
      dijetcore::PrettyPrint1D(runid_nhit_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_nhit", "", "Run ID",
                               "#Delta <N_{hit}> / <N_{hit}>");

      dijetcore::PrintWithRatio(runid_nhitposs, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_nhitposs", "Run ID",
                                "<N_{hitposs}>");
      dijetcore::PrettyPrint1D(runid_nhitposs_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_nhitposs", "",
                               "Run ID",
                               "#Delta <N_{hitposs}> / <N_{hitposs}>");

      dijetcore::PrintWithRatio(runid_nhitfrac, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_nhitfrac", "Run ID",
                                "<N_{hit fraction}>");
      dijetcore::PrettyPrint1D(
          runid_nhit_diff, hopts, coptsNoLeg, "", FLAGS_outputDir,
          "runid_rel_nhitfrac", "", "Run ID",
          "#Delta <N_{hit fraction}> / <N_{hit fraction}>");

      dijetcore::PrintWithRatio(runid_eta, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_eta", "Run ID",
                                "<#eta>");
      dijetcore::PrettyPrint1D(runid_eta_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_eta", "", "Run ID",
                               "#Delta <#eta> / <#eta>");

      dijetcore::PrintWithRatio(runid_phi, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_phi", "Run ID",
                                "<#phi>");
      dijetcore::PrettyPrint1D(runid_phi_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_phi", "", "Run ID",
                               "#Delta <#phi> / <#phi>");
    }
  }

  // start tower QA if info is present in both datasets
  if (do_tower_qa) {
    std::vector<TH1D *> tower_e;
    std::vector<TH1D *> tower_et;
    std::vector<TH1D *> tower_adc;
    std::vector<TH1D *> zdc_e_ave;
    std::vector<TH1D *> zdc_et_ave;
    std::vector<TH1D *> zdc_adc_ave;

    for (auto &prefix : prefixes) {
      tower_e.push_back(e_et[prefix]->ProjectionX());
      tower_et.push_back(e_et[prefix]->ProjectionY());
      tower_adc.push_back(zdc_adc[prefix]->ProjectionY());
      zdc_e_ave.push_back(zdc_e[prefix]->ProfileX());
      zdc_et_ave.push_back(zdc_et[prefix]->ProfileX());
      zdc_adc_ave.push_back(zdc_adc[prefix]->ProfileX());
    }

    TH1D *zdc_rel_e = RelativeDiff(zdc_e_ave);
    TH1D *zdc_rel_et = RelativeDiff(zdc_et_ave);
    TH1D *zdc_rel_adc = RelativeDiff(zdc_adc_ave);

    dijetcore::PrintWithRatio(tower_e, legend_labels, hopts, coptslogy,
                              FLAGS_outputDir, "tower_e", "E [GeV]",
                              "fraction");
    dijetcore::PrintWithRatio(tower_et, legend_labels, hopts, coptslogy,
                              FLAGS_outputDir, "tower_et", "E_{T} [GeV]",
                              "fraction");
    dijetcore::PrintWithRatio(tower_adc, legend_labels, hopts, coptslogy,
                              FLAGS_outputDir, "tower_adc", "ADC", "fraction");

    dijetcore::PrintWithRatio(zdc_e_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_e_ave", "zdcX [kHz]",
                              "<E> [GeV]");
    dijetcore::PrettyPrint1D(zdc_rel_e, hopts, coptsNoLeg, "", FLAGS_outputDir,
                             "zdc_rel_e", "", "zdcX [kHz]", "#Delta <E> / <E>");

    dijetcore::PrintWithRatio(zdc_et_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_et_ave", "zdcX [kHz]",
                              "<E_{T}> [GeV]");
    dijetcore::PrettyPrint1D(zdc_rel_et, hopts, coptsNoLeg, "", FLAGS_outputDir,
                             "zdc_rel_et", "", "zdcX [kHz]",
                             "#Delta <E_{T}> / <E_{T}>");

    dijetcore::PrintWithRatio(zdc_adc_ave, legend_labels, hopts, copts,
                              FLAGS_outputDir, "zdc_adc_ave", "zdcX [kHz]",
                              "<ADC>");
    dijetcore::PrettyPrint1D(zdc_rel_adc, hopts, coptsNoLeg, "",
                             FLAGS_outputDir, "zdc_rel_adc", "", "zdcX [kHz]",
                             "#Delta <ADC> / <ADC>");

    // do runid/towerid based QA if its present
    if (do_runid_qa) {
      //   dijetcore_map<string, THnSparseS *> run_id_tower_e;
      //   dijetcore_map<string, THnSparseS *> run_id_tower_et;
      //   dijetcore_map<string, THnSparseS *> run_id_tower_adc;
      std::vector<TProfile *> runid_e_ave;
      std::vector<TProfile *> runid_et_ave;
      std::vector<TProfile *> runid_adc_ave;
      std::vector<TProfile *> towerid_e_ave;
      std::vector<TProfile *> towerid_et_ave;
      std::vector<TProfile *> towerid_adc_ave;

      for (auto &prefix : prefixes) {
        runid_e_ave.push_back(
            ((TH2D *)run_id_tower_e[prefix]->Projection(2, 0))->ProfileX());
        runid_et_ave.push_back(
            ((TH2D *)run_id_tower_et[prefix]->Projection(2, 0))->ProfileX());
        runid_adc_ave.push_back(
            ((TH2D *)run_id_tower_adc[prefix]->Projection(2, 0))->ProfileX());
        towerid_e_ave.push_back(
            ((TH2D *)run_id_tower_e[prefix]->Projection(2, 1))->ProfileX());
        towerid_et_ave.push_back(
            ((TH2D *)run_id_tower_et[prefix]->Projection(2, 1))->ProfileX());
        towerid_adc_ave.push_back(
            ((TH2D *)run_id_tower_adc[prefix]->Projection(2, 1))->ProfileX());
      }

      TH1D *runid_e_diff = RelativeDiff(runid_e_ave);
      TH1D *runid_et_diff = RelativeDiff(runid_et_ave);
      TH1D *runid_adc_diff = RelativeDiff(runid_adc_ave);
      TH1D *towerid_e_diff = RelativeDiff(towerid_e_ave);
      TH1D *towerid_et_diff = RelativeDiff(towerid_et_ave);
      TH1D *towerid_adc_diff = RelativeDiff(towerid_adc_ave);

      dijetcore::PrintWithRatio(runid_e_ave, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_e_ave", "Run ID",
                                "<E> [GeV]");
      dijetcore::PrettyPrint1D(runid_e_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_e", "", "Run ID",
                               "#Delta <E> / <E>");

      dijetcore::PrintWithRatio(runid_et_ave, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_et_ave", "Run ID",
                                "<E_{T}> [GeV]");
      dijetcore::PrettyPrint1D(runid_et_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_et", "", "Run ID",
                               "#Delta <E_{T}> / <E_{T}>");

      dijetcore::PrintWithRatio(runid_adc_ave, legend_labels, hopts, copts,
                                FLAGS_outputDir, "runid_adc_ave", "Run ID",
                                "<ADC>");
      dijetcore::PrettyPrint1D(runid_adc_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "runid_rel_adc", "", "Run ID",
                               "#Delta <ADC> / <ADC>");

      dijetcore::PrintWithRatio(towerid_e_ave, legend_labels, hopts, copts,
                                FLAGS_outputDir, "towerid_e_ave", "Tower ID",
                                "<E> [GeV]");
      dijetcore::PrettyPrint1D(towerid_e_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "towerid_rel_e", "", "Tower ID",
                               "#Delta <E> / <E>");

      dijetcore::PrintWithRatio(towerid_et_ave, legend_labels, hopts, copts,
                                FLAGS_outputDir, "towerid_et_ave", "Tower ID",
                                "<E_{T}> [GeV]");
      dijetcore::PrettyPrint1D(towerid_et_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "towerid_rel_et", "", "Tower ID",
                               "#Delta <E_{T}> / <E_{T}>");

      dijetcore::PrintWithRatio(towerid_adc_ave, legend_labels, hopts, copts,
                                FLAGS_outputDir, "towerid_adc_ave", "Tower ID",
                                "<ADC>");
      dijetcore::PrettyPrint1D(towerid_adc_diff, hopts, coptsNoLeg, "",
                               FLAGS_outputDir, "towerid_rel_adc", "", "Tower ID",
                               "#Delta <ADC> / <ADC>");
    }
  }

  google::ShutDownCommandLineFlags();
  return 0;
}