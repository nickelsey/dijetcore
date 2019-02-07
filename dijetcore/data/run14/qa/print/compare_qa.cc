#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/memory.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/root/root_print_utils.h"
#include "dijetcore/util/root/histogram_utils.h"

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

template<class H>
H* Norm(H* h) {
  h->Scale(1.0 / h->Integral());
  return h;
}

template<class H1, class H2>
TH1D* Diff(H1* h1, H2* h2) {
  if (h1->GetNbins() != h2->GetNbins() ||
      h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax() ||
      h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin()) {
    LOG(ERROR) << "can't subtract histograms with different bins";
  }
  string name = MakeString(h1->GetName(), "minus", h2->GetName());
  TH1D* diff = new TH1D(name.c_str(), "", h1->GetNbins(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());

  for (int i = 1; i <= h1->GetNbins(); ++i) {
    double binContent = h1->GetBinContent(i) - h2->GetBinContent(i);
    double error1 = h1->GetBinError(i);
    double error2 = h2->GetBinError(i);
    double binError = sqrt(pow(error1, 2.0) + pow(error2, 2.0));

    diff->SetBinContent(i, binContent);
    diff->SetBinError(i, binError);
  }
  return diff;
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
    //dvz.push_back(Norm(zdc_dvz[prefix]->ProjectionY()));
    zdc.push_back(Norm(zdc_refmult[prefix]->ProjectionX()));
    nvertices.push_back(Norm((TH1F*)n_vertices[prefix]->Clone(MakeString(prefix, "nvert_clone").c_str())));
    refmult.push_back(Norm(zdc_refmult[prefix]->ProjectionY()));
    grefmult.push_back(Norm(zdc_grefmult[prefix]->ProjectionY()));
  }

  // print comparisons and ratios between the two datasets
  dijetcore::PrintWithRatio(vx, legend_labels, hopts, coptslogy, FLAGS_outputDir,
              "vx", "V_{x} [cm]", "fraction");
  dijetcore::PrintWithRatio(vy, legend_labels, hopts, coptslogy, FLAGS_outputDir,
              "vy", "V_{y} [cm]", "fraction");
  dijetcore::PrintWithRatio(vz, legend_labels, hopts, copts, FLAGS_outputDir,
              "vz", "V_{z} [cm]", "fraction");
  dijetcore::PrintWithRatio(zdc, legend_labels, hopts, copts, FLAGS_outputDir,
              "zdc", "zdcX [kHz]", "fraction");
  dijetcore::PrintWithRatio(nvertices, legend_labels, hopts, copts, FLAGS_outputDir,
              "nvertices", "N_{vertices}", "fraction");
  dijetcore::PrintWithRatio(refmult, legend_labels, hopts, coptslogy, FLAGS_outputDir,
              "refmult", "refmult", "fraction");
  dijetcore::PrintWithRatio(grefmult, legend_labels, hopts, coptslogy, FLAGS_outputDir,
              "grefmult", "grefmult", "fraction");
  // PrintRatios(dvz, legend_labels, hopts, copts, FLAGS_outputDir,
  //             "dvz", "#DeltaV_{z}^{TPC, VPD} [cm]", "fraction");


  

  google::ShutDownCommandLineFlags();
  return 0;
}