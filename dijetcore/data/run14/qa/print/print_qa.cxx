#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/root/root_print_routines.h"

#include <set>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TStyle.h"

DIJETCORE_DEFINE_string(input, "", "root file with results from run_qa");
DIJETCORE_DEFINE_string(outputDir, "results", "directory for output");
DIJETCORE_DEFINE_string(histPrefixes, "",
                        "histogram prefixes, separated by commas");

using dijetcore::dijetcore_map;
using std::string;

int main(int argc, char *argv[]) {
  string usage = "Run 7 differential di-jet imbalance print routine";

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

  TFile input_file(FLAGS_input.c_str(), "READ");

  if (!input_file.IsOpen()) {
    LOG(ERROR) << "Could not open input file: " << FLAGS_input << ", exiting";
    return 1;
  }

  // first, parse all histogram prefixes
  std::set<string> prefixes;
  dijetcore::SplitString(FLAGS_histPrefixes, prefixes, ",");
  prefixes.insert("");

  // we'll store each histogram in a map - key will be the prefix
  dijetcore_map<string, TH2F *> run_id_refmult;
  dijetcore_map<string, TH2F *> run_id_grefmult;
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
  dijetcore_map<string, TH3F *> run_id_tower_e;
  dijetcore_map<string, TH3F *> run_id_tower_et;
  dijetcore_map<string, TH3F *> run_id_tower_adc;

  // event QA
  dijetcore_map<string, TH2F *> zdc_vz;
  dijetcore_map<string, TH2F *> vz_vx;
  dijetcore_map<string, TH2F *> vz_vy;
  dijetcore_map<string, TH2F *> zdc_refmult;
  dijetcore_map<string, TH2F *> bbc_refmult;
  dijetcore_map<string, TH1F *> n_vertices;

  // track QA
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

  // tower QA
  dijetcore_map<string, TH2F *> e_et;
  dijetcore_map<string, TH2F *> zdc_e;
  dijetcore_map<string, TH2F *> zdc_et;
  dijetcore_map<string, TH2F *> zdc_adc;
  dijetcore_map<string, TH2F *> tow_eta_phi;

  // load histograms for each prefix
  for (auto &prefix : prefixes) {
    run_id_refmult[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidref").c_str());
    run_id_grefmult[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidgref").c_str());
    run_id_nprim_nglob[prefix] = (TH3F*) input_file.Get(dijetcore::MakeString(prefix, "runidnprimnglob").c_str());
    run_id_zdc[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidzdc").c_str());
    run_id_bbc[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidbbc").c_str());
    run_id_vzvpdvz[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runiddvz").c_str());
    run_id_vz[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidvz").c_str());
    run_id_vx_vy[prefix] = (TH3F*) input_file.Get(dijetcore::MakeString(prefix, "runidvxvy").c_str());
    run_id_track_pt[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidpt").c_str());
    run_id_track_px[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidpx").c_str());
    run_id_track_py[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidpy").c_str());
    run_id_track_pz[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidpz").c_str());
    run_id_track_eta[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runideta").c_str());
    run_id_track_phi[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidphi").c_str());
    run_id_track_dca[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runiddca").c_str());
    run_id_track_nhit[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidnhit").c_str());
    run_id_track_nhitposs[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidnhitposs").c_str());
    run_id_track_nhitfrac[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "runidnhitfrac").c_str());
    run_id_tower_e[prefix] = (TH3F*) input_file.Get(dijetcore::MakeString(prefix, "runide").c_str());
    run_id_tower_et[prefix] = (TH3F*) input_file.Get(dijetcore::MakeString(prefix, "runidet").c_str());
    run_id_tower_adc[prefix] = (TH3F*) input_file.Get(dijetcore::MakeString(prefix, "runidadc").c_str());
    zdc_vz[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcvz").c_str());
    vz_vx[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "vzvx").c_str());
    vz_vy[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "vzvy").c_str());
    zdc_refmult[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcref").c_str());
    bbc_refmult[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "bbcref").c_str());
    n_vertices[prefix] = (TH1F*) input_file.Get(dijetcore::MakeString(prefix, "nvertex").c_str());
    px_py[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "pxpy").c_str());
    pz_px[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "pzpx").c_str());
    pz_py[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "pzpy").c_str());
    zdc_px[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcpx").c_str());
    zdc_py[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcpy").c_str());
    zdc_pz[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcpz").c_str());
    zdc_pt[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcpt").c_str());
    zdc_dca[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcdca").c_str());
    zdc_nhit[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcnhit").c_str());
    zdc_nhitposs[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcnhitposs").c_str());
    zdc_nhitfrac[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcnhitfrac").c_str());
    zdc_eta[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdceta").c_str());
    zdc_phi[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcphi").c_str());
    eta_phi[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "tracketaphi").c_str());
    e_et[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "eet").c_str());
    zdc_e[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdce").c_str());
    zdc_et[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcet").c_str());
    zdc_adc[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "zdcadc").c_str());
    tow_eta_phi[prefix] = (TH2F*) input_file.Get(dijetcore::MakeString(prefix, "towetaphi").c_str());
  }

  gflags::ShutDownCommandLineFlags();
  return 0;
}
