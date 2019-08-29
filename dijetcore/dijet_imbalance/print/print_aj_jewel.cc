#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/root/root_print_utils.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "sct/lib/root_ext/TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector2.h"

#include "dijetcore/lib/map.h"
#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "",
                        "input root file with processed jewel di-jet trees");
DIJETCORE_DEFINE_string(outputDir, "results", "directory for output");
DIJETCORE_DEFINE_bool(
    setScanning, true,
    "fixes the initial radius, and scans through matched radii");
DIJETCORE_DEFINE_double(initRadius, 0.2, "initial radius when set to scan");
DIJETCORE_DEFINE_string(radii, "0.2,0.25,0.3,0.35,0.4", "radii to put in grid");
DIJETCORE_DEFINE_string(constPt, "1.0,1.5,2.0,2.5,3.0", "radii to put in grid");
DIJETCORE_DEFINE_string(centrality, "0-20%",
                        "string for centrality definitions (cosmetic)");
DIJETCORE_DEFINE_double(trigger, 5.4, "trigger ET threshold");

// adds to the map all TTrees that conform to
// the DijetWorker's naming convention
void GetTreesFromFile(const TFile &file,
                      std::unordered_map<string, TTree *> &map) {
  TKey *key;
  TIter next(file.GetListOfKeys());

  while ((key = (TKey *)next())) {
    // check if its a TTree
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TTree"))
      continue;

    // check if its produced by the DijetWorker
    // (in a rather handwavy fashion)
    string tmp(key->GetName());
    if (tmp.find("LEAD_INIT") == string::npos ||
        tmp.find("SUB_INIT") == string::npos)
      continue;

    map.insert({tmp, (TTree *)key->ReadObj()});
  }
}

int main(int argc, char *argv[]) {
  string usage = "JEWEL differential di-jet imbalance print routine";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // set drawing preferences for histograms and graphs
  // gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  // turn off print messages
  gErrorIgnoreLevel = kInfo + 1;

  // check to make sure we have valid inputs
  if (!boost::filesystem::exists(FLAGS_input)) {
    LOG(ERROR) << "input file " << FLAGS_input << "doesn't exist: exiting"
               << std::endl;
    return 1;
  }

  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);

  // read in the file
  TFile input_file(FLAGS_input.c_str(), "READ");

  // and the radii & constituent pT we want to use
  std::vector<double> radii =
      dijetcore::ParseArgStringToVec<double>(FLAGS_radii);
  std::sort(radii.begin(), radii.end());
  std::vector<string> radii_string;
  for (auto &val : radii) {
    std::stringstream stream;
    stream << val;
    radii_string.push_back(stream.str());
  }
  std::vector<double> constpt =
      dijetcore::ParseArgStringToVec<double>(FLAGS_constPt);
  std::sort(constpt.begin(), constpt.end());
  std::vector<string> constpt_string;
  for (auto &val : constpt) {
    std::stringstream stream;
    stream << val;
    constpt_string.push_back(stream.str());
  }
  if (FLAGS_setScanning) {
    LOG(INFO) << "match R x const Pt matrix selected";
    LOG(INFO) << "match R: " << radii_string;
    LOG(INFO) << "constPt: " << constpt_string;
  } else {
    LOG(INFO) << "R x const Pt matrix selected";
    LOG(INFO) << "Radius: " << radii_string;
    LOG(INFO) << "constPt: " << constpt_string;
  }

  // get largest radii - this will set our eta cut for all jets
  double max_radii = *radii.rbegin();
  LOG(INFO) << "maximum radius = " << max_radii;

  // set our max jet eta
  double max_eta = 1.0 - max_radii;

  LOG(INFO) << "so our maximum jet eta = " << max_eta;

  if (max_eta < 0 || max_eta > 1.0) {
    LOG(ERROR) << "something weird is going on with max jet radius...";
    return 1;
  }

  // now we'll get the trees from the files, ignoring any objects
  // in the file that don't conform to the naming conventions from
  // the DijetWorker. There are also coincidence histograms to save
  std::vector<string> keys;
  std::unordered_map<std::string, dijetcore::DijetKey> parsed_keys;
  std::unordered_map<string, TTree *> trees;

  GetTreesFromFile(input_file, trees);

  // parse keys
  for (auto entry : trees) {
    keys.push_back(entry.first);
    parsed_keys[entry.first] = dijetcore::ParseStringToDijetKey(entry.first);
  }

  // now sort the keys into their grid spots
  std::vector<std::vector<string>> grid_keys;
  std::vector<std::vector<dijetcore::DijetKey>> grid_key_params;
  for (int rad = 0; rad < radii.size(); ++rad) {
    grid_keys.push_back(std::vector<string>());
    grid_key_params.push_back(std::vector<dijetcore::DijetKey>());
    for (int pt = 0; pt < constpt.size(); ++pt) {
      // find if the requested tree is present
      for (auto &key_string : keys) {
        dijetcore::DijetKey &key = parsed_keys[key_string];

        if (FLAGS_setScanning && key.lead_init_r == FLAGS_initRadius &&
            key.sub_init_r == FLAGS_initRadius &&
            key.lead_match_r == radii[rad] && key.sub_match_r == radii[rad] &&
            key.lead_init_const_pt == constpt[pt] &&
            key.sub_init_const_pt == constpt[pt]) {
          grid_keys[rad].push_back(key_string);
          grid_key_params[rad].push_back(key);
        } else if (!FLAGS_setScanning && key.lead_init_r == radii[rad] &&
                   key.sub_init_r == radii[rad] &&
                   key.lead_match_r == radii[rad] &&
                   key.sub_match_r == radii[rad] &&
                   key.lead_init_const_pt == constpt[pt] &&
                   key.sub_init_const_pt == constpt[pt]) {
          grid_keys[rad].push_back(key_string);
          grid_key_params[rad].push_back(key);
        }
      }

      if (grid_keys[rad].size() < pt + 1) {
        LOG(INFO) << "could not find a key for R=" << radii[rad]
                  << " const Pt=" << constpt[pt];
        return 1;
      }
      if (grid_keys[rad].size() > pt + 1) {
        LOG(INFO) << "found multiple keys for R=" << radii[rad]
                  << " const Pt=" << constpt[pt];
        return 1;
      }
    }
  }

  // save all histograms so we can do comparisons
  // between different keys if we want

  std::unordered_map<string, TH1D *> event_vx;
  std::unordered_map<string, TH1D *> event_vy;
  std::unordered_map<string, TH2D *> event_vxvy;

  std::unordered_map<string, TH1D *> event_vx_rotated;
  std::unordered_map<string, TH1D *> event_vy_rotated;
  std::unordered_map<string, TH2D *> event_vxvy_rotated;

  std::unordered_map<string, TH1D *> xsec_dist;
  std::unordered_map<string, TH1D *> parton_pt;

  std::unordered_map<string, TH1D *> hard_lead_eta;
  std::unordered_map<string, TH1D *> hard_lead_phi;
  std::unordered_map<string, TH1D *> hard_lead_pt;
  std::unordered_map<string, TH1D *> hard_lead_const;
  std::unordered_map<string, TH1D *> hard_lead_rho;
  std::unordered_map<string, TH1D *> hard_lead_sig;
  std::unordered_map<string, TH1D *> match_lead_eta;
  std::unordered_map<string, TH1D *> match_lead_phi;
  std::unordered_map<string, TH1D *> match_lead_pt;
  std::unordered_map<string, TH1D *> match_lead_const;
  std::unordered_map<string, TH1D *> match_lead_rho;
  std::unordered_map<string, TH1D *> match_lead_sig;

  std::unordered_map<string, TH1D *> hard_sub_eta;
  std::unordered_map<string, TH1D *> hard_sub_phi;
  std::unordered_map<string, TH1D *> hard_sub_pt;
  std::unordered_map<string, TH1D *> hard_sub_const;
  std::unordered_map<string, TH1D *> hard_sub_rho;
  std::unordered_map<string, TH1D *> hard_sub_sig;
  std::unordered_map<string, TH1D *> match_sub_eta;
  std::unordered_map<string, TH1D *> match_sub_phi;
  std::unordered_map<string, TH1D *> match_sub_pt;
  std::unordered_map<string, TH1D *> match_sub_const;
  std::unordered_map<string, TH1D *> match_sub_rho;
  std::unordered_map<string, TH1D *> match_sub_sig;

  std::unordered_map<string, TH1D *> npart;
  std::unordered_map<string, TH1D *> hard_dphi;
  std::unordered_map<string, TH1D *> match_dphi;
  std::unordered_map<string, TH1D *> lead_dr;
  std::unordered_map<string, TH1D *> sub_dr;
  std::unordered_map<string, TH1D *> lead_dpt;
  std::unordered_map<string, TH1D *> sub_dpt;
  std::unordered_map<string, TH1D *> lead_dpt_frac;
  std::unordered_map<string, TH1D *> sub_dpt_frac;

  std::unordered_map<string, TH1D *> hard_aj;
  std::unordered_map<string, TH1D *> match_aj;

  std::unordered_map<string, TH1D *> hard_pp_only_aj;
  std::unordered_map<string, TH1D *> match_pp_only_aj;
  std::unordered_map<string, TH1D *> hard_lead_pp_dpt;
  std::unordered_map<string, TH1D *> match_lead_pp_dpt;
  std::unordered_map<string, TH1D *> hard_sub_pp_dpt;
  std::unordered_map<string, TH1D *> match_sub_pp_dpt;
  std::unordered_map<string, TH1D *> hard_lead_pp_dpt_frac;
  std::unordered_map<string, TH1D *> match_lead_pp_dpt_frac;
  std::unordered_map<string, TH1D *> hard_sub_pp_dpt_frac;
  std::unordered_map<string, TH1D *> match_sub_pp_dpt_frac;
  std::unordered_map<string, TH1D *> pp_lead_hard_match_fraction;
  std::unordered_map<string, TH1D *> pp_sub_hard_match_fraction;

  // create our histogram and canvas options
  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;
  copts.leg_left_bound = 0.6;
  dijetcore::canvasOpts coptslogz;
  coptslogz.leg_left_bound = 0.6;
  coptslogz.log_z = true;
  dijetcore::canvasOpts coptslogy;
  coptslogy.do_legend = false;
  coptslogy.leg_left_bound = 0.6;
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

  // count which entry we're on
  int entry = -1;
  for (auto &key : keys) {
    // increment entry counter
    entry++;
    LOG(INFO) << "loading: " << key;
    // make our output directory
    string out_loc = FLAGS_outputDir + "/" + key + "/jewel";
    boost::filesystem::path dir(out_loc.c_str());
    boost::filesystem::create_directories(dir);

    // load our reader
    TTree *current_tree = trees[key];
    TTreeReader reader(current_tree);

    // create readervalues for auau first
    dijetcore::unique_ptr<TTreeReaderValue<double>> vx =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "vx");
    dijetcore::unique_ptr<TTreeReaderValue<double>> vy =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "vy");
    dijetcore::unique_ptr<TTreeReaderValue<double>> w =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "w");
    dijetcore::unique_ptr<TTreeReaderValue<double>> trigger_e =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "e");
    dijetcore::unique_ptr<TTreeReaderValue<double>> xsec =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "xsec");

    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pl =
        dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "pl");
    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> ps =
        dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "ps");

    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jl =
        dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "jl");
    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> js =
        dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "js");
    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jlm =
        dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "jlm");
    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jsm =
        dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "jsm");
    dijetcore::unique_ptr<TTreeReaderValue<int>> jlconst =
        dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jlconst");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jlrho =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jlrho");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jlsig =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jlsig");
    dijetcore::unique_ptr<TTreeReaderValue<int>> jlmconst =
        dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jlmconst");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jlmrho =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jlmrho");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jlmsig =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jlmsig");
    dijetcore::unique_ptr<TTreeReaderValue<int>> jsconst =
        dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jsconst");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jsrho =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jsrho");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jssig =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jssig");
    dijetcore::unique_ptr<TTreeReaderValue<int>> jsmconst =
        dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jsmconst");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jsmrho =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jsmrho");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jsmsig =
        dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jsmsig");

    // build prefix for names so histograms don't get confused
    // (root fuckin sucks)
    std::string hist_prefix = "key_" + std::to_string(entry);

    event_vx[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "eventvx").c_str(), "", 100,
                 -10, 10);
    event_vy[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "eventvy").c_str(), "", 100,
                 -10, 10);
    event_vxvy[key] =
        new TH2D(dijetcore::MakeString(hist_prefix, "eventvxvy").c_str(), "",
                 100, -10, 10, 100, -10, 10);

    event_vx_rotated[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "eventvxrot").c_str(), "",
                 100, -10, 10);
    event_vy_rotated[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "eventvyrot").c_str(), "",
                 100, -10, 10);
    event_vxvy_rotated[key] =
        new TH2D(dijetcore::MakeString(hist_prefix, "eventvxvyrot").c_str(), "",
                 100, -10, 10, 100, -10, 10);
    
    xsec_dist[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "xsec").c_str(), "",
                 100, 0, 4e6);
    parton_pt[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "parton_pt").c_str(), "",
                 100, 0, 100);

    // create the histograms
    hard_lead_eta[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardleadeta").c_str(), "",
                 50, -1, 1);
    hard_lead_phi[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardleadphi").c_str(), "",
                 50, -TMath::Pi(), TMath::Pi());
    hard_lead_pt[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardleadpt").c_str(), "",
                 50, 0, 50);
    hard_lead_const[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardleadconst").c_str(),
                 "", 50, 0, 100);
    hard_lead_rho[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardleadrho").c_str(), "",
                 50, 0, 100);
    hard_lead_sig[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardleadsig").c_str(), "",
                 50, 0, 20);
    match_lead_eta[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchleadeta").c_str(), "",
                 50, -1, 1);
    match_lead_phi[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchleadphi").c_str(), "",
                 50, -TMath::Pi(), TMath::Pi());
    match_lead_pt[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchleadpt").c_str(), "",
                 50, 0, 50);
    match_lead_const[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchleadconst").c_str(),
                 "", 50, 0, 100);
    match_lead_rho[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchleadrho").c_str(), "",
                 50, 0, 100);
    match_lead_sig[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchleadsig").c_str(), "",
                 50, 0, 20);
    hard_sub_eta[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardsubeta").c_str(), "",
                 50, -1, 1);
    hard_sub_phi[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardsubphi").c_str(), "",
                 50, -TMath::Pi(), TMath::Pi());
    hard_sub_pt[key] = new TH1D(
        dijetcore::MakeString(hist_prefix, "hardsubpt").c_str(), "", 50, 0, 50);
    hard_sub_const[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardsubconst").c_str(), "",
                 50, 0, 100);
    hard_sub_rho[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardsubrho").c_str(), "",
                 50, 0, 100);
    hard_sub_sig[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardsubsig").c_str(), "",
                 50, 0, 20);
    match_sub_eta[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchsubeta").c_str(), "",
                 50, -1, 1);
    match_sub_phi[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchsubphi").c_str(), "",
                 50, -TMath::Pi(), TMath::Pi());
    match_sub_pt[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchsubpt").c_str(), "",
                 50, 0, 50);
    match_sub_const[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchsubconst").c_str(),
                 "", 50, 0, 100);
    match_sub_rho[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchsubrho").c_str(), "",
                 50, 0, 100);
    match_sub_sig[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchsubsig").c_str(), "",
                 50, 0, 20);
    npart[key] = new TH1D(dijetcore::MakeString(hist_prefix, "npart").c_str(),
                          "", 100, 0, 1200);
    hard_dphi[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "harddphi").c_str(), "", 50,
                 0, TMath::Pi());
    match_dphi[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchdphi").c_str(), "",
                 50, 0, TMath::Pi());
    lead_dr[key] = new TH1D(
        dijetcore::MakeString(hist_prefix, "leaddr").c_str(), "", 50, 0, 0.5);
    sub_dr[key] = new TH1D(dijetcore::MakeString(hist_prefix, "subdr").c_str(),
                           "", 50, 0, 0.5);
    lead_dpt[key] = new TH1D(
        dijetcore::MakeString(hist_prefix, "leaddpt").c_str(), "", 50, -20, 20);
    sub_dpt[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "subdptfrac").c_str(), "",
                 50, -20, 20);
    lead_dpt_frac[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "leaddptfrac").c_str(), "",
                 50, -2.0, 2.0);
    sub_dpt_frac[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "subdpt").c_str(), "", 50,
                 -2.0, 2.0);
    hard_aj[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "hardaj").c_str(), "", 15,
                 0.0001, 0.9);
    match_aj[key] =
        new TH1D(dijetcore::MakeString(hist_prefix, "matchaj").c_str(), "", 15,
                 0.0001, 0.9);

    size_t total_events = reader.GetTree()->GetEntries();
    LOG(INFO) << "number of events: " << total_events;
    size_t event_counter = 0;
    while (reader.Next()) {
      if (event_counter % (total_events / 10) == 0)
        LOG(INFO) << "event: " << event_counter;
      event_counter++;

      // check jet eta
      if (FLAGS_setScanning &&
          (fabs((*jl)->Eta()) > max_eta || fabs((*js)->Eta()) > max_eta))
        continue;

      // check for trigger pT
      if (**trigger_e < FLAGS_trigger)
        continue;

      if (**xsec > 4000000)
        continue;

      event_vx[key]->Fill(**vx);
      event_vy[key]->Fill(**vy);
      event_vxvy[key]->Fill(**vx, **vy);

      // to fill the rotated entries, we first have to find the angle of
      // rotation
      TVector2 event_vector(**vx, **vy);
      TVector2 jet_vector((**jl).Px(), (**jl).Py());
      TVector2 event_vector_rotated =
          event_vector.Rotate(TMath::Pi() - jet_vector.Phi());

      event_vx_rotated[key]->Fill(event_vector_rotated.X());
      event_vy_rotated[key]->Fill(event_vector_rotated.Y());
      event_vxvy_rotated[key]->Fill(event_vector_rotated.X(),
                                    event_vector_rotated.Y(), // 1.0);
                                    **xsec);
      xsec_dist[key]->Fill(**xsec);
      parton_pt[key]->Fill((*pl)->Pt(), **xsec);

      hard_lead_eta[key]->Fill((*jl)->Eta());
      hard_lead_phi[key]->Fill((*jl)->Phi());
      hard_lead_pt[key]->Fill((*jl)->Pt());
      hard_lead_const[key]->Fill(**jlconst);
      hard_lead_rho[key]->Fill(**jlrho);
      hard_lead_sig[key]->Fill(**jlsig);
      match_lead_eta[key]->Fill((*jlm)->Eta());
      match_lead_phi[key]->Fill((*jlm)->Phi());
      match_lead_pt[key]->Fill((*jlm)->Pt());
      match_lead_const[key]->Fill(**jlmconst);
      match_lead_rho[key]->Fill(**jlmrho);
      match_lead_sig[key]->Fill(**jlmsig);

      hard_sub_eta[key]->Fill((*js)->Eta());
      hard_sub_phi[key]->Fill((*js)->Phi());
      hard_sub_pt[key]->Fill((*js)->Pt());
      hard_sub_const[key]->Fill(**jsconst);
      hard_sub_rho[key]->Fill(**jsrho);
      hard_sub_sig[key]->Fill(**jssig);
      match_sub_eta[key]->Fill((*jsm)->Eta());
      match_sub_phi[key]->Fill((*jsm)->Phi());
      match_sub_pt[key]->Fill((*jsm)->Pt());
      match_sub_const[key]->Fill(**jsmconst);
      match_sub_rho[key]->Fill(**jsmrho);
      match_sub_sig[key]->Fill(**jsmsig);

      hard_dphi[key]->Fill(fabs((*jl)->DeltaPhi(**js)));
      match_dphi[key]->Fill(fabs((*jlm)->DeltaPhi(**jsm)));
      lead_dr[key]->Fill((*jl)->DeltaR(**jlm));
      sub_dr[key]->Fill((*js)->DeltaR(**jsm));
      lead_dpt[key]->Fill((*jl)->Pt() - (*jlm)->Pt());
      sub_dpt[key]->Fill((*js)->Pt() - (*jsm)->Pt());
      lead_dpt_frac[key]->Fill(((*jl)->Pt() - (*jlm)->Pt()) / (*jl)->Pt());
      sub_dpt_frac[key]->Fill(((*js)->Pt() - (*jsm)->Pt()) / (*js)->Pt());

      hard_aj[key]->Fill(fabs((*jl)->Pt() - (*js)->Pt()) /
                         ((*jl)->Pt() + (*js)->Pt()));
      match_aj[key]->Fill(fabs((*jlm)->Pt() - (*jsm)->Pt()) /
                          ((*jlm)->Pt() + (*jsm)->Pt()));
    }

    // pave text for hard & matched jets
    dijetcore::DijetKey params = parsed_keys[key];
    std::stringstream streamHard;
    std::stringstream streamMatch;
    std::stringstream streamLeadPt;
    std::stringstream streamSubPt;
    std::stringstream streamConstPt;
    streamHard << params.lead_init_r;
    streamMatch << params.lead_match_r;
    streamLeadPt << params.lead_init_pt;
    streamSubPt << params.sub_init_pt;
    streamConstPt << params.lead_init_const_pt;

    TPaveText hardPave(0.68, 0.3, 0.88, 0.6, "NB NDC");
    hardPave.SetFillStyle(0);
    hardPave.SetBorderSize(0);
    hardPave
        .AddText(dijetcore::MakeString("JEWEL, ", FLAGS_centrality).c_str())
        ->SetTextSize(0.038);
    hardPave
        .AddText(
            dijetcore::MakeString("anti-k_{T}, R = ", streamHard.str()).c_str())
        ->SetTextSize(0.038);
    hardPave
        .AddText(dijetcore::MakeString("p_{T}^{hard const} > ",
                                       streamConstPt.str(), "GeV/c")
                     .c_str())
        ->SetTextSize(0.038);
    hardPave
        .AddText(dijetcore::MakeString("p_{T}^{lead} > ", streamLeadPt.str(),
                                       "GeV/c")
                     .c_str())
        ->SetTextSize(0.038);
    hardPave
        .AddText(dijetcore::MakeString("p_{T}^{sublead} > ", streamSubPt.str(),
                                       "GeV/c")
                     .c_str())
        ->SetTextSize(0.038);
    dijetcore::Print2DSimple(event_vxvy_rotated[key], hopts, copts, out_loc,
                             "rotated_vxvy", "", "x [fm]", "y [fm]", "CONTZ");
    dijetcore::PrettyPrint1D(xsec_dist[key], hopts, coptslogy, "", out_loc, "xsec", "", "xsec", "counts");
    dijetcore::PrettyPrint1D(parton_pt[key], hopts, coptslogy, "", out_loc, "parton_pt", "", "p_{T}", "counts");
  }

  // now we will print out the average shift away from zero
  std::vector<TGraphErrors> shift_x(radii.size(), TGraphErrors(constpt.size()));
  std::vector<TGraphErrors> shift_x2(radii.size(),
                                     TGraphErrors(constpt.size()));
  std::vector<TGraphErrors> shift_y(radii.size(), TGraphErrors(constpt.size()));
  std::vector<TGraphErrors> shift_y2(radii.size(),
                                     TGraphErrors(constpt.size()));
  std::vector<TGraphErrors> shift_r(radii.size(), TGraphErrors(constpt.size()));

  for (int rad = 0; rad < radii.size(); ++rad) {
    for (int pt = 0; pt < constpt.size(); ++pt) {
      auto &key = grid_keys[rad][pt];
      TH2D *h = event_vxvy_rotated[key];

      shift_x[rad].SetPoint(pt, constpt[pt], h->GetMean(1));
      shift_x2[rad].SetPoint(pt, constpt[pt], h->GetStdDev(1));
      shift_y[rad].SetPoint(pt, constpt[pt], h->GetMean(2));
      shift_y2[rad].SetPoint(pt, constpt[pt], h->GetStdDev(2));
      double mean_r = sqrt(pow(h->GetMean(1), 2) + pow(h->GetMean(2), 2));
      shift_r[rad].SetPoint(pt, constpt[pt], mean_r);
    }
  }

  TCanvas c1;
  shift_r[0].SetTitle("");
  shift_r[0].GetXaxis()->SetTitle("p_T^{const} [GeV/c]");
  shift_r[0].GetYaxis()->SetTitle("<r> [fm]");
  shift_r[0].GetYaxis()->SetRangeUser(0, 1);
  shift_r[0].Draw();
  shift_r[0].SetLineColor(kBlack);
  TLegend *legend1 = new TLegend(0.6, 0.7, 0.88, 0.88);
  legend1->AddEntry(&shift_r[0], radii_string[0].c_str(), "lep");
  for (int i = 1; i < shift_r.size(); ++i) {
    shift_r[i].SetLineColor(i + 1);
    shift_r[i].SetMarkerColor(i + 1);
    shift_r[i].Draw("SAME");
    legend1->AddEntry(&shift_r[i], radii_string[i].c_str(), "lep");
  }
  legend1->Draw();
  c1.SaveAs(dijetcore::MakeString(FLAGS_outputDir, "/r.pdf").c_str());

  TCanvas c2;
  shift_y2[0].SetTitle("");
  shift_y2[0].GetXaxis()->SetTitle("p_T^{const} [GeV/c]");
  shift_y2[0].GetYaxis()->SetTitle("<y^2> [fm^2]");
  shift_y2[0].GetYaxis()->SetRangeUser(2.7, 3.4);
  shift_y2[0].Draw();
  shift_y2[0].SetLineColor(kBlack);
  TLegend *legend2 = new TLegend(0.6, 0.7, 0.88, 0.88);
  legend2->AddEntry(&shift_y2[0], radii_string[0].c_str(), "lep");
  for (int i = 1; i < shift_y2.size(); ++i) {
    shift_y2[i].SetLineColor(i + 1);
    shift_y2[i].SetMarkerColor(i + 1);
    shift_y2[i].Draw("SAME");
    legend2->AddEntry(&shift_y2[i], radii_string[i].c_str(), "lep");
  }
  legend2->Draw();
  c2.SaveAs(dijetcore::MakeString(FLAGS_outputDir, "/y2.pdf").c_str());

  TCanvas c3;
  shift_x[0].SetTitle("");
  shift_x[0].GetXaxis()->SetTitle("p_T^{const} [GeV/c]");
  shift_x[0].GetYaxis()->SetTitle("<x> [fm]");
  shift_x[0].GetYaxis()->SetRangeUser(-1, 0);
  shift_x[0].Draw();
  shift_x[0].SetLineColor(kBlack);
  TLegend *legend3 = new TLegend(0.6, 0.7, 0.88, 0.88);
  legend2->AddEntry(&shift_x[0], radii_string[0].c_str(), "lep");
  for (int i = 0; i < shift_x.size(); ++i) {
    shift_x[i].SetLineColor(i + 1);
    shift_x[i].SetMarkerColor(i + 1);
    shift_x[i].Draw("SAME");
    legend3->AddEntry(&shift_x[i], radii_string[i].c_str(), "lep");
  }
  legend3->Draw();
  c3.SaveAs(dijetcore::MakeString(FLAGS_outputDir, "/x.pdf").c_str());

  return 0;
}

// loop over all keys
//   for (auto &key : keys) {
//     // increment entry counter
//     entry++;
//     LOG(INFO) << "loading: " << key;
//     // make our output directory
//     string key_loc = FLAGS_outputDir + "/" + key;
//     boost::filesystem::path dir(key_loc.c_str());
//     boost::filesystem::create_directories(dir);

//     // now iterate over all our datasets
//     for (auto type : TYPEMAP) {
//       auto &type_enum = type.first;
//       auto &data_index = type.second;

//       // make our output directory
//       string out_loc = key_loc + "/" + TYPENAME[type_enum];
//       boost::filesystem::path dir(out_loc.c_str());
//       boost::filesystem::create_directories(dir);

//       // load our reader
//       TTree *current_tree = trees[data_index][key];
//       TTreeReader reader(current_tree);

//       // create readervalues for auau first
//       dijetcore::unique_ptr<TTreeReaderValue<int>> runid =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader, "runid");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> eventid =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader, "eventid");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> vz =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader, "vz");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> refmult =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader, "refmult");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> grefmult =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader,
//           "grefmult");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> refmultcorr =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//                                                            "refmultcorr");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> grefmultcorr =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//                                                            "grefmultcorr");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> cent =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader, "cent");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> zdcrate =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "zdcrate");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> rp =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader, "rp");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> nglobal =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader, "nglobal");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> nprt =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader, "npart");
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jl =
//           dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
//                                                                    "jl");
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> js =
//           dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
//                                                                    "js");
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jlm =
//           dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
//                                                                    "jlm");
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jsm =
//           dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
//                                                                    "jsm");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> jlconst =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jlconst");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jlrho =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "jlrho");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jlsig =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "jlsig");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> jlmconst =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader,
//           "jlmconst");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jlmrho =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "jlmrho");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jlmsig =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "jlmsig");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> jsconst =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jsconst");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jsrho =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "jsrho");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jssig =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "jssig");
//       dijetcore::unique_ptr<TTreeReaderValue<int>> jsmconst =
//           dijetcore::make_unique<TTreeReaderValue<int>>(reader,
//           "jsmconst");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jsmrho =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "jsmrho");
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jsmsig =
//           dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//           "jsmsig");

//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jloa;
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jsoa;
//       dijetcore::unique_ptr<TTreeReaderValue<int>> jloaconst;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jloarho;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jloasig;
//       dijetcore::unique_ptr<TTreeReaderValue<int>> jsoaconst;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jsoarho;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> jsoasig;

//       // build all embedding branches (only used if pp_embedded is true)
//       dijetcore::unique_ptr<TTreeReaderValue<int>> embed_eventid;
//       dijetcore::unique_ptr<TTreeReaderValue<int>> embed_runid;
//       dijetcore::unique_ptr<TTreeReaderValue<int>> embed_refmult;
//       dijetcore::unique_ptr<TTreeReaderValue<int>> embed_grefmult;
//       dijetcore::unique_ptr<TTreeReaderValue<int>> embed_nprt;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> embed_refmultcorr;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> embed_grefmultcorr;
//       dijetcore::unique_ptr<TTreeReaderValue<int>> embed_cent;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> embed_rp;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> embed_zdcrate;
//       dijetcore::unique_ptr<TTreeReaderValue<double>> embed_vz;

//       dijetcore::unique_ptr<TTreeReaderValue<bool>> only_found_match;
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pp_only_jl;
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pp_only_js;
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pp_only_jlm;
//       dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pp_only_jsm;

//       // check if we have off-axis entries (for auau)
//       bool auau_off_axis_present = false;
//       if (data_index == TYPEMAP[DATATYPE::AUAU] &&
//           current_tree->GetBranch("jloa"))
//         auau_off_axis_present = true;

//       // check if we have pp-only cluster results
//       bool pp_only_present = false;
//       if (data_index == TYPEMAP[DATATYPE::PP] &&
//           current_tree->GetBranch("foundpp"))
//         pp_only_present = true;

//       // and check if we have embedding (for pp tree)
//       bool pp_embedded_present = false;
//       if (data_index == TYPEMAP[DATATYPE::PP] &&
//           current_tree->GetBranch("embed_eventid"))
//         pp_embedded_present = true;

//       // initialize source dependent reader values

//       if (auau_off_axis_present) {
//         jloa =
//         dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
//                                                                         "jloa");
//         jsoa =
//         dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
//                                                                         "jsoa");
//         jloaconst =
//             dijetcore::make_unique<TTreeReaderValue<int>>(reader,
//             "jloaconst");
//         jloarho =
//             dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//             "jloarho");
//         jloasig =
//             dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//             "jloasig");
//         jsoaconst =
//             dijetcore::make_unique<TTreeReaderValue<int>>(reader,
//             "jsoaconst");
//         jsoarho =
//             dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//             "jsoarho");
//         jsoasig =
//             dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//             "jsoasig");
//       }

//       if (pp_embedded_present) {
//         embed_eventid = dijetcore::make_unique<TTreeReaderValue<int>>(
//             reader, "embed_eventid");
//         embed_runid = dijetcore::make_unique<TTreeReaderValue<int>>(
//             reader, "embed_runid");
//         embed_refmult = dijetcore::make_unique<TTreeReaderValue<int>>(
//             reader, "embed_refmult");
//         embed_grefmult = dijetcore::make_unique<TTreeReaderValue<int>>(
//             reader, "embed_grefmult");
//         embed_nprt = dijetcore::make_unique<TTreeReaderValue<int>>(
//             reader, "embed_npart");
//         embed_refmultcorr =
//         dijetcore::make_unique<TTreeReaderValue<double>>(
//             reader, "embed_refmultcorr");
//         embed_grefmultcorr =
//         dijetcore::make_unique<TTreeReaderValue<double>>(
//             reader, "embed_grefmultcorr");
//         embed_cent =
//             dijetcore::make_unique<TTreeReaderValue<int>>(reader,
//             "embed_cent");
//         embed_rp = dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//                                                                     "embed_rp");
//         embed_zdcrate = dijetcore::make_unique<TTreeReaderValue<double>>(
//             reader, "embed_zdcrate");
//         embed_vz = dijetcore::make_unique<TTreeReaderValue<double>>(reader,
//                                                                     "embed_vz");
//       }

//       if (pp_only_present) {
//         only_found_match =
//             dijetcore::make_unique<TTreeReaderValue<bool>>(reader,
//             "foundpp");
//         pp_only_jl =
//         dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(
//             reader, "ppjl");
//         pp_only_js =
//         dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(
//             reader, "ppjs");
//         pp_only_jlm =
//         dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(
//             reader, "ppjlm");
//         pp_only_jsm =
//         dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(
//             reader, "ppjsm");
//       }

//       // build prefix for names so histograms don't get confused
//       // (root fuckin sucks)
//       std::string key_prefix = "key_" + std::to_string(entry) + "_";
//       std::string datatype_prefix = TYPENAME[type_enum];
//       std::string hist_prefix = key_prefix + datatype_prefix;

//       // create the histograms
//       hard_lead_eta[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardleadeta").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, -1, 1);
//       hard_lead_phi[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardleadphi").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, -TMath::Pi(),
//                    TMath::Pi());
//       hard_lead_pt[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardleadpt").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, 50);
//       hard_lead_const[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardleadconst").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 100);
//       hard_lead_rho[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardleadrho").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 100);
//       hard_lead_sig[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardleadsig").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 20);

//       match_lead_eta[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchleadeta").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, -1, 1);
//       match_lead_phi[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchleadphi").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, -TMath::Pi(),
//                    TMath::Pi());
//       match_lead_pt[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchleadpt").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 50);
//       match_lead_const[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchleadconst").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 100);
//       match_lead_rho[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchleadrho").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 100);
//       match_lead_sig[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchleadsig").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 20);

//       hard_sub_eta[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardsubeta").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, -1, 1);
//       hard_sub_phi[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardsubphi").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, -TMath::Pi(), TMath::Pi());
//       hard_sub_pt[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "hardsubpt").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, 50);
//       hard_sub_const[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardsubconst").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 100);
//       hard_sub_rho[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardsubrho").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, 100);
//       hard_sub_sig[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardsubsig").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, 20);

//       match_sub_eta[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchsubeta").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, -1, 1);
//       match_sub_phi[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchsubphi").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, -TMath::Pi(),
//                    TMath::Pi());
//       match_sub_pt[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchsubpt").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, 50);
//       match_sub_const[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchsubconst").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 100);
//       match_sub_rho[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchsubrho").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 100);
//       match_sub_sig[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchsubsig").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, 0, 20);

//       npart[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "npart").c_str(), "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 100, 0, 1200);
//       hard_dphi[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "harddphi").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, TMath::Pi());
//       match_dphi[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "matchdphi").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, TMath::Pi());
//       lead_dr[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "leaddr").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, 0.5);
//       sub_dr[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "subdr").c_str(), "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, 0, 0.5);
//       lead_dpt[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "leaddpt").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, -20, 20);
//       sub_dpt[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "subdptfrac").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, -20, 20);
//       lead_dpt_frac[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "leaddptfrac").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 50, -2.0, 2.0);
//       sub_dpt_frac[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "subdpt").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 50, -2.0, 2.0);

//       hard_aj[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "hardaj").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 15, 0.0001, 0.9);
//       match_aj[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix, "matchaj").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 15, 0.0001, 0.9);
//       hard_aj_test[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "hardajtest").c_str(),
//           "",
//                    cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                    0.5, 10000, 0.0001, 0.9);
//       match_aj_test[data_index][key] =
//           new TH2D(dijetcore::MakeString(hist_prefix,
//           "matchajtest").c_str(),
//                    "", cent_boundaries.size(), -0.5,
//                    cent_boundaries.size() - 0.5, 10000, 0.0001, 0.9);

//       if (auau_off_axis_present) {
//         off_axis_aj[data_index][key] =
//             new TH2D(dijetcore::MakeString(hist_prefix,
//             "offaxisaj").c_str(),
//                      "", cent_boundaries.size(), -0.5,
//                      cent_boundaries.size() - 0.5, 15, 0.0001, 0.9);
//         off_axis_aj_test[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix, "offaxisajtest").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             10000, 0.0001, 0.9);
//         off_axis_rho[data_index][key] =
//             new TH2D(dijetcore::MakeString(hist_prefix,
//             "offaxisrho").c_str(),
//                      "", cent_boundaries.size(), -0.5,
//                      cent_boundaries.size() - 0.5, 50, 0, 100);
//         off_axis_sig[data_index][key] =
//             new TH2D(dijetcore::MakeString(hist_prefix,
//             "offaxissig").c_str(),
//                      "", cent_boundaries.size(), -0.5,
//                      cent_boundaries.size() - 0.5, 50, 0, 20);
//       }
//       if (pp_only_present) {
//         hard_pp_only_aj[data_index][key] =
//             new TH2D(dijetcore::MakeString(hist_prefix,
//             "pphardaj").c_str(),
//             "",
//                      cent_boundaries.size(), -0.5, cent_boundaries.size() -
//                      0.5, 15, 0.0001, 0.9);
//         match_pp_only_aj[data_index][key] =
//             new TH2D(dijetcore::MakeString(hist_prefix,
//             "ppmatchaj").c_str(),
//                      "", cent_boundaries.size(), -0.5,
//                      cent_boundaries.size() - 0.5, 15, 0.0001, 0.9);
//         hard_lead_pp_dpt[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix, "pphardleaddpt").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             50, 0, 30);
//         match_lead_pp_dpt[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix, "ppmatchleaddpt").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             50, 0, 30);
//         hard_sub_pp_dpt[data_index][key] =
//             new TH2D(dijetcore::MakeString(hist_prefix,
//             "pphardsubdpt").c_str(),
//                      "", cent_boundaries.size(), -0.5,
//                      cent_boundaries.size() - 0.5, 50, 0, 30);
//         match_sub_pp_dpt[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix, "ppmatchsubdpt").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             50, 0, 30);
//         hard_lead_pp_dpt_frac[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix,
//             "pphardleaddptfrac").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             50, 0, 30);
//         match_lead_pp_dpt_frac[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix,
//             "ppmatchleaddptfrac").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             50, 0, 30);
//         hard_sub_pp_dpt_frac[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix, "pphardsubdptfrac").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             50, 0, 30);
//         match_sub_pp_dpt_frac[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix,
//             "ppmatchsubdptfrac").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             50, 0, 30);
//         pp_lead_hard_match_fraction[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix,
//             "ppleadhardmatchfrac").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             2, -0.5, 1.5);
//         pp_sub_hard_match_fraction[data_index][key] = new TH2D(
//             dijetcore::MakeString(hist_prefix,
//             "ppsubhardmatchfrac").c_str(),
//             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
//             2, -0.5, 1.5);
//       }

//       // event loop

//       while (reader.Next()) {
//         int cent_bin = -1;
//         for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
//           auto &pair = cent_bin_boundaries[i];
//           if (**cent >= pair.first && **cent <= pair.second) {
//             cent_bin = i;
//             break;
//           }
//         }

//         if (cent_bin == -1)
//           continue;

//         // check jet eta
//         if (FLAGS_setScanning &&
//             (fabs((*jl)->Eta()) > max_eta || fabs((*js)->Eta()) > max_eta))
//           continue;

//         hard_lead_eta[data_index][key]->Fill(cent_bin, (*jl)->Eta());
//         hard_lead_phi[data_index][key]->Fill(cent_bin, (*jl)->Phi());
//         hard_lead_pt[data_index][key]->Fill(cent_bin, (*jl)->Pt());
//         hard_lead_const[data_index][key]->Fill(cent_bin, **jlconst);
//         hard_lead_rho[data_index][key]->Fill(cent_bin, **jlrho);
//         hard_lead_sig[data_index][key]->Fill(cent_bin, **jlsig);
//         match_lead_eta[data_index][key]->Fill(cent_bin, (*jlm)->Eta());
//         match_lead_phi[data_index][key]->Fill(cent_bin, (*jlm)->Phi());
//         match_lead_pt[data_index][key]->Fill(cent_bin, (*jlm)->Pt());
//         match_lead_const[data_index][key]->Fill(cent_bin, **jlmconst);
//         match_lead_rho[data_index][key]->Fill(cent_bin, **jlmrho);
//         match_lead_sig[data_index][key]->Fill(cent_bin, **jlmsig);

//         hard_sub_eta[data_index][key]->Fill(cent_bin, (*js)->Eta());
//         hard_sub_phi[data_index][key]->Fill(cent_bin, (*js)->Phi());
//         hard_sub_pt[data_index][key]->Fill(cent_bin, (*js)->Pt());
//         hard_sub_const[data_index][key]->Fill(cent_bin, **jsconst);
//         hard_sub_rho[data_index][key]->Fill(cent_bin, **jsrho);
//         hard_sub_sig[data_index][key]->Fill(cent_bin, **jssig);
//         match_sub_eta[data_index][key]->Fill(cent_bin, (*jsm)->Eta());
//         match_sub_phi[data_index][key]->Fill(cent_bin, (*jsm)->Phi());
//         match_sub_pt[data_index][key]->Fill(cent_bin, (*jsm)->Pt());
//         match_sub_const[data_index][key]->Fill(cent_bin, **jsmconst);
//         match_sub_rho[data_index][key]->Fill(cent_bin, **jsmrho);
//         match_sub_sig[data_index][key]->Fill(cent_bin, **jsmsig);

//         npart[data_index][key]->Fill(cent_bin, **nprt);
//         hard_dphi[data_index][key]->Fill(cent_bin,
//         fabs((*jl)->DeltaPhi(**js)));
//         match_dphi[data_index][key]->Fill(cent_bin,
//                                           fabs((*jlm)->DeltaPhi(**jsm)));
//         lead_dr[data_index][key]->Fill(cent_bin, (*jl)->DeltaR(**jlm));
//         sub_dr[data_index][key]->Fill(cent_bin, (*js)->DeltaR(**jsm));
//         lead_dpt[data_index][key]->Fill(cent_bin, (*jl)->Pt() -
//         (*jlm)->Pt()); sub_dpt[data_index][key]->Fill(cent_bin, (*js)->Pt()
//         -
//         (*jsm)->Pt()); lead_dpt_frac[data_index][key]->Fill(
//             cent_bin, ((*jl)->Pt() - (*jlm)->Pt()) / (*jl)->Pt());
//         sub_dpt_frac[data_index][key]->Fill(
//             cent_bin, ((*js)->Pt() - (*jsm)->Pt()) / (*js)->Pt());

//         hard_aj[data_index][key]->Fill(cent_bin,
//                                        fabs((*jl)->Pt() - (*js)->Pt()) /
//                                            ((*jl)->Pt() + (*js)->Pt()));
//         match_aj[data_index][key]->Fill(cent_bin,
//                                         fabs((*jlm)->Pt() - (*jsm)->Pt()) /
//                                             ((*jlm)->Pt() + (*jsm)->Pt()));
//         hard_aj_test[data_index][key]->Fill(cent_bin,
//                                             fabs((*jl)->Pt() - (*js)->Pt())
//                                             /
//                                                 ((*jl)->Pt() +
//                                                 (*js)->Pt()));
//         match_aj_test[data_index][key]->Fill(cent_bin,
//                                              fabs((*jlm)->Pt() -
//                                              (*jsm)->Pt()) /
//                                                  ((*jlm)->Pt() +
//                                                  (*jsm)->Pt()));

//         if (auau_off_axis_present) {
//           off_axis_aj[data_index][key]->Fill(
//               cent_bin, fabs((*jloa)->Pt() - (*jsoa)->Pt()) /
//                             ((*jloa)->Pt() + (*jsoa)->Pt()));
//           off_axis_aj_test[data_index][key]->Fill(
//               cent_bin, fabs((*jloa)->Pt() - (*jsoa)->Pt()) /
//                             ((*jloa)->Pt() + (*jsoa)->Pt()));
//           off_axis_rho[data_index][key]->Fill(cent_bin, **jloarho);
//           off_axis_sig[data_index][key]->Fill(cent_bin, **jloasig);
//         }

//         if (pp_only_present) {
//           if ((*pp_only_jl)->Pt() > 0) {
//             pp_lead_hard_match_fraction[data_index][key]->Fill(cent_bin,
//             1); if ((*pp_only_js)->Pt() > 0) {
//               pp_sub_hard_match_fraction[data_index][key]->Fill(cent_bin,
//               1); hard_pp_only_aj[data_index][key]->Fill(
//                   cent_bin, fabs((*pp_only_jl)->Pt() - (*pp_only_js)->Pt())
//                   /
//                                 ((*pp_only_jl)->Pt() +
//                                 (*pp_only_js)->Pt()));
//               match_pp_only_aj[data_index][key]->Fill(
//                   cent_bin, fabs((*pp_only_jlm)->Pt() -
//                   (*pp_only_jsm)->Pt())
//                   /
//                                 ((*pp_only_jlm)->Pt() +
//                                 (*pp_only_jsm)->Pt()));
//               hard_lead_pp_dpt[data_index][key]->Fill(
//                   cent_bin, (*jl)->Pt() - (*pp_only_jl)->Pt());
//               match_lead_pp_dpt[data_index][key]->Fill(
//                   cent_bin, (*jlm)->Pt() - (*pp_only_jlm)->Pt());
//               hard_sub_pp_dpt[data_index][key]->Fill(
//                   cent_bin, (*js)->Pt() - (*pp_only_js)->Pt());
//               match_sub_pp_dpt[data_index][key]->Fill(
//                   cent_bin, (*jsm)->Pt() - (*pp_only_jsm)->Pt());
//               hard_lead_pp_dpt_frac[data_index][key]->Fill(
//                   cent_bin, ((*jl)->Pt() - (*pp_only_jl)->Pt()) /
//                   (*jl)->Pt());
//               match_lead_pp_dpt_frac[data_index][key]->Fill(
//                   cent_bin,
//                   ((*jlm)->Pt() - (*pp_only_jlm)->Pt()) / (*jlm)->Pt());
//               hard_sub_pp_dpt_frac[data_index][key]->Fill(
//                   cent_bin, ((*js)->Pt() - (*pp_only_js)->Pt()) /
//                   (*js)->Pt());
//               match_sub_pp_dpt_frac[data_index][key]->Fill(
//                   cent_bin,
//                   ((*jsm)->Pt() - (*pp_only_jsm)->Pt()) / (*jsm)->Pt());
//             } else {
//               pp_sub_hard_match_fraction[data_index][key]->Fill(cent_bin,
//               0);
//             }
//           } else
//             pp_lead_hard_match_fraction[data_index][key]->Fill(cent_bin,
//             0);
//         }
//       }
//     }
//   } // datatype

//   // now we will get systematics and print out comparisons - do everything
//   in
//   // bins of centrality
//   for (int i = 0; i < refcent_string.size(); ++i) {
//     // make our output directory
//     string out_loc = key_loc + "/cent_" + refcent_string_fs[i];
//     boost::filesystem::path dir(out_loc.c_str());
//     boost::filesystem::create_directories(dir);

//     // get systematic errors
//     systematic_errors_hard[key].push_back(
//         GetSystematic(hard_aj_cent[pp_index][key][i],
//                       hard_aj_cent[TYPEMAP[DATATYPE::TOWP]][key][i],
//                       hard_aj_cent[TYPEMAP[DATATYPE::TOWM]][key][i],
//                       hard_aj_cent[TYPEMAP[DATATYPE::TRACKP]][key][i],
//                       hard_aj_cent[TYPEMAP[DATATYPE::TRACKM]][key][i]));
//     systematic_errors_match[key].push_back(
//         GetSystematic(match_aj_cent[pp_index][key][i],
//                       match_aj_cent[TYPEMAP[DATATYPE::TOWP]][key][i],
//                       match_aj_cent[TYPEMAP[DATATYPE::TOWM]][key][i],
//                       match_aj_cent[TYPEMAP[DATATYPE::TRACKP]][key][i],
//                       match_aj_cent[TYPEMAP[DATATYPE::TRACKM]][key][i]));

//     // pave text for hard & matched jets
//     dijetcore::DijetKey params = parsed_keys[key];
//     std::stringstream streamHard;
//     std::stringstream streamMatch;
//     std::stringstream streamLeadPt;
//     std::stringstream streamSubPt;
//     std::stringstream streamConstPt;
//     streamHard << params.lead_init_r;
//     streamMatch << params.lead_match_r;
//     streamLeadPt << params.lead_init_pt;
//     streamSubPt << params.sub_init_pt;
//     streamConstPt << params.lead_init_const_pt;

//     TPaveText hardPave(0.68, 0.3, 0.88, 0.6, "NB NDC");
//     hardPave.SetFillStyle(0);
//     hardPave.SetBorderSize(0);
//     hardPave
//         .AddText(dijetcore::MakeString("Au+Au, ",
//         refcent_string[i]).c_str())
//         ->SetTextSize(0.038);
//     hardPave
//         .AddText(
//             dijetcore::MakeString("anti-k_{T}, R = ",
//             streamHard.str()).c_str())
//         ->SetTextSize(0.038);
//     hardPave
//         .AddText(dijetcore::MakeString("p_{T}^{hard const} > ",
//                                        streamConstPt.str(), "GeV/c")
//                      .c_str())
//         ->SetTextSize(0.038);
//     hardPave
//         .AddText(dijetcore::MakeString("p_{T}^{lead} > ",
//         streamLeadPt.str(),
//                                        "GeV/c")
//                      .c_str())
//         ->SetTextSize(0.038);
//     hardPave
//         .AddText(dijetcore::MakeString("p_{T}^{sublead} > ",
//         streamSubPt.str(),
//                                        "GeV/c")
//                      .c_str())
//         ->SetTextSize(0.038);
//     TPaveText matchPave(0.65, 0.3, 0.88, 0.6, "NB NDC");
//     matchPave.SetFillStyle(0);
//     matchPave.SetBorderSize(0);
//     matchPave
//         .AddText(dijetcore::MakeString("Au+Au, ",
//         refcent_string[i]).c_str())
//         ->SetTextSize(0.038);
//     matchPave
//         .AddText(dijetcore::MakeString("anti-k_{T}, R = ",
//         streamMatch.str())
//                      .c_str())
//         ->SetTextSize(0.038);
//     matchPave
//         .AddText(dijetcore::MakeString("p_{T}^{hard const} > ",
//                                        streamConstPt.str(), "GeV/c")
//                      .c_str())
//         ->SetTextSize(0.038);
//     matchPave
//         .AddText(
//             dijetcore::MakeString("p_{T}^{match const} > 0.2
//             GeV/c").c_str())
//         ->SetTextSize(0.038);

//     // print aj
//     AjPrintout(hard_aj_cent[auau_index][key][i],
//     hard_aj_cent[pp_index][key][i],
//                systematic_errors_hard[key][i], 0.0, 0.25, 0.0001, 0.9,
//                hardPave, "Au+Au HT", "p+p HT #oplus Au+Au MB", hopts,
//                copts, out_loc, "aj_hard", "", "|A_{J}|", "event fraction");
//     AjPrintout(match_aj_cent[auau_index][key][i],
//                match_aj_cent[pp_index][key][i],
//                systematic_errors_match[key][i], 0.0, 0.3, 0.0001, 0.9,
//                matchPave, "Au+Au HT", "p+p HT #oplus Au+Au MB", hopts,
//                copts, out_loc, "aj_match", "",
//                "|A_{J}|", "event fraction");
//     AjPrintout(match_aj_cent[auau_index][key][i],
//                off_axis_aj_cent[auau_index][key][i], nullptr, 0.0, 0.3,
//                0.0, 0.9, matchPave, "Au+Au HT", "Au+Au HC embedded", hopts,
//                copts, out_loc, "aj_embed", "", "|A_{J}|", "event
//                fraction");

//     // now print off-axis AJ with the matched
//     std::vector<TH1D *> off_axis_match_compare{
//         match_aj_cent[auau_index][key][i], match_aj_cent[pp_index][key][i],
//         off_axis_aj_cent[auau_index][key][i]};
//     std::vector<string> off_axis_match_string{"Au+Au", "embedded p+p",
//                                               "embedded HC Au+Au"};

//     // now we will compare all the other observables
//     Overlay1D(npart_cent[auau_index][key][i], npart_cent[pp_index][key][i],
//               "Au+Au", "p+p embedded", hopts, copts, out_loc, "npart", "",
//               "N_{part}", "event fraction");
//     Overlay1D(hard_dphi_cent[auau_index][key][i],
//               hard_dphi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "hard_dphi", "", "d#phi", "event
//               fraction");
//     Overlay1D(match_dphi_cent[auau_index][key][i],
//               match_dphi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "match_dphi", "", "d#phi", "event
//               fraction");
//     Overlay1D(lead_dr_cent[auau_index][key][i],
//     lead_dr_cent[pp_index][key][i],
//               "Au+Au", "p+p embedded", hopts, copts, out_loc, "lead_dr",
//               "", "d#phi", "event fraction");
//     Overlay1D(sub_dr_cent[auau_index][key][i],
//     sub_dr_cent[pp_index][key][i],
//               "Au+Au", "p+p embedded", hopts, copts, out_loc, "sub_dr", "",
//               "d#phi", "event fraction");
//     Overlay1D(lead_dpt_cent[auau_index][key][i],
//               lead_dpt_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "lead_dpt", "", "dp_{T}", "event
//               fraction");
//     Overlay1D(sub_dpt_cent[auau_index][key][i],
//     sub_dpt_cent[pp_index][key][i],
//               "Au+Au", "p+p embedded", hopts, copts, out_loc, "sub_dpt",
//               "", "dp_{T}", "event fraction");
//     Overlay1D(lead_dpt_frac_cent[auau_index][key][i],
//               lead_dpt_frac_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "lead_dpt_frac", "",
//               "dp_{T}", "event fraction");
//     Overlay1D(sub_dpt_frac_cent[auau_index][key][i],
//               sub_dpt_frac_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "sub_dpt_frac", "", "dp_{T}",
//               "event fraction");
//     Overlay1D(hard_lead_eta_cent[auau_index][key][i],
//               hard_lead_eta_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "hard_lead_eta", "",
//               "#eta", "event fraction");
//     Overlay1D(match_lead_eta_cent[auau_index][key][i],
//               match_lead_eta_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_lead_eta", "",
//               "#eta", "event fraction");
//     Overlay1D(hard_sub_eta_cent[auau_index][key][i],
//               hard_sub_eta_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "hard_sub_eta", "", "#eta",
//               "event fraction");
//     Overlay1D(match_sub_eta_cent[auau_index][key][i],
//               match_sub_eta_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_sub_eta", "",
//               "#eta", "event fraction");
//     Overlay1D(hard_lead_phi_cent[auau_index][key][i],
//               hard_lead_phi_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "hard_lead_phi", "",
//               "#phi", "event fraction");
//     Overlay1D(match_lead_phi_cent[auau_index][key][i],
//               match_lead_phi_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_lead_phi", "",
//               "#phi", "event fraction");
//     Overlay1D(hard_sub_phi_cent[auau_index][key][i],
//               hard_sub_phi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "hard_sub_phi", "", "#phi",
//               "event fraction");
//     Overlay1D(match_sub_phi_cent[auau_index][key][i],
//               match_sub_phi_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_sub_phi", "",
//               "#phi", "event fraction");
//     Overlay1D(hard_lead_pt_cent[auau_index][key][i],
//               hard_lead_pt_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "hard_lead_pt", "", "p_{T}",
//               "event fraction");
//     Overlay1D(match_lead_pt_cent[auau_index][key][i],
//               match_lead_pt_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_lead_pt", "",
//               "p_{T}", "event fraction");
//     Overlay1D(hard_sub_pt_cent[auau_index][key][i],
//               hard_sub_pt_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "hard_sub_pt", "", "p_{T}",
//               "event fraction");
//     Overlay1D(match_sub_pt_cent[auau_index][key][i],
//               match_sub_pt_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "match_sub_pt", "", "p_{T}",
//               "event fraction");
//     Overlay1D(hard_lead_const_cent[auau_index][key][i],
//               hard_lead_const_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "hard_lead_const", "",
//               "N_{part}", "event fraction");
//     Overlay1D(match_lead_const_cent[auau_index][key][i],
//               match_lead_const_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_lead_const", "",
//               "N_{part}", "event fraction");
//     Overlay1D(hard_sub_const_cent[auau_index][key][i],
//               hard_sub_const_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "hard_sub_const", "",
//               "N_{part}", "event fraction");
//     Overlay1D(match_sub_const_cent[auau_index][key][i],
//               match_sub_const_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_sub_const", "",
//               "N_{part}", "event fraction");
//     Overlay1D(hard_lead_rho_cent[auau_index][key][i],
//               hard_lead_rho_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "hard_lead_rho", "",
//               "#rho", "event fraction");
//     Overlay1D(match_lead_rho_cent[auau_index][key][i],
//               match_lead_rho_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_lead_rho", "",
//               "#rho", "event fraction");
//     Overlay1D(hard_sub_rho_cent[auau_index][key][i],
//               hard_sub_rho_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "hard_sub_rho", "", "#rho",
//               "event fraction");
//     Overlay1D(match_sub_rho_cent[auau_index][key][i],
//               match_sub_rho_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_sub_rho", "",
//               "#rho", "event fraction");
//     Overlay1D(hard_lead_sig_cent[auau_index][key][i],
//               hard_lead_sig_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "hard_lead_sig", "",
//               "#sigma", "event fraction");
//     Overlay1D(match_lead_sig_cent[auau_index][key][i],
//               match_lead_sig_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_lead_sig", "",
//               "#sigma", "event fraction");
//     Overlay1D(hard_sub_sig_cent[auau_index][key][i],
//               hard_sub_sig_cent[pp_index][key][i], "Au+Au", "p+p embedded",
//               hopts, copts, out_loc, "hard_sub_sig", "", "#sigma",
//               "event fraction");
//     Overlay1D(match_sub_sig_cent[auau_index][key][i],
//               match_sub_sig_cent[pp_index][key][i], "Au+Au", "p+p
//               embedded", hopts, copts, out_loc, "match_sub_sig", "",
//               "#sigma", "event fraction");

//   } // centrality

// } // key

// for (int cent = 0; cent < cent_boundaries.size(); ++cent) {
//   // make a directory for our grid outputs
//   string out_loc_grid =
//       FLAGS_outputDir + "/grid_cent_" + refcent_string_fs[cent];
//   LOG(INFO) << "output directory: " << out_loc_grid;
//   // build output directory if it doesn't exist, using boost::filesystem
//   boost::filesystem::path dir_cent(out_loc_grid.c_str());
//   boost::filesystem::create_directories(dir_cent);

//   // now build histograms
//   TH2D *p_values_hard =
//       new TH2D(dijetcore::MakeString("p_values_hard_", cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);
//   TH2D *ks_values_hard =
//       new TH2D(dijetcore::MakeString("ks_values_hard_", cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);
//   TH2D *p_values_match =
//       new TH2D(dijetcore::MakeString("p_values_match_", cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);
//   TH2D *ks_values_match =
//       new TH2D(dijetcore::MakeString("ks_values_match_", cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);
//   TH2D *p_values_bkg =
//       new TH2D(dijetcore::MakeString("p_values_bkg_", cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);
//   TH2D *ks_values_bkg =
//       new TH2D(dijetcore::MakeString("ks_values_bkg_", cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);

//   TH2D *p_values_hard_error =
//       new TH2D(dijetcore::MakeString("p_values_hard_error_", cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);
//   TH2D *ks_values_hard_error =
//       new TH2D(dijetcore::MakeString("ks_values_hard_error_",
//       cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);
//   TH2D *p_values_match_error =
//       new TH2D(dijetcore::MakeString("p_values_match_error_",
//       cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);
//   TH2D *ks_values_match_error =
//       new TH2D(dijetcore::MakeString("ks_values_match_error_",
//       cent).c_str(),
//                ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
//                constpt.size(), -0.5, constpt.size() - 0.5);

//   for (int j = 0; j < radii_string.size(); ++j) {
//     p_values_hard->GetXaxis()->SetBinLabel(j + 1, radii_string[j].c_str());
//     ks_values_hard->GetXaxis()->SetBinLabel(j + 1,
//     radii_string[j].c_str()); p_values_match->GetXaxis()->SetBinLabel(j +
//     1, radii_string[j].c_str()); ks_values_match->GetXaxis()->SetBinLabel(j
//     + 1, radii_string[j].c_str()); p_values_bkg->GetXaxis()->SetBinLabel(j
//     + 1, radii_string[j].c_str()); ks_values_bkg->GetXaxis()->SetBinLabel(j
//     + 1, radii_string[j].c_str());
//     p_values_hard_error->GetXaxis()->SetBinLabel(j + 1,
//                                                  radii_string[j].c_str());
//     ks_values_hard_error->GetXaxis()->SetBinLabel(j + 1,
//                                                   radii_string[j].c_str());
//     p_values_match_error->GetXaxis()->SetBinLabel(j + 1,
//                                                   radii_string[j].c_str());
//     ks_values_match_error->GetXaxis()->SetBinLabel(j + 1,
//                                                    radii_string[j].c_str());
//   }
//   for (int j = 0; j < constpt_string.size(); ++j) {
//     p_values_hard->GetYaxis()->SetBinLabel(j + 1,
//     constpt_string[j].c_str()); ks_values_hard->GetYaxis()->SetBinLabel(j +
//     1, constpt_string[j].c_str());
//     p_values_match->GetYaxis()->SetBinLabel(j + 1,
//     constpt_string[j].c_str()); ks_values_match->GetYaxis()->SetBinLabel(j
//     + 1, constpt_string[j].c_str());
//     p_values_bkg->GetYaxis()->SetBinLabel(j
//     + 1, constpt_string[j].c_str());
//     ks_values_bkg->GetYaxis()->SetBinLabel(j
//     + 1, constpt_string[j].c_str());
//     p_values_hard_error->GetYaxis()->SetBinLabel(j + 1,
//                                                  constpt_string[j].c_str());
//     ks_values_hard_error->GetYaxis()->SetBinLabel(j + 1,
//                                                   constpt_string[j].c_str());
//     p_values_match_error->GetYaxis()->SetBinLabel(j + 1,
//                                                   constpt_string[j].c_str());
//     ks_values_match_error->GetYaxis()->SetBinLabel(j + 1,
//                                                    constpt_string[j].c_str());
//   }

//   TCanvas *c_hard =
//       new TCanvas(dijetcore::MakeString("hard_canvas_", cent).c_str());
//   std::vector<std::vector<TPad *>> hard_pads = dijetcore::CanvasPartition(
//       c_hard, radii.size(), constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_match =
//       new TCanvas(dijetcore::MakeString("match_canvas_", cent).c_str());
//   std::vector<std::vector<TPad *>> matched_pads =
//   dijetcore::CanvasPartition(
//       c_match, radii.size(), constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_match_oa =
//       new TCanvas(dijetcore::MakeString("match_canvas_oa_", cent).c_str());
//   std::vector<std::vector<TPad *>> matched_oa_pads =
//   dijetcore::CanvasPartition(
//       c_match_oa, radii.size(), constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_hard_tow_sys_unc =
//       new TCanvas(dijetcore::MakeString("hard_tow_sys_unc_",
//       cent).c_str());
//   std::vector<std::vector<TPad *>> hard_tow_sys_unc_pads =
//       dijetcore::CanvasPartition(c_hard_tow_sys_unc, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_hard_trk_sys_unc =
//       new TCanvas(dijetcore::MakeString("hard_trk_sys_unc_pad_",
//       cent).c_str());
//   std::vector<std::vector<TPad *>> hard_trk_sys_unc_pads =
//       dijetcore::CanvasPartition(c_hard_trk_sys_unc, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_match_tow_sys_unc = new TCanvas(
//       dijetcore::MakeString("match_tow_sys_unc_pad_", cent).c_str());
//   std::vector<std::vector<TPad *>> match_tow_sys_unc_pads =
//       dijetcore::CanvasPartition(c_match_tow_sys_unc, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_match_trk_sys_unc = new TCanvas(
//       dijetcore::MakeString("match_trk_sys_unc_pad_", cent).c_str());
//   std::vector<std::vector<TPad *>> match_trk_sys_unc_pads =
//       dijetcore::CanvasPartition(c_match_trk_sys_unc, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_hard_tow_sys_unc_ratio = new TCanvas(
//       dijetcore::MakeString("hard_tow_sys_unc_ratio_pad_", cent).c_str());
//   std::vector<std::vector<TPad *>> hard_tow_sys_unc_ratio_pads =
//       dijetcore::CanvasPartition(c_hard_tow_sys_unc_ratio, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_hard_trk_sys_unc_ratio = new TCanvas(
//       dijetcore::MakeString("hard_trk_sys_unc_ratio_pad_", cent).c_str());
//   std::vector<std::vector<TPad *>> hard_trk_sys_unc_ratio_pads =
//       dijetcore::CanvasPartition(c_hard_trk_sys_unc_ratio, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_match_tow_sys_unc_ratio = new TCanvas(
//       dijetcore::MakeString("match_tow_sys_unc_ratio_pad_", cent).c_str());
//   std::vector<std::vector<TPad *>> match_tow_sys_unc_ratio_pads =
//       dijetcore::CanvasPartition(c_match_tow_sys_unc_ratio, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_match_trk_sys_unc_ratio = new TCanvas(
//       dijetcore::MakeString("match_trk_sys_unc_ratio_pad_", cent).c_str());
//   std::vector<std::vector<TPad *>> match_trk_sys_unc_ratio_pads =
//       dijetcore::CanvasPartition(c_match_trk_sys_unc_ratio, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_hard_sys_unc_ratio = new TCanvas(
//       dijetcore::MakeString("hard_sys_unc_ratio_pad_", cent).c_str());
//   std::vector<std::vector<TPad *>> hard_sys_unc_ratio_pads =
//       dijetcore::CanvasPartition(c_hard_sys_unc_ratio, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   TCanvas *c_match_sys_unc_ratio = new TCanvas(
//       dijetcore::MakeString("match_sys_unc_ratio_pad_", cent).c_str());
//   std::vector<std::vector<TPad *>> match_sys_unc_ratio_pads =
//       dijetcore::CanvasPartition(c_match_sys_unc_ratio, radii.size(),
//                                  constpt.size(), 0.10, 0.10, 0.12, 0.05);

//   // for drawing x/y axis labels
//   TPad *invis = new TPad("invis_pad", "", 0, 0, 1, 1);
//   invis->SetFillStyle(4000);
//   TPad *invis_error = new TPad("invis_err_pad", "", 0, 0, 1, 1);
//   invis_error->SetFillStyle(4000);
//   TLatex *x_label_text = new TLatex(0.47, 0.03, "|A_{J}|");
//   x_label_text->SetTextSize(0.05);
//   TLatex *y_label_text = new TLatex(0.05, 0.4, "event fraction");
//   y_label_text->SetTextSize(0.05);
//   y_label_text->SetTextAngle(90);
//   TLatex *y_label_text_error_frac = new TLatex(0.05, 0.4, "percent error");
//   y_label_text_error_frac->SetTextSize(0.05);
//   y_label_text_error_frac->SetTextAngle(90);

//   for (int rad = 0; rad < radii.size(); ++rad) {
//     for (int pt = 0; pt < constpt.size(); ++pt) {
//       int x_bin = p_values_hard->GetXaxis()->FindBin(rad);
//       int y_bin = p_values_hard->GetYaxis()->FindBin(pt);
//       string key = grid_keys[rad][pt];
//       dijetcore::DijetKey key_params = grid_key_params[rad][pt];

//       // get all histograms that will be used in producing the various
//       grids

//       TH1D *aj_hard_test = hard_aj_test_cent[auau_index][key][cent];
//       TH1D *aj_match_test = match_aj_test_cent[auau_index][key][cent];

//       TH1D *aj_hard_pp_test = hard_aj_test_cent[pp_index][key][cent];
//       TH1D *aj_match_pp_test = match_aj_test_cent[pp_index][key][cent];

//       TH1D *aj_hard = hard_aj_cent[auau_index][key][cent];
//       TH1D *aj_match = match_aj_cent[auau_index][key][cent];

//       TH1D *aj_hard_pp = hard_aj_cent[pp_index][key][cent];
//       TH1D *aj_match_pp = match_aj_cent[pp_index][key][cent];

//       TGraphErrors *hard_sys_error = systematic_errors_hard[key][cent];
//       TGraphErrors *match_sys_error = systematic_errors_match[key][cent];

//       TH1D *aj_off_axis = off_axis_aj_cent[auau_index][key][cent];
//       TH1D *aj_off_axis_test =
//       off_axis_aj_test_cent[auau_index][key][cent];

//       // including systematics, and fractional systematics
//       TH1D *hard_tow_p_aj = hard_aj_cent[towp_index][key][cent];
//       hopts.SetHistogram(hard_tow_p_aj);
//       TH1D *hard_tow_p_aj_ratio =
//           GetErrorFractional(aj_hard_pp, hard_tow_p_aj,
//                              dijetcore::MakeString("hardtowperr", cent));
//       hopts.SetHistogram(hard_tow_p_aj_ratio);

//       TH1D *hard_tow_m_aj = hard_aj_cent[towm_index][key][cent];
//       hopts.SetHistogram(hard_tow_m_aj);
//       TH1D *hard_tow_m_aj_ratio =
//           GetErrorFractional(aj_hard_pp, hard_tow_m_aj,
//                              dijetcore::MakeString("hardtowmerr", cent));
//       hopts.SetHistogram(hard_tow_m_aj_ratio);

//       TH1D *hard_trk_p_aj = hard_aj_cent[trkp_index][key][cent];
//       hopts.SetHistogram(hard_trk_p_aj);
//       TH1D *hard_trk_p_aj_ratio =
//           GetErrorFractional(aj_hard_pp, hard_trk_p_aj,
//                              dijetcore::MakeString("hardtrkperr", cent));
//       hopts.SetHistogram(hard_trk_p_aj_ratio);

//       TH1D *hard_trk_m_aj = hard_aj_cent[trkm_index][key][cent];
//       hopts.SetHistogram(hard_trk_m_aj);
//       TH1D *hard_trk_m_aj_ratio =
//           GetErrorFractional(aj_hard_pp, hard_trk_m_aj,
//                              dijetcore::MakeString("hardtrkmerr", cent));
//       hopts.SetHistogram(hard_trk_m_aj_ratio);

//       TH1D *match_tow_p_aj = match_aj_cent[towp_index][key][cent];
//       hopts.SetHistogram(match_tow_p_aj);
//       TH1D *match_tow_p_aj_ratio =
//           GetErrorFractional(aj_match_pp, match_tow_p_aj,
//                              dijetcore::MakeString("matchtowperr", cent));
//       hopts.SetHistogram(match_tow_p_aj_ratio);

//       TH1D *match_tow_m_aj = match_aj_cent[towm_index][key][cent];
//       hopts.SetHistogram(match_tow_m_aj);
//       TH1D *match_tow_m_aj_ratio =
//           GetErrorFractional(aj_match_pp, match_tow_m_aj,
//                              dijetcore::MakeString("matchtowmerr", cent));
//       hopts.SetHistogram(match_tow_m_aj_ratio);

//       TH1D *match_trk_p_aj = match_aj_cent[trkp_index][key][cent];
//       hopts.SetHistogram(match_trk_p_aj);
//       TH1D *match_trk_p_aj_ratio =
//           GetErrorFractional(aj_match_pp, match_trk_p_aj,
//                              dijetcore::MakeString("matchtrkperr", cent));
//       hopts.SetHistogram(match_trk_p_aj_ratio);

//       TH1D *match_trk_m_aj = match_aj_cent[trkm_index][key][cent];
//       hopts.SetHistogram(match_trk_m_aj);
//       TH1D *match_trk_m_aj_ratio =
//           GetErrorFractional(aj_match_pp, match_trk_m_aj,
//                              dijetcore::MakeString("matchtrkmerr", cent));
//       hopts.SetHistogram(match_trk_m_aj_ratio);

//       // // full systematic ratios for hard and matched aj
//       TH1D *hard_aj_sys =
//           GetErrorFractional(aj_hard_pp, hard_sys_error,
//                              dijetcore::MakeString("hardtotalerr", cent));
//       hopts.SetHistogram(hard_aj_sys);
//       TH1D *match_aj_sys =
//           GetErrorFractional(aj_match_pp, match_sys_error,
//                              dijetcore::MakeString("matchtotalerr", cent));
//       hopts.SetHistogram(match_aj_sys);

//       // print

//       // first, the hard aj
//       c_hard->cd(0);
//       aj_hard->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       aj_hard->GetXaxis()->SetTitleOffset(1.15);
//       aj_hard->GetXaxis()->CenterTitle();
//       aj_hard->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       aj_hard->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         aj_hard->GetYaxis()->SetTitleSize(aj_hard->GetYaxis()->GetTitleSize()
//         *
//                                           1.05);
//       if (pt == 0)
//         aj_hard->GetYaxis()->SetTitleSize(aj_hard->GetYaxis()->GetTitleSize()
//         *
//                                           0.95);
//       aj_hard->GetYaxis()->CenterTitle();
//       aj_hard->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       aj_hard->SetMarkerSize(0.2);
//       aj_hard_pp->SetMarkerSize(0.2);
//       hard_sys_error->SetFillColorAlpha(
//           hard_aj_cent[pp_index][key][cent]->GetLineColor(), 0.4);
//       hard_sys_error->SetFillStyle(1001);
//       hard_sys_error->SetLineWidth(0);
//       hard_sys_error->SetMarkerSize(0);
//       TPad *hard_pad = hard_pads[rad][pt];
//       double scale_factor =
//           hard_pads[0][0]->GetAbsHNDC() / hard_pad->GetAbsHNDC();
//       hard_pad->Draw();
//       hard_pad->SetFillStyle(4000);
//       hard_pad->SetFrameFillStyle(4000);
//       hard_pad->cd();
//       aj_hard->GetXaxis()->SetNdivisions(305);
//       aj_hard->GetXaxis()->SetLabelSize(0.08);
//       aj_hard->GetYaxis()->SetNdivisions(305);
//       aj_hard->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         aj_hard->GetYaxis()->SetLabelSize(0.11);

//       aj_hard->Draw();
//       aj_hard_pp->Draw("SAME");
//       hard_sys_error->Draw("9e2SAME");

//       // draw the STAR Prelim
//       TLegend *leg1 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg1->AddEntry(aj_hard, "Au+Au HT", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg1->AddEntry(aj_hard_pp, "p+p HT #oplus Au+Au MB", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg1->Draw();
//       c_hard->cd(0);

//       // next, plot matched aj

//       c_match->cd(0);
//       aj_match->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       aj_match->GetXaxis()->SetTitleOffset(1.15);
//       aj_match->GetXaxis()->CenterTitle();
//       aj_match->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       aj_match->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         aj_match->GetYaxis()->SetTitleSize(
//             aj_match->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         aj_match->GetYaxis()->SetTitleSize(
//             aj_match->GetYaxis()->GetTitleSize() * 0.95);
//       aj_match->GetYaxis()->CenterTitle();

//       aj_match->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       aj_match->SetMarkerSize(0.2);
//       aj_match_pp->SetMarkerSize(0.2);
//       match_sys_error->SetFillColorAlpha(
//           match_aj_cent[pp_index][key][cent]->GetLineColor(), 0.4);
//       match_sys_error->SetFillStyle(1001);
//       match_sys_error->SetLineWidth(0);
//       match_sys_error->SetMarkerSize(0);
//       TPad *match_pad = matched_pads[rad][pt];
//       match_pad->SetFillStyle(4000);
//       match_pad->SetFrameFillStyle(4000);
//       match_pad->cd();
//       aj_match->GetXaxis()->SetNdivisions(305);
//       aj_match->GetXaxis()->SetLabelSize(0.08);
//       aj_match->GetYaxis()->SetNdivisions(305);
//       aj_match->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         aj_match->GetYaxis()->SetLabelSize(0.11);

//       aj_match->Draw();
//       aj_match_pp->Draw("SAME");
//       match_sys_error->Draw("9e2SAME");

//       // draw STAR prelim
//       TLegend *leg2 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg2->AddEntry(aj_match, "Au+Au HT", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg2->AddEntry(aj_match_pp, "p+p HT #oplus Au+Au MB", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg2->Draw();
//       c_match->cd(0);

//       // next, plot off-axis aj
//       c_match_oa->cd(0);
//       aj_off_axis->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       aj_off_axis->GetXaxis()->CenterTitle();
//       aj_off_axis->GetXaxis()->SetTitleOffset(1.15);
//       aj_off_axis->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       aj_off_axis->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         aj_off_axis->GetYaxis()->SetTitleSize(
//             aj_off_axis->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         aj_off_axis->GetYaxis()->SetTitleSize(
//             aj_off_axis->GetYaxis()->GetTitleSize() * 0.95);

//       aj_off_axis->GetYaxis()->CenterTitle();

//       aj_off_axis->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       aj_off_axis->SetMarkerSize(0.2);
//       aj_off_axis->SetMarkerStyle(23);
//       aj_off_axis->SetMarkerColor(kAzure);
//       aj_off_axis->SetLineColor(kAzure);
//       aj_off_axis->SetLineWidth(2);
//       aj_off_axis->SetFillColorAlpha(kAzure, 0.65);
//       aj_off_axis->SetLineWidth(0);
//       aj_off_axis->SetMarkerSize(0);
//       TPad *match_oa_pad = matched_oa_pads[rad][pt];
//       match_oa_pad->SetFillStyle(4000);
//       match_oa_pad->SetFrameFillStyle(4000);
//       match_oa_pad->cd();
//       aj_off_axis->GetXaxis()->SetNdivisions(305);
//       aj_off_axis->GetXaxis()->SetLabelSize(0.08);
//       aj_off_axis->GetYaxis()->SetNdivisions(305);
//       aj_off_axis->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         aj_off_axis->GetYaxis()->SetLabelSize(0.11);

//       aj_off_axis->Draw("9e3");
//       aj_match->Draw("9SAME");

//       // draw STAR prelim
//       TLegend *leg3 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg3->AddEntry(aj_match, "Au+Au HT", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg3->AddEntry(aj_off_axis, "Au+Au HT #oplus Au+Au MB", "f")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg3->Draw();

//       c_match_oa->cd(0);

//       // now, tower systematic uncertainty for hard aj
//       c_hard_tow_sys_unc->cd(0);
//       hard_tow_p_aj->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       hard_tow_p_aj->GetXaxis()->SetTitleOffset(1.15);
//       hard_tow_p_aj->GetXaxis()->CenterTitle();
//       hard_tow_p_aj->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       hard_tow_p_aj->GetYaxis()->SetTitleOffset(1.07);

//       if (pt != 0 && pt != constpt.size() - 1)
//         hard_tow_p_aj->GetYaxis()->SetTitleSize(
//             hard_tow_p_aj->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         hard_tow_p_aj->GetYaxis()->SetTitleSize(
//             hard_tow_p_aj->GetYaxis()->GetTitleSize() * 0.95);
//       hard_tow_p_aj->GetYaxis()->CenterTitle();

//       hard_tow_p_aj->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       hard_tow_p_aj->SetMarkerSize(0.2);
//       hard_tow_p_aj->SetLineColor(kBlue);
//       hard_tow_p_aj->SetMarkerColor(kBlue);
//       hard_tow_m_aj->SetMarkerSize(0.2);
//       hard_tow_m_aj->SetLineColor(kMagenta);
//       hard_tow_m_aj->SetMarkerColor(kMagenta);
//       TPad *hard_sys_pad = hard_tow_sys_unc_pads[rad][pt];
//       hard_sys_pad->SetFillStyle(4000);
//       hard_sys_pad->SetFrameFillStyle(4000);
//       hard_sys_pad->cd();

//       hard_tow_p_aj->GetXaxis()->SetNdivisions(305);
//       hard_tow_p_aj->GetXaxis()->SetLabelSize(0.08);
//       hard_tow_p_aj->GetYaxis()->SetNdivisions(305);
//       hard_tow_p_aj->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         hard_tow_p_aj->GetYaxis()->SetLabelSize(0.11);

//       hard_tow_p_aj->Draw();
//       aj_hard_pp->Draw("SAME");
//       hard_tow_m_aj->Draw("SAME");

//       // draw STAR prelim
//       TLegend *leg4 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg4->AddEntry(aj_hard_pp, "p+p HT #oplus Au+Au MB", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg4->AddEntry(hard_tow_p_aj, "p+p increased tower E", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg4->AddEntry(hard_tow_m_aj, "p+p decreased tower E", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg4->Draw();

//       c_hard_tow_sys_unc->cd(0);

//       // repeat for matched tower systematic uncertainties

//       c_match_tow_sys_unc->cd(0);
//       match_tow_p_aj->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       match_tow_p_aj->GetXaxis()->SetTitleOffset(1.15);
//       match_tow_p_aj->GetXaxis()->CenterTitle();
//       match_tow_p_aj->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());

//       match_tow_p_aj->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         match_tow_p_aj->GetYaxis()->SetTitleSize(
//             match_tow_p_aj->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         match_tow_p_aj->GetYaxis()->SetTitleSize(
//             match_tow_p_aj->GetYaxis()->GetTitleSize() * 0.95);
//       match_tow_p_aj->GetYaxis()->CenterTitle();

//       match_tow_p_aj->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       match_tow_p_aj->SetMarkerSize(0.2);
//       match_tow_p_aj->SetLineColor(kBlue);
//       match_tow_p_aj->SetMarkerColor(kBlue);

//       match_tow_m_aj->SetMarkerSize(0.2);
//       match_tow_m_aj->SetLineColor(kMagenta);
//       match_tow_m_aj->SetMarkerColor(kMagenta);
//       TPad *match_tow_sys_pad = match_tow_sys_unc_pads[rad][pt];
//       match_tow_sys_pad->SetFillStyle(4000);
//       match_tow_sys_pad->SetFrameFillStyle(4000);
//       match_tow_sys_pad->cd();
//       match_tow_p_aj->GetXaxis()->SetNdivisions(305);
//       match_tow_p_aj->GetXaxis()->SetLabelSize(0.08);
//       match_tow_p_aj->GetYaxis()->SetNdivisions(305);
//       match_tow_p_aj->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         match_tow_p_aj->GetYaxis()->SetLabelSize(0.11);

//       match_tow_p_aj->Draw();
//       aj_match_pp->Draw("SAME");
//       match_tow_m_aj->Draw("SAME");

//       // draw STAR prelim
//       TLegend *leg5 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg5->AddEntry(aj_match_pp, "p+p HT #oplus Au+Au MB", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg5->AddEntry(match_tow_p_aj, "p+p increased tower E", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg5->AddEntry(match_tow_m_aj, "p+p decreased tower E", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg5->Draw();

//       c_match_tow_sys_unc->cd(0);
//       //
//       ----------------------------------------------------------------------
//       // now, tracking systematic uncertainty for hard aj

//       c_hard_trk_sys_unc->cd(0);
//       hard_trk_p_aj->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       hard_trk_p_aj->GetXaxis()->SetTitleOffset(1.15);
//       hard_trk_p_aj->GetXaxis()->CenterTitle();
//       hard_trk_p_aj->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       hard_trk_p_aj->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         hard_trk_p_aj->GetYaxis()->SetTitleSize(
//             hard_trk_p_aj->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         hard_trk_p_aj->GetYaxis()->SetTitleSize(
//             hard_trk_p_aj->GetYaxis()->GetTitleSize() * 0.95);
//       hard_trk_p_aj->GetYaxis()->CenterTitle();

//       hard_trk_p_aj->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       hard_trk_p_aj->SetMarkerSize(0.2);
//       hard_trk_p_aj->SetLineColor(kBlue);
//       hard_trk_p_aj->SetMarkerColor(kBlue);
//       hard_trk_m_aj->SetMarkerSize(0.2);
//       hard_trk_m_aj->SetLineColor(kMagenta);
//       hard_trk_m_aj->SetMarkerColor(kMagenta);
//       TPad *hard_trk_sys_pad = hard_trk_sys_unc_pads[rad][pt];
//       hard_trk_sys_pad->SetFillStyle(4000);
//       hard_trk_sys_pad->SetFrameFillStyle(4000);
//       hard_trk_sys_pad->cd();
//       hard_trk_p_aj->GetXaxis()->SetNdivisions(305);
//       hard_trk_p_aj->GetXaxis()->SetLabelSize(0.08);
//       hard_trk_p_aj->GetYaxis()->SetNdivisions(305);
//       hard_trk_p_aj->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         hard_trk_p_aj->GetYaxis()->SetLabelSize(0.11);

//       hard_trk_p_aj->Draw();
//       aj_hard_pp->Draw("SAME");
//       hard_trk_m_aj->Draw("SAME");

//       // draw STAR prelim
//       TLegend *leg6 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg6->AddEntry(aj_hard_pp, "p+p HT #oplus Au+Au MB", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg6->AddEntry(hard_trk_p_aj, "p+p increased tracking eff", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg6->AddEntry(hard_trk_m_aj, "p+p decreased tracking eff", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg6->Draw();

//       c_hard_trk_sys_unc->cd(0);

//       // repeat for matched tracking systematic uncertainties
//       c_match_trk_sys_unc->cd(0);
//       match_trk_p_aj->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       match_trk_p_aj->GetXaxis()->SetTitleOffset(1.15);
//       match_trk_p_aj->GetXaxis()->CenterTitle();
//       match_trk_p_aj->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       match_trk_p_aj->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         match_trk_p_aj->GetYaxis()->SetTitleSize(
//             match_trk_p_aj->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         match_trk_p_aj->GetYaxis()->SetTitleSize(
//             match_trk_p_aj->GetYaxis()->GetTitleSize() * 0.95);
//       match_trk_p_aj->GetYaxis()->CenterTitle();

//       match_trk_p_aj->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       match_trk_p_aj->SetMarkerSize(0.2);
//       match_trk_p_aj->SetLineColor(kBlue);
//       match_trk_p_aj->SetMarkerColor(kBlue);
//       match_trk_m_aj->SetMarkerSize(0.2);
//       match_trk_m_aj->SetLineColor(kMagenta);
//       match_trk_m_aj->SetMarkerColor(kMagenta);
//       TPad *match_trk_sys_pad = match_trk_sys_unc_pads[rad][pt];
//       match_trk_sys_pad->SetFillStyle(4000);
//       match_trk_sys_pad->SetFrameFillStyle(4000);
//       match_trk_sys_pad->cd();
//       match_trk_p_aj->GetXaxis()->SetNdivisions(305);
//       match_trk_p_aj->GetXaxis()->SetLabelSize(0.08);
//       match_trk_p_aj->GetYaxis()->SetNdivisions(305);
//       match_trk_p_aj->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         match_trk_p_aj->GetYaxis()->SetLabelSize(0.11);

//       match_trk_p_aj->Draw();
//       aj_match_pp->Draw("SAME");
//       match_trk_m_aj->Draw("SAME");

//       // draw STAR prelim
//       TLegend *leg7 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg7->AddEntry(aj_match_pp, "p+p HT #oplus Au+Au MB", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg7->AddEntry(match_trk_p_aj, "p+p increased tracking eff", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg7->AddEntry(match_trk_m_aj, "p+p decreased tracking eff", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg7->Draw();

//       c_match_trk_sys_unc->cd(0);

//       // now do systematic ratios
//       // first, hard aj tower variation
//       c_hard_tow_sys_unc_ratio->cd(0);
//       hard_tow_p_aj_ratio->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       hard_tow_p_aj_ratio->GetXaxis()->SetTitleOffset(1.15);
//       hard_tow_p_aj_ratio->GetXaxis()->CenterTitle();
//       hard_tow_p_aj_ratio->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       hard_tow_p_aj_ratio->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         hard_tow_p_aj_ratio->GetYaxis()->SetTitleSize(
//             hard_tow_p_aj_ratio->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         hard_tow_p_aj_ratio->GetYaxis()->SetTitleSize(
//             hard_tow_p_aj_ratio->GetYaxis()->GetTitleSize() * 0.95);
//       hard_tow_p_aj_ratio->GetYaxis()->CenterTitle();

//       hard_tow_p_aj_ratio->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       hard_tow_p_aj_ratio->SetMarkerSize(0.2);
//       hard_tow_p_aj_ratio->SetLineColor(kBlue);
//       hard_tow_p_aj_ratio->SetMarkerColor(kBlue);
//       hard_tow_m_aj_ratio->SetMarkerSize(0.2);
//       hard_tow_m_aj_ratio->SetLineColor(kMagenta);
//       hard_tow_m_aj_ratio->SetMarkerColor(kMagenta);
//       TPad *hard_tow_sys_ratio_pad = hard_tow_sys_unc_ratio_pads[rad][pt];
//       hard_tow_sys_ratio_pad->SetFillStyle(4000);
//       hard_tow_sys_ratio_pad->SetFrameFillStyle(4000);
//       hard_tow_sys_ratio_pad->cd();
//       hard_tow_p_aj_ratio->GetXaxis()->SetNdivisions(305);
//       hard_tow_p_aj_ratio->GetXaxis()->SetLabelSize(0.08);
//       hard_tow_p_aj_ratio->GetYaxis()->SetNdivisions(305);
//       hard_tow_p_aj_ratio->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         hard_tow_p_aj_ratio->GetYaxis()->SetLabelSize(0.11);

//       hard_tow_p_aj_ratio->Draw();
//       hard_tow_m_aj_ratio->Draw("SAME");

//       // draw STAR prelim
//       TLegend *leg8 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg8->AddEntry(hard_tow_p_aj_ratio, "p+p increased tower E", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg8->AddEntry(hard_tow_m_aj_ratio, "p+p decreased tower E", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg8->Draw();

//       c_hard_tow_sys_unc_ratio->cd(0);

//       // now, hard aj tracking variation
//       c_hard_trk_sys_unc_ratio->cd(0);

//       hard_trk_p_aj_ratio->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       hard_trk_p_aj_ratio->GetXaxis()->SetTitleOffset(1.15);
//       hard_trk_p_aj_ratio->GetXaxis()->CenterTitle();
//       hard_trk_p_aj_ratio->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       hard_trk_p_aj_ratio->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         hard_trk_p_aj_ratio->GetYaxis()->SetTitleSize(
//             hard_trk_p_aj_ratio->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         hard_trk_p_aj_ratio->GetYaxis()->SetTitleSize(
//             hard_trk_p_aj_ratio->GetYaxis()->GetTitleSize() * 0.95);
//       hard_trk_p_aj_ratio->GetYaxis()->CenterTitle();

//       hard_trk_p_aj_ratio->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       hard_trk_p_aj_ratio->SetMarkerSize(0.2);
//       hard_trk_p_aj_ratio->SetLineColor(kBlue);
//       hard_trk_p_aj_ratio->SetMarkerColor(kBlue);
//       hard_trk_m_aj_ratio->SetMarkerSize(0.2);
//       hard_trk_m_aj_ratio->SetLineColor(kMagenta);
//       hard_trk_m_aj_ratio->SetMarkerColor(kMagenta);
//       TPad *hard_trk_sys_ratio_pad = hard_trk_sys_unc_ratio_pads[rad][pt];
//       hard_trk_sys_ratio_pad->SetFillStyle(4000);
//       hard_trk_sys_ratio_pad->SetFrameFillStyle(4000);
//       hard_trk_sys_ratio_pad->cd();
//       hard_trk_p_aj_ratio->GetXaxis()->SetNdivisions(305);
//       hard_trk_p_aj_ratio->GetXaxis()->SetLabelSize(0.08);
//       hard_trk_p_aj_ratio->GetYaxis()->SetNdivisions(305);
//       hard_trk_p_aj_ratio->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         hard_trk_p_aj_ratio->GetYaxis()->SetLabelSize(0.11);

//       hard_trk_p_aj_ratio->Draw();
//       hard_trk_m_aj_ratio->Draw("SAME");

//       // draw STAR prelim
//       TLegend *leg9 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg9->AddEntry(hard_trk_p_aj_ratio, "p+p increased tracking
//       efficiency",
//                      "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg9->AddEntry(hard_trk_m_aj_ratio, "p+p decreased tracking
//       efficiency",
//                      "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg9->Draw();

//       c_hard_trk_sys_unc_ratio->cd(0);

//       // next, matched aj tower variation
//       c_match_tow_sys_unc_ratio->cd(0);

//       match_tow_p_aj_ratio->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       match_tow_p_aj_ratio->GetXaxis()->SetTitleOffset(1.15);
//       match_tow_p_aj_ratio->GetXaxis()->CenterTitle();
//       match_tow_p_aj_ratio->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       match_tow_p_aj_ratio->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         match_tow_p_aj_ratio->GetYaxis()->SetTitleSize(
//             match_tow_p_aj_ratio->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         match_tow_p_aj_ratio->GetYaxis()->SetTitleSize(
//             match_tow_p_aj_ratio->GetYaxis()->GetTitleSize() * 0.95);
//       match_tow_p_aj_ratio->GetYaxis()->CenterTitle();

//       match_tow_p_aj_ratio->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       match_tow_p_aj_ratio->SetMarkerSize(0.2);
//       match_tow_p_aj_ratio->SetLineColor(kBlue);
//       match_tow_p_aj_ratio->SetMarkerColor(kBlue);
//       match_tow_m_aj_ratio->SetMarkerSize(0.2);
//       match_tow_m_aj_ratio->SetLineColor(kMagenta);
//       match_tow_m_aj_ratio->SetMarkerColor(kMagenta);
//       TPad *match_tow_sys_ratio_pad =
//       match_tow_sys_unc_ratio_pads[rad][pt];
//       match_tow_sys_ratio_pad->SetFillStyle(4000);
//       match_tow_sys_ratio_pad->SetFrameFillStyle(4000);
//       match_tow_sys_ratio_pad->cd();
//       match_tow_p_aj_ratio->GetXaxis()->SetNdivisions(305);
//       match_tow_p_aj_ratio->GetXaxis()->SetLabelSize(0.08);
//       match_tow_p_aj_ratio->GetYaxis()->SetNdivisions(305);
//       match_tow_p_aj_ratio->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         match_tow_p_aj_ratio->GetYaxis()->SetLabelSize(0.11);

//       match_tow_p_aj_ratio->Draw();
//       match_tow_m_aj_ratio->Draw("SAME");

//       // draw STAR prelim
//       TLegend *leg10 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg10->AddEntry(match_tow_p_aj_ratio, "p+p increased tower E", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg10->AddEntry(match_tow_m_aj_ratio, "p+p decreased tower E", "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg10->Draw();

//       c_match_tow_sys_unc_ratio->cd(0);

//       // next, matched aj tracking variation
//       c_match_trk_sys_unc_ratio->cd(0);

//       match_trk_p_aj_ratio->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       match_trk_p_aj_ratio->GetXaxis()->SetTitleOffset(1.15);
//       match_trk_p_aj_ratio->GetXaxis()->CenterTitle();
//       match_trk_p_aj_ratio->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       match_trk_p_aj_ratio->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         match_trk_p_aj_ratio->GetYaxis()->SetTitleSize(
//             match_trk_p_aj_ratio->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         match_trk_p_aj_ratio->GetYaxis()->SetTitleSize(
//             match_trk_p_aj_ratio->GetYaxis()->GetTitleSize() * 0.95);
//       match_trk_p_aj_ratio->GetYaxis()->CenterTitle();

//       match_trk_p_aj_ratio->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       match_trk_p_aj_ratio->SetMarkerSize(0.2);
//       match_trk_p_aj_ratio->SetLineColor(kBlue);
//       match_trk_p_aj_ratio->SetMarkerColor(kBlue);
//       match_trk_m_aj_ratio->SetMarkerSize(0.2);
//       match_trk_m_aj_ratio->SetLineColor(kMagenta);
//       match_trk_m_aj_ratio->SetMarkerColor(kMagenta);
//       TPad *match_trk_sys_ratio_pad =
//       match_trk_sys_unc_ratio_pads[rad][pt];
//       match_trk_sys_ratio_pad->SetFillStyle(4000);
//       match_trk_sys_ratio_pad->SetFrameFillStyle(4000);
//       match_trk_sys_ratio_pad->cd();
//       match_trk_p_aj_ratio->GetXaxis()->SetNdivisions(305);
//       match_trk_p_aj_ratio->GetXaxis()->SetLabelSize(0.08);
//       match_trk_p_aj_ratio->GetYaxis()->SetNdivisions(305);
//       match_trk_p_aj_ratio->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         match_trk_p_aj_ratio->GetYaxis()->SetLabelSize(0.11);

//       match_trk_p_aj_ratio->Draw();
//       match_trk_m_aj_ratio->Draw("SAME");

//       // draw STAR prelim
//       TLegend *leg11 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg11
//           ->AddEntry(match_trk_p_aj_ratio, "p+p increased tracking
//           efficiency",
//                      "lep")
//           ->SetTextSize(.04 * scale_factor);
//       leg11
//           ->AddEntry(match_trk_m_aj_ratio, "p+p decreased tracking
//           efficiency",
//                      "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg11->Draw();

//       c_match_trk_sys_unc_ratio->cd(0);

//       // now, total systematic uncertainty for hard aj
//       c_hard_sys_unc_ratio->cd(0);

//       hard_aj_sys->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       hard_aj_sys->GetXaxis()->SetTitleOffset(1.15);
//       hard_aj_sys->GetXaxis()->CenterTitle();
//       hard_aj_sys->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       hard_aj_sys->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         hard_aj_sys->GetYaxis()->SetTitleSize(
//             hard_aj_sys->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         hard_aj_sys->GetYaxis()->SetTitleSize(
//             hard_aj_sys->GetYaxis()->GetTitleSize() * 0.95);
//       hard_aj_sys->GetYaxis()->CenterTitle();

//       hard_aj_sys->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       hard_aj_sys->SetMarkerSize(0.2);
//       hard_aj_sys->SetLineColor(kBlue);
//       hard_aj_sys->SetMarkerColor(kBlue);
//       TPad *hard_sys_ratio_pad = hard_sys_unc_ratio_pads[rad][pt];
//       hard_sys_ratio_pad->SetFillStyle(4000);
//       hard_sys_ratio_pad->SetFrameFillStyle(4000);
//       hard_sys_ratio_pad->cd();
//       hard_aj_sys->GetXaxis()->SetNdivisions(305);
//       hard_aj_sys->GetXaxis()->SetLabelSize(0.08);
//       hard_aj_sys->GetYaxis()->SetNdivisions(305);
//       hard_aj_sys->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         hard_aj_sys->GetYaxis()->SetLabelSize(0.11);

//       hard_aj_sys->Draw();

//       // draw STAR prelim
//       TLegend *leg12 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg12->AddEntry(hard_aj_sys, "p+p total systematic uncertainty",
//       "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg12->Draw();

//       c_hard_sys_unc_ratio->cd(0);

//       // total systematic uncertainty for matched aj
//       c_match_sys_unc_ratio->cd(0);

//       match_aj_sys->GetXaxis()->SetTitle(
//           dijetcore::MakeString("R=", radii_string[rad]).c_str());
//       match_aj_sys->GetXaxis()->SetTitleOffset(1.15);
//       match_aj_sys->GetXaxis()->CenterTitle();
//       match_aj_sys->GetYaxis()->SetTitle(
//           dijetcore::MakeString("p_{T}^{const}=",
//           constpt_string[pt]).c_str());
//       match_aj_sys->GetYaxis()->SetTitleOffset(1.07);
//       if (pt != 0 && pt != constpt.size() - 1)
//         match_aj_sys->GetYaxis()->SetTitleSize(
//             match_aj_sys->GetYaxis()->GetTitleSize() * 1.05);
//       if (pt == 0)
//         match_aj_sys->GetYaxis()->SetTitleSize(
//             match_aj_sys->GetYaxis()->GetTitleSize() * 0.95);
//       match_aj_sys->GetYaxis()->CenterTitle();

//       match_aj_sys->GetYaxis()->SetRangeUser(0.00001, 0.2499);
//       match_aj_sys->SetMarkerSize(0.2);
//       match_aj_sys->SetLineColor(kBlue);
//       match_aj_sys->SetMarkerColor(kBlue);
//       TPad *match_sys_ratio_pad = match_sys_unc_ratio_pads[rad][pt];
//       match_sys_ratio_pad->SetFillStyle(4000);
//       match_sys_ratio_pad->SetFrameFillStyle(4000);
//       match_sys_ratio_pad->cd();
//       match_aj_sys->GetXaxis()->SetNdivisions(305);
//       match_aj_sys->GetXaxis()->SetLabelSize(0.08);
//       match_aj_sys->GetYaxis()->SetNdivisions(305);
//       match_aj_sys->GetYaxis()->SetLabelSize(0.08);
//       if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
//         match_aj_sys->GetYaxis()->SetLabelSize(0.11);

//       match_aj_sys->Draw();

//       // draw STAR prelim
//       TLegend *leg13 =
//           GetLegend(rad, pt, radii_string.size() - 1, constpt_string.size()
//           - 1,
//                     radii[rad], constpt[pt], scale_factor);
//       leg13->AddEntry(match_aj_sys, "p+p total systematic uncertainty",
//       "lep")
//           ->SetTextSize(.04 * scale_factor);
//       if (pt == constpt.size() - 1 && rad == 0)
//         leg13->Draw();

//       c_match_sys_unc_ratio->cd(0);

//       // finally, calculate the variation in the statistical
//       // tests when varying the pp by the systematic errors, as an estimate
//       // of the robustness of the result

//       std::vector<TH1D *> aj_hard_pp_var{
//           hard_aj_cent[TYPEMAP[DATATYPE::TOWP]][key][cent],
//           hard_aj_cent[TYPEMAP[DATATYPE::TOWM]][key][cent],
//           hard_aj_cent[TYPEMAP[DATATYPE::TRACKP]][key][cent],
//           hard_aj_cent[TYPEMAP[DATATYPE::TRACKM]][key][cent]};
//       std::vector<TH1D *> aj_match_pp_var{
//           match_aj_cent[TYPEMAP[DATATYPE::TOWP]][key][cent],
//           match_aj_cent[TYPEMAP[DATATYPE::TOWM]][key][cent],
//           match_aj_cent[TYPEMAP[DATATYPE::TRACKP]][key][cent],
//           match_aj_cent[TYPEMAP[DATATYPE::TRACKM]][key][cent]};

//       std::vector<TH1D *> aj_hard_pp_var_test{
//           hard_aj_test_cent[TYPEMAP[DATATYPE::TOWP]][key][cent],
//           hard_aj_test_cent[TYPEMAP[DATATYPE::TOWM]][key][cent],
//           hard_aj_test_cent[TYPEMAP[DATATYPE::TRACKP]][key][cent],
//           hard_aj_test_cent[TYPEMAP[DATATYPE::TRACKM]][key][cent]};
//       std::vector<TH1D *> aj_match_pp_var_test{
//           match_aj_test_cent[TYPEMAP[DATATYPE::TOWP]][key][cent],
//           match_aj_test_cent[TYPEMAP[DATATYPE::TOWM]][key][cent],
//           match_aj_test_cent[TYPEMAP[DATATYPE::TRACKP]][key][cent],
//           match_aj_test_cent[TYPEMAP[DATATYPE::TRACKM]][key][cent]};

//       double p_value_hard = aj_hard->Chi2Test(aj_hard_pp, "UU NORM");
//       double p_value_match = aj_match->Chi2Test(aj_match_pp, "UU NORM");
//       double p_value_bkg = aj_match->Chi2Test(aj_off_axis, "UU NORM");
//       double ks_value_hard = aj_hard_test->KolmogorovTest(aj_hard_pp_test);
//       double ks_value_match =
//       aj_match_test->KolmogorovTest(aj_match_pp_test); double ks_value_bkg
//       = aj_match_test->KolmogorovTest(aj_off_axis_test);

//       double max_hard_p_deviation = 0;
//       double max_match_p_deviation = 0;
//       double max_hard_ks_deviation = 0;
//       double max_match_ks_deviation = 0;

//       for (int i = 0; i < aj_hard_pp_var.size(); ++i) {
//         double hard_p_deviation =
//             aj_hard->Chi2Test(aj_hard_pp_var[i], "UU NORM");
//         double hard_ks_deviation =
//             aj_hard_test->KolmogorovTest(aj_hard_pp_var_test[i]);
//         double match_p_deviation =
//             aj_match->Chi2Test(aj_match_pp_var[i], "UU NORM");
//         double match_ks_deviation =
//             aj_match_test->KolmogorovTest(aj_match_pp_var_test[i]);

//         if (fabs(p_value_hard - hard_p_deviation) > max_hard_p_deviation)
//           max_hard_p_deviation = fabs(hard_p_deviation);
//         if (fabs(ks_value_hard - hard_ks_deviation) >
//         max_hard_ks_deviation)
//           max_hard_ks_deviation = fabs(hard_ks_deviation);
//         if (fabs(p_value_match - match_p_deviation) >
//         max_match_p_deviation)
//           max_match_p_deviation = fabs(match_p_deviation);
//         if (fabs(ks_value_match - match_ks_deviation) >
//         max_match_ks_deviation)
//           max_match_ks_deviation = fabs(match_ks_deviation);
//       }

//       // set bins
//       p_values_hard->SetBinContent(x_bin, y_bin, p_value_hard);
//       p_values_match->SetBinContent(x_bin, y_bin, p_value_match);
//       p_values_bkg->SetBinContent(x_bin, y_bin, p_value_bkg);
//       ks_values_hard->SetBinContent(x_bin, y_bin, ks_value_hard);
//       ks_values_match->SetBinContent(x_bin, y_bin, ks_value_match);
//       ks_values_bkg->SetBinContent(x_bin, y_bin, ks_value_bkg);
//       p_values_hard_error->SetBinContent(x_bin, y_bin,
//       max_hard_p_deviation); p_values_match_error->SetBinContent(x_bin,
//       y_bin, max_match_p_deviation);
//       ks_values_hard_error->SetBinContent(x_bin, y_bin,
//       max_hard_ks_deviation); ks_values_match_error->SetBinContent(x_bin,
//       y_bin,
//                                            max_match_ks_deviation);
//     }
//   }

//   Print2DSimple(p_values_hard, hopts, copts, out_loc_grid, "ajphard", "",
//   "R",
//                 "p_{T}^{const}", "TEXT COLZ");
//   Print2DSimple(p_values_match, hopts, copts, out_loc_grid, "ajpmatch", "",
//   "R",
//                 "p_{T}^{const}", "TEXT COLZ");
//   Print2DSimple(p_values_bkg, hopts, copts, out_loc_grid, "ajpbkg", "",
//   "R",
//                 "p_{T}^{const}", "TEXT COLZ");
//   Print2DSimple(ks_values_hard, hopts, copts, out_loc_grid, "ajkshard", "",
//   "R",
//                 "p_{T}^{const}", "TEXT COLZ");
//   Print2DSimple(ks_values_match, hopts, copts, out_loc_grid, "ajksmatch",
//   "",
//                 "R", "p_{T}^{const}", "TEXT COLZ");
//   Print2DSimple(ks_values_bkg, hopts, copts, out_loc_grid, "ajksbkg", "",
//   "R",
//                 "p_{T}^{const}", "TEXT COLZ");

//   Print2DSimple(p_values_hard_error, hopts, copts, out_loc_grid,
//   "ajpharderror",
//                 "", "R", "p_{T}^{const}", "TEXT COLZ");
//   Print2DSimple(p_values_match_error, hopts, copts, out_loc_grid,
//                 "ajpmatcherror", "", "R", "p_{T}^{const}", "TEXT COLZ");
//   Print2DSimple(ks_values_hard_error, hopts, copts, out_loc_grid,
//                 "ajksharderror", "", "R", "p_{T}^{const}", "TEXT COLZ");
//   Print2DSimple(ks_values_match_error, hopts, copts, out_loc_grid,
//                 "ajksmatcherror", "", "R", "p_{T}^{const}", "TEXT COLZ");

//   // draw axis labels
//   c_hard->cd();
//   invis->Draw("clone");
//   invis->cd();
//   x_label_text->Draw();
//   y_label_text->Draw();
//   c_hard->SaveAs(
//       dijetcore::MakeString(out_loc_grid, "/ajhardgrid.pdf").c_str());
//   c_match->cd();
//   // draw axis labels
//   invis->Draw("clone");
//   invis->cd();
//   x_label_text->Draw();
//   y_label_text->Draw();
//   c_match->SaveAs(
//       dijetcore::MakeString(out_loc_grid, "/ajmatchgrid.pdf").c_str());
//   c_match_oa->cd();
//   // draw axis labels
//   invis->Draw("clone");
//   invis->cd();
//   x_label_text->Draw();
//   y_label_text->Draw();
//   c_match_oa->SaveAs(
//       dijetcore::MakeString(out_loc_grid,
//       "/ajmatchgridwithoa.pdf").c_str());
//   // now for systematics
//   // hard tower
//   c_hard_tow_sys_unc->cd();
//   invis->Draw("clone");
//   invis->cd();
//   x_label_text->Draw();
//   y_label_text->Draw();
//   c_hard_tow_sys_unc->SaveAs(
//       dijetcore::MakeString(out_loc_grid, "/hard_tow_sys.pdf").c_str());
//   // hard tracking
//   c_hard_trk_sys_unc->cd();
//   invis->Draw("clone");
//   invis->cd();
//   x_label_text->Draw();
//   y_label_text->Draw();
//   c_hard_trk_sys_unc->SaveAs(
//       dijetcore::MakeString(out_loc_grid, "/hard_trk_sys.pdf").c_str());
//   // match tower
//   c_match_tow_sys_unc->cd();
//   invis->Draw("clone");
//   invis->cd();
//   x_label_text->Draw();
//   y_label_text->Draw();
//   c_match_tow_sys_unc->SaveAs(
//       dijetcore::MakeString(out_loc_grid, "/match_tow_sys.pdf").c_str());
//   // match tracking
//   c_match_trk_sys_unc->cd();
//   invis->Draw("clone");
//   invis->cd();
//   x_label_text->Draw();
//   y_label_text->Draw();
//   c_match_trk_sys_unc->SaveAs(
//       dijetcore::MakeString(out_loc_grid, "/match_trk_sys.pdf").c_str());

//   // hard tower ratio
//   c_hard_tow_sys_unc_ratio->cd();
//   invis_error->Draw("clone");
//   invis_error->cd();
//   x_label_text->Draw();
//   y_label_text_error_frac->Draw();
//   c_hard_tow_sys_unc_ratio->SaveAs(
//       dijetcore::MakeString(out_loc_grid,
//       "/hard_tow_sys_ratio.pdf").c_str());
//   // hard tracking ratio
//   c_hard_trk_sys_unc_ratio->cd();
//   invis_error->Draw("clone");
//   invis_error->cd();
//   x_label_text->Draw();
//   y_label_text_error_frac->Draw();
//   c_hard_trk_sys_unc_ratio->SaveAs(
//       dijetcore::MakeString(out_loc_grid,
//       "/hard_trk_sys_ratio.pdf").c_str());
//   // match tower ratio
//   c_match_tow_sys_unc_ratio->cd();
//   invis_error->Draw("clone");
//   invis_error->cd();
//   x_label_text->Draw();
//   y_label_text_error_frac->Draw();
//   c_match_tow_sys_unc_ratio->SaveAs(
//       dijetcore::MakeString(out_loc_grid,
//       "/match_tow_sys_ratio.pdf").c_str());
//   // match tracking ratio
//   c_match_trk_sys_unc_ratio->cd();
//   invis_error->Draw("clone");
//   invis_error->cd();
//   x_label_text->Draw();
//   y_label_text_error_frac->Draw();
//   c_match_trk_sys_unc_ratio->SaveAs(
//       dijetcore::MakeString(out_loc_grid,
//       "/match_trk_sys_ratio.pdf").c_str());

//   // total hard systmatic uncertainty ratio
//   c_hard_sys_unc_ratio->cd();
//   invis_error->Draw("clone");
//   invis_error->cd();
//   x_label_text->Draw();
//   y_label_text_error_frac->Draw();
//   c_hard_sys_unc_ratio->SaveAs(
//       dijetcore::MakeString(out_loc_grid, "/hard_sys_ratio.pdf").c_str());

//   // total matched systmatic uncertainty ratio
//   c_match_sys_unc_ratio->cd();
//   invis_error->Draw("clone");
//   invis_error->cd();
//   x_label_text->Draw();
//   y_label_text_error_frac->Draw();
//   c_match_sys_unc_ratio->SaveAs(
//       dijetcore::MakeString(out_loc_grid, "/match_sys_ratio.pdf").c_str());
// } // centrality

// return 0;
// }
