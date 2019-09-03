#include "dijetcore/lib/containers.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/dijet_imbalance/generic_reader.h"
#include "dijetcore/util/dijet_imbalance/print_routines.h"
#include "dijetcore/util/dijet_imbalance/systematics.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/root/grid_print_utils.h"
#include "dijetcore/util/root/histogram_utils.h"
#include "dijetcore/util/root/root_print_utils.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
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

#include "dijetcore/lib/map.h"
#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(auau, "", "input file for Au+Au");
DIJETCORE_DEFINE_string(ppDir, "", "input directory for p+p");
DIJETCORE_DEFINE_string(outputDir, "results", "directory for output");
DIJETCORE_DEFINE_bool(
    setScanning, true,
    "fixes the initial radius, and scans through matched radii")
    DIJETCORE_DEFINE_double(initRadius, 0.2, "initial radius when set to scan");
DIJETCORE_DEFINE_string(radii, "0.2,0.25,0.3,0.35,0.4", "radii to put in grid");
DIJETCORE_DEFINE_string(constPt, "1.0,1.5,2.0,2.5,3.0", "radii to put in grid");
DIJETCORE_DEFINE_bool(useSingleCentrality, true,
                      "do things in 3 centrality bins or 1");

// enumerate the different file types we have to read in - AuAu, PP, and the 4
// systematic variations of the pp
enum DATATYPE {
  AUAU = 0,
  PP = 1,
  TOWP = 2,
  TOWM = 3,
  TRACKP = 4,
  TRACKM = 5,
  DATATYPE_SIZE = 6
};

std::unordered_map<DATATYPE, int, dijetcore::EnumClassHash> TYPEMAP{
    {DATATYPE::AUAU, static_cast<int>(DATATYPE::AUAU)},
    {DATATYPE::PP, static_cast<int>(DATATYPE::PP)},
    {DATATYPE::TOWP, static_cast<int>(DATATYPE::TOWP)},
    {DATATYPE::TOWM, static_cast<int>(DATATYPE::TOWM)},
    {DATATYPE::TRACKP, static_cast<int>(DATATYPE::TRACKP)},
    {DATATYPE::TRACKM, static_cast<int>(DATATYPE::TRACKM)}};

std::unordered_map<DATATYPE, string, dijetcore::EnumClassHash> TYPENAME{
    {DATATYPE::AUAU, "auau"},     {DATATYPE::PP, "pp"},
    {DATATYPE::TOWP, "towp"},     {DATATYPE::TOWM, "towm"},
    {DATATYPE::TRACKP, "trackp"}, {DATATYPE::TRACKM, "trackm"}};

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

  // histograms are to calculate errors by default
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  // turn off print messages
  gErrorIgnoreLevel = kInfo + 1;

  // check to make sure we have valid inputs for both AuAu and pp, because
  // the routine will break otherwise.
  std::vector<string> inputs{FLAGS_auau, FLAGS_ppDir + "/tow_0_track_0.root"};
  for (auto &file : inputs) {
    if (!boost::filesystem::exists(file)) {
      LOG(INFO) << "input file " << file;
      LOG(INFO) << "doesn't exist: exiting" << std::endl;
      return 1;
    }
  }

  TFile auau_file(FLAGS_auau.c_str(), "READ");
  TFile pp_file((FLAGS_ppDir + "/tow_0_track_0.root").c_str(), "READ");
  TFile tow_p_file((FLAGS_ppDir + "/tow_1_track_0.root").c_str(), "READ");
  TFile tow_m_file((FLAGS_ppDir + "/tow_-1_track_0.root").c_str(), "READ");
  TFile track_p_file((FLAGS_ppDir + "/tow_0_track_1.root").c_str(), "READ");
  TFile track_m_file((FLAGS_ppDir + "/tow_0_track_-1.root").c_str(), "READ");

  // if the output directory doesn't exist, we create it, since ROOT will not
  // create new directories for output
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);

  // define centralities - single centrality is the y7 0-20%,
  // otherwise, use 0-5, 5-10, 10-20%
  std::vector<unsigned> cent_boundaries;
  std::vector<std::pair<int, int>> cent_bin_boundaries;
  std::vector<string> refcent_string;
  std::vector<string> refcent_string_fs;

  if (FLAGS_useSingleCentrality) {
    cent_boundaries = {269};
    cent_bin_boundaries = {{0, 2}};
    refcent_string = {"0-20%"};
    refcent_string_fs = {"0_20"};
  } else {
    cent_boundaries = {485, 399, 269};
    cent_bin_boundaries = {{0, 0}, {1, 1}, {2, 2}};
    refcent_string = {"0-5%", "5-10%", "10-20%"};
    refcent_string_fs = {"0_5", "5_10", "10_20"};
  }

  // we need to parse the grid values to floating point so that they can be
  // compared to the key values
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
  double max_eta = 1.0 - max_radii;

  LOG(INFO) << "maximum radius = " << max_radii;
  LOG(INFO) << "so our maximum jet eta = " << max_eta;

  if (max_eta < 0 || max_eta > 1.0) {
    LOG(ERROR) << "something weird is going on with max jet radius...";
    return 1;
  }

  // now we'll get the trees from the files, ignoring any objects
  // in the file that don't conform to the naming conventions from
  // the DijetWorker.
  std::vector<string> keys;
  std::unordered_map<std::string, dijetcore::DijetKey> parsed_keys;
  std::vector<std::unordered_map<string, TTree *>> trees(TYPEMAP.size());

  GetTreesFromFile(auau_file, trees[TYPEMAP[DATATYPE::AUAU]]);
  GetTreesFromFile(pp_file, trees[TYPEMAP[DATATYPE::PP]]);
  GetTreesFromFile(tow_p_file, trees[TYPEMAP[DATATYPE::TOWP]]);
  GetTreesFromFile(tow_m_file, trees[TYPEMAP[DATATYPE::TOWM]]);
  GetTreesFromFile(track_p_file, trees[TYPEMAP[DATATYPE::TRACKP]]);
  GetTreesFromFile(track_m_file, trees[TYPEMAP[DATATYPE::TRACKM]]);

  // match keys
  for (auto entry : trees[AUAU]) {
    if (trees[PP].find(entry.first) != trees[PP].end()) {
      keys.push_back(entry.first);
      parsed_keys[entry.first] = dijetcore::ParseStringToDijetKey(entry.first);
    }
  }

  LOG(INFO) << "number of matched keys found: " << keys.size();

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

  // save all histograms so we can do comparisons between different keys if we
  // want

  const int num_datatypes = TYPEMAP.size();

  std::vector<std::unordered_map<string, TH2D *>> hard_aj(num_datatypes);
  std::vector<std::unordered_map<string, TH2D *>> match_aj(num_datatypes);
  std::vector<std::unordered_map<string, TH2D *>> hard_aj_test(num_datatypes);
  std::vector<std::unordered_map<string, TH2D *>> match_aj_test(num_datatypes);

  std::vector<std::unordered_map<string, TH2D *>> hard_aj_full(num_datatypes);
  std::vector<std::unordered_map<string, TH2D *>> match_aj_full(num_datatypes);
  std::vector<std::unordered_map<string, TH2D *>> hard_aj_full_test(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D *>> match_aj_full_test(
      num_datatypes);

  std::vector<std::unordered_map<string, TH2D *>> off_axis_aj(num_datatypes);
  std::vector<std::unordered_map<string, TH2D *>> off_axis_aj_test(
      num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D *>>> hard_aj_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D *>>> match_aj_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D *>>>
      hard_aj_test_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D *>>>
      match_aj_test_cent(num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D *>>>
      hard_aj_full_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D *>>>
      match_aj_full_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D *>>>
      hard_aj_full_test_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D *>>>
      match_aj_full_test_cent(num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D *>>> off_axis_aj_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D *>>>
      off_axis_aj_test_cent(num_datatypes);

  std::unordered_map<string, std::vector<TGraphErrors *>>
      systematic_errors_hard;
  std::unordered_map<string, std::vector<TGraphErrors *>>
      systematic_errors_match;

  std::unordered_map<string, std::vector<TGraphErrors *>>
      systematic_errors_hard_full;
  std::unordered_map<string, std::vector<TGraphErrors *>>
      systematic_errors_match_full;

  std::unordered_map<string, std::vector<TH1D *>> hard_aj_cent_err_frac(
      num_datatypes);
  std::unordered_map<string, std::vector<TH1D *>> match_aj_cent_err_frac(
      num_datatypes);

  // slightly different structure for the ks test results - the single index is
  // centrality, while each position in the 2D histogram corresponds to a single
  // key
  std::vector<TH2D *> ks_hard;
  std::vector<TH2D *> ks_match;
  std::vector<TH2D *> ks_hard_err;
  std::vector<TH2D *> ks_hard_rel_err;
  std::vector<TH2D *> ks_match_err;
  std::vector<TH2D *> ks_match_rel_err;
  std::vector<TH2D *> ks_oa;
  std::vector<TH2D *> ks_hard_full;
  std::vector<TH2D *> ks_match_full;

  for (int cent = 0; cent < refcent_string.size(); ++cent) {
    ks_hard.push_back(new TH2D(
        dijetcore::MakeString("ks_values_hard_", refcent_string[cent]).c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
    ks_match.push_back(new TH2D(
        dijetcore::MakeString("ks_values_match_", refcent_string[cent]).c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
    ks_hard_err.push_back(new TH2D(
        dijetcore::MakeString("ks_values_hard_err_", refcent_string[cent])
            .c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
    ks_hard_rel_err.push_back(new TH2D(
        dijetcore::MakeString("ks_values_hard_rel_err_", refcent_string[cent])
            .c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
    ks_match_err.push_back(new TH2D(
        dijetcore::MakeString("ks_values_match_err_", refcent_string[cent])
            .c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
    ks_match_rel_err.push_back(new TH2D(
        dijetcore::MakeString("ks_values_match_rel_err_", refcent_string[cent])
            .c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
    ks_oa.push_back(new TH2D(
        dijetcore::MakeString("ks_values_oa_", refcent_string[cent]).c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
    ks_hard_full.push_back(new TH2D(
        dijetcore::MakeString("ks_values_hard_full_", refcent_string[cent])
            .c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
    ks_match_full.push_back(new TH2D(
        dijetcore::MakeString("ks_values_match_full_", refcent_string[cent])
            .c_str(),
        ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
        constpt.size(), -0.5, constpt.size() - 0.5));
  }

  // create our histogram and canvas options
  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;
  copts.leg_left_bound = 0.6;
  dijetcore::canvasOpts coptslogz;
  coptslogz.leg_left_bound = 0.6;
  coptslogz.log_z = true;
  dijetcore::canvasOpts coptslogy;
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

  // start loop over all keys
  // start entry at -1 so that the indexing starts at 0 :)
  int entry = -1;
  for (auto &key : keys) {
    // increment entry counter
    entry++;
    LOG(INFO) << "loading: " << key;

    // make our output directory (again, ROOT won't do it for us)
    string key_loc = FLAGS_outputDir + "/" + key;
    boost::filesystem::path dir(key_loc.c_str());
    boost::filesystem::create_directories(dir);

    // read in each datatype individually, and calculate their distributions
    for (auto type : TYPEMAP) {
      auto &type_enum = type.first;
      auto &data_index = type.second;

      string out_loc = key_loc + "/" + TYPENAME[type_enum];
      boost::filesystem::path dir(out_loc.c_str());
      boost::filesystem::create_directories(dir);

      TTree *current_tree = trees[data_index][key];

      // check if we have off-axis entries (for auau)
      bool auau_off_axis_present = false;
      if (data_index == TYPEMAP[DATATYPE::AUAU] &&
          current_tree->GetBranch("jloa"))
        auau_off_axis_present = true;

      // check if we have pp-only cluster results
      bool pp_only_present = false;
      if (data_index == TYPEMAP[DATATYPE::PP] &&
          current_tree->GetBranch("foundpp"))
        pp_only_present = true;

      // and check if we have embedding (for pp tree)
      bool pp_embedded_present = false;
      if (data_index == TYPEMAP[DATATYPE::PP] &&
          current_tree->GetBranch("embed_eventid"))
        pp_embedded_present = true;

      dijetcore::GenericReader reader(current_tree, auau_off_axis_present,
                                      pp_only_present, pp_embedded_present);

      // build prefix for histogram names so ROOT doesn't get confused
      // (root fuckin sucks)
      std::string key_prefix = "key_" + std::to_string(entry) + "_";
      std::string datatype_prefix = TYPENAME[type_enum];
      std::string hist_prefix = key_prefix + datatype_prefix;

      hard_aj[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardaj").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   15, 0.0001, 0.9);
      match_aj[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchaj").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   15, 0.0001, 0.9);
      hard_aj_test[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardajtest").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   10000, 0.0001, 0.9);
      match_aj_test[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchajtest").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 10000, 0.0001, 0.9);

      hard_aj_full[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardajfull").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   25, -0.4999, 0.9);
      match_aj_full[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchajfull").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 25, -0.4999, 0.9);
      hard_aj_full_test[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardajfulltest").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 10000, -0.4999, 0.9);
      match_aj_full_test[data_index][key] = new TH2D(
          dijetcore::MakeString(hist_prefix, "matchajfulltest").c_str(), "",
          cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000,
          -0.4999, 0.9);

      if (auau_off_axis_present) {
        off_axis_aj[data_index][key] =
            new TH2D(dijetcore::MakeString(hist_prefix, "offaxisaj").c_str(),
                     "", cent_boundaries.size(), -0.5,
                     cent_boundaries.size() - 0.5, 15, 0.0001, 0.9);
        off_axis_aj_test[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "offaxisajtest").c_str(), "",
            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000,
            0.0001, 0.9);
      }

      // event loop
      while (reader.Next()) {
        int cent_bin = -1;
        for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
          auto &pair = cent_bin_boundaries[i];
          if (reader.cent >= pair.first && reader.cent <= pair.second) {
            cent_bin = i;
            break;
          }
        }

        if (cent_bin == -1)
          continue;

        // check jet eta
        if (FLAGS_setScanning && (fabs(reader.jl->Eta()) > max_eta ||
                                  fabs(reader.js->Eta()) > max_eta))
          continue;

        double jl_pt = reader.jl->Pt();
        double jlm_pt = reader.jlm->Pt();
        double js_pt = reader.js->Pt();
        double jsm_pt = reader.jsm->Pt();

        hard_aj[data_index][key]->Fill(cent_bin,
                                       fabs((jl_pt - js_pt) / (jl_pt + js_pt)));
        match_aj[data_index][key]->Fill(
            cent_bin, fabs((jlm_pt - jsm_pt) / (jlm_pt + jsm_pt)));
        hard_aj_test[data_index][key]->Fill(
            cent_bin, fabs((jl_pt - js_pt) / (jl_pt + js_pt)));
        match_aj_test[data_index][key]->Fill(
            cent_bin, fabs((jlm_pt - jsm_pt) / (jlm_pt + jsm_pt)));

        hard_aj_full[data_index][key]->Fill(cent_bin,
                                            (jl_pt - js_pt) / (jl_pt + js_pt));
        match_aj_full[data_index][key]->Fill(cent_bin, (jlm_pt - jsm_pt) /
                                                           (jlm_pt + jsm_pt));
        hard_aj_full_test[data_index][key]->Fill(cent_bin, (jl_pt - js_pt) /
                                                               (jl_pt + js_pt));
        match_aj_full_test[data_index][key]->Fill(
            cent_bin, (jlm_pt - jsm_pt) / (jlm_pt + jsm_pt));

        if (auau_off_axis_present) {
          double jloa_pt = reader.jloa->Pt();
          double jsoa_pt = reader.jsoa->Pt();

          off_axis_aj[data_index][key]->Fill(
              cent_bin, fabs((jloa_pt - jsoa_pt) / (jloa_pt + jsoa_pt)));
          off_axis_aj_test[data_index][key]->Fill(
              cent_bin, fabs((jloa_pt - jsoa_pt) / (jloa_pt + jsoa_pt)));
        }
      }

      hard_aj_cent[data_index][key] =
          dijetcore::Split2DByBinNormalized(hard_aj[data_index][key]);
      match_aj_cent[data_index][key] =
          dijetcore::Split2DByBinNormalized(match_aj[data_index][key]);
      hard_aj_test_cent[data_index][key] =
          dijetcore::Split2DByBinNormalized(hard_aj_test[data_index][key]);
      match_aj_test_cent[data_index][key] =
          dijetcore::Split2DByBinNormalized(match_aj_test[data_index][key]);

      hard_aj_full_cent[data_index][key] =
          dijetcore::Split2DByBinNormalized(hard_aj_full[data_index][key]);
      match_aj_full_cent[data_index][key] =
          dijetcore::Split2DByBinNormalized(match_aj_full[data_index][key]);
      hard_aj_full_test_cent[data_index][key] =
          dijetcore::Split2DByBinNormalized(hard_aj_full_test[data_index][key]);
      match_aj_full_test_cent[data_index][key] =
          dijetcore::Split2DByBinNormalized(
              match_aj_full_test[data_index][key]);

      if (auau_off_axis_present) {
        off_axis_aj_cent[data_index][key] =
            dijetcore::Split2DByBinNormalized(off_axis_aj[data_index][key]);
        off_axis_aj_test_cent[data_index][key] =
            dijetcore::Split2DByBinNormalized(
                off_axis_aj_test[data_index][key]);
        Overlay1D(off_axis_aj_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "off_axis_aj", "", "|A_{J}|",
                  "event fraction", "Centrality");
      }
    } // datatype

    // now we will get systematics and print out comparisons - do everything in
    // bins of centrality
    LOG(INFO) << "printing systematics";
    for (int i = 0; i < refcent_string.size(); ++i) {
      // make our output directory
      string out_loc = key_loc + "/cent_" + refcent_string_fs[i];
      boost::filesystem::path dir(out_loc.c_str());
      boost::filesystem::create_directories(dir);

      // get systematic errors by comparing the PP nominal result to the 4
      // systematic variations
      systematic_errors_hard[key].push_back(dijetcore::GetSystematic(
          hard_aj_cent[PP][key][i], hard_aj_cent[TOWP][key][i],
          hard_aj_cent[TOWM][key][i], hard_aj_cent[TRACKP][key][i],
          hard_aj_cent[TRACKM][key][i]));
      systematic_errors_match[key].push_back(dijetcore::GetSystematic(
          match_aj_cent[PP][key][i], match_aj_cent[TOWP][key][i],
          match_aj_cent[TOWM][key][i], match_aj_cent[TRACKP][key][i],
          match_aj_cent[TRACKM][key][i]));

      systematic_errors_hard_full[key].push_back(dijetcore::GetSystematic(
          hard_aj_full_cent[PP][key][i], hard_aj_full_cent[TOWP][key][i],
          hard_aj_full_cent[TOWM][key][i], hard_aj_full_cent[TRACKP][key][i],
          hard_aj_full_cent[TRACKM][key][i]));
      systematic_errors_match_full[key].push_back(dijetcore::GetSystematic(
          match_aj_full_cent[PP][key][i], match_aj_full_cent[TOWP][key][i],
          match_aj_full_cent[TOWM][key][i], match_aj_full_cent[TRACKP][key][i],
          match_aj_full_cent[TRACKM][key][i]));

      // fill fractional errors
      hard_aj_cent_err_frac[key].push_back(dijetcore::GetFractionalError(
          hard_aj_cent[PP][key][i], systematic_errors_hard[key][i]));
      match_aj_cent_err_frac[key].push_back(dijetcore::GetFractionalError(
          match_aj_cent[PP][key][i], systematic_errors_match[key][i]));

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
          .AddText(dijetcore::MakeString("Au+Au, ", refcent_string[i]).c_str())
          ->SetTextSize(0.038);
      hardPave
          .AddText(dijetcore::MakeString("anti-k_{T}, R = ", streamHard.str())
                       .c_str())
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
          .AddText(dijetcore::MakeString("p_{T}^{sublead} > ",
                                         streamSubPt.str(), "GeV/c")
                       .c_str())
          ->SetTextSize(0.038);
      TPaveText matchPave(0.65, 0.3, 0.88, 0.6, "NB NDC");
      matchPave.SetFillStyle(0);
      matchPave.SetBorderSize(0);
      matchPave
          .AddText(dijetcore::MakeString("Au+Au, ", refcent_string[i]).c_str())
          ->SetTextSize(0.038);
      matchPave
          .AddText(dijetcore::MakeString("anti-k_{T}, R = ", streamMatch.str())
                       .c_str())
          ->SetTextSize(0.038);
      matchPave
          .AddText(dijetcore::MakeString("p_{T}^{hard const} > ",
                                         streamConstPt.str(), "GeV/c")
                       .c_str())
          ->SetTextSize(0.038);
      matchPave
          .AddText(
              dijetcore::MakeString("p_{T}^{match const} > 0.2 GeV/c").c_str())
          ->SetTextSize(0.038);

      // print aj
      dijetcore::AjPrintout(
          hard_aj_cent[AUAU][key][i], hard_aj_cent[PP][key][i],
          systematic_errors_hard[key][i], 0.0, 0.25, 0.0001, 0.9, hardPave,
          "Au+Au HT", "p+p HT #oplus Au+Au MB", hopts, copts, out_loc,
          "aj_hard", "", "|A_{J}|", "event fraction");
      dijetcore::AjPrintout(
          match_aj_cent[AUAU][key][i], match_aj_cent[PP][key][i],
          systematic_errors_match[key][i], 0.0, 0.3, 0.0001, 0.9, matchPave,
          "Au+Au HT", "p+p HT #oplus Au+Au MB", hopts, copts, out_loc,
          "aj_match", "", "|A_{J}|", "event fraction");
      if (off_axis_aj_cent[AUAU].count(key)) {
        dijetcore::AjPrintout(match_aj_cent[AUAU][key][i],
                              off_axis_aj_cent[AUAU][key][i], nullptr, 0.0, 0.3,
                              0.0, 0.9, matchPave, "Au+Au HT",
                              "Au+Au HC embedded", hopts, copts, out_loc,
                              "aj_embed", "", "|A_{J}|", "event fraction");
      }
      dijetcore::AjPrintout(
          hard_aj_full_cent[AUAU][key][i], hard_aj_full_cent[PP][key][i],
          systematic_errors_hard_full[key][i], 0.0, 0.25, -0.4999, 0.9,
          hardPave, "Au+Au HT", "p+p HT #oplus Au+Au MB", hopts, copts, out_loc,
          "aj_hard_full", "", "A_{J}", "event fraction");
      dijetcore::AjPrintout(
          match_aj_full_cent[AUAU][key][i], match_aj_full_cent[PP][key][i],
          systematic_errors_match_full[key][i], 0.0, 0.3, -0.4999, 0.9,
          matchPave, "Au+Au HT", "p+p HT #oplus Au+Au MB", hopts, copts,
          out_loc, "aj_match_full", "", "A_{J}", "event fraction");

      // run statistical tests - get ks values for each
      double ks_hard_val = hard_aj_test_cent[AUAU][key][i]->KolmogorovTest(
          hard_aj_test_cent[PP][key][i]);
      double ks_match_val = match_aj_test_cent[AUAU][key][i]->KolmogorovTest(
          match_aj_test_cent[PP][key][i]);
      double ks_oa_val = 0;
      if (off_axis_aj[AUAU][key])
        ks_oa_val = match_aj_test_cent[AUAU][key][i]->KolmogorovTest(
            off_axis_aj_test_cent[AUAU][key][i]);
      double ks_hard_full_val =
          hard_aj_full_test_cent[AUAU][key][i]->KolmogorovTest(
              hard_aj_full_test_cent[PP][key][i]);
      double ks_match_full_val =
          match_aj_full_test_cent[AUAU][key][i]->KolmogorovTest(
              match_aj_full_test_cent[PP][key][i]);

      // get the error on the ks test value for hard and matched di-jets by
      // comparing AuAu to the systematic variations of PP
      double max_hard_variation = 0;
      double max_match_variation = 0;
      double max_hard_rel_var = 0.0;
      double max_match_rel_var = 0.0;

      for (int var = TOWP; var < DATATYPE_SIZE; ++var) {
        double var_hard = hard_aj_test_cent[AUAU][key][i]->KolmogorovTest(
            hard_aj_test_cent[var][key][i]);
        double var_match = match_aj_test_cent[AUAU][key][i]->KolmogorovTest(
            match_aj_test_cent[var][key][i]);

        if (max_hard_variation == 0.0 ||
            fabs(var_hard - ks_hard_val) >
                fabs(max_hard_variation - ks_hard_val)) {
          max_hard_variation = var_hard;
          max_hard_rel_var = (var_hard - ks_hard_val) / ks_hard_val;
        }
        if (max_match_variation == 0.0 ||
            fabs(var_match - ks_match_val) >
                fabs(max_match_variation - ks_match_val)) {
          max_match_variation = var_match;
          max_match_rel_var = (var_match - ks_match_val) / ks_match_val;
        }
      }

      // get the parsed key
      dijetcore::DijetKey key_vals = parsed_keys[key];

      size_t rad_idx = dijetcore::FindFirst(radii, key_vals.lead_match_r);
      size_t pt_idx =
          dijetcore::FindFirst(constpt, key_vals.lead_init_const_pt);

      int x_bin = ks_hard[i]->GetXaxis()->FindBin(rad_idx);
      int y_bin = ks_hard[i]->GetYaxis()->FindBin(pt_idx);

      ks_hard[i]->SetBinContent(x_bin, y_bin, ks_hard_val);
      ks_match[i]->SetBinContent(x_bin, y_bin, ks_match_val);
      ks_oa[i]->SetBinContent(x_bin, y_bin, ks_oa_val);
      ks_hard_err[i]->SetBinContent(x_bin, y_bin, max_hard_variation);
      ks_hard_rel_err[i]->SetBinContent(x_bin, y_bin, max_hard_rel_var);
      ks_match_err[i]->SetBinContent(x_bin, y_bin, max_match_variation);
      ks_match_rel_err[i]->SetBinContent(x_bin, y_bin, max_match_rel_var);
      ks_hard_full[i]->SetBinContent(x_bin, y_bin, ks_hard_full_val);
      ks_match_full[i]->SetBinContent(x_bin, y_bin, ks_match_full_val);

    } // centrality

  } // key

  // ------------------------------------------------------------------------------------
  // printout the ks test results
  LOG(INFO) << "print test results";
  for (int cent = 0; cent < ks_hard.size(); ++cent) {

    boost::filesystem::path test_dir(FLAGS_outputDir);
    test_dir /= "test_results";
    test_dir /= refcent_string[cent];
    boost::filesystem::create_directories(test_dir);

    string out_dir = test_dir.string();

    for (int rad = 0; rad < radii_string.size(); ++rad) {
      ks_hard[cent]->GetXaxis()->SetBinLabel(rad + 1,
                                             radii_string[rad].c_str());
      ks_match[cent]->GetXaxis()->SetBinLabel(rad + 1,
                                              radii_string[rad].c_str());
      ks_oa[cent]->GetXaxis()->SetBinLabel(rad + 1, radii_string[rad].c_str());
      ks_hard_err[cent]->GetXaxis()->SetBinLabel(rad + 1,
                                                 radii_string[rad].c_str());
      ks_match_err[cent]->GetXaxis()->SetBinLabel(rad + 1,
                                                  radii_string[rad].c_str());
      ks_hard_full[cent]->GetXaxis()->SetBinLabel(rad + 1,
                                                  radii_string[rad].c_str());
      ks_match_full[cent]->GetXaxis()->SetBinLabel(rad + 1,
                                                   radii_string[rad].c_str());
    }

    for (int pt = 0; pt < constpt_string.size(); ++pt) {
      ks_hard[cent]->GetYaxis()->SetBinLabel(pt + 1,
                                             constpt_string[pt].c_str());
      ks_match[cent]->GetYaxis()->SetBinLabel(pt + 1,
                                              constpt_string[pt].c_str());
      ks_oa[cent]->GetYaxis()->SetBinLabel(pt + 1, constpt_string[pt].c_str());
      ks_hard_err[cent]->GetYaxis()->SetBinLabel(pt + 1,
                                                 constpt_string[pt].c_str());
      ks_match_err[cent]->GetYaxis()->SetBinLabel(pt + 1,
                                                  constpt_string[pt].c_str());
      ks_hard_full[cent]->GetYaxis()->SetBinLabel(pt + 1,
                                                  constpt_string[pt].c_str());
      ks_match_full[cent]->GetYaxis()->SetBinLabel(pt + 1,
                                                   constpt_string[pt].c_str());
    }
    LOG(INFO) << "printing";
    Print2DSimple(ks_hard[cent], hopts, copts, out_dir, "ks_hard", "", "R",
                  "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_match[cent], hopts, copts, out_dir, "ks_match", "", "R",
                  "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_oa[cent], hopts, copts, out_dir, "ks_oa", "", "R",
                  "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_hard_err[cent], hopts, copts, out_dir, "ks_hard_err", "",
                  "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_hard_rel_err[cent], hopts, copts, out_dir,
                  "ks_hard_rel_err", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_match_err[cent], hopts, copts, out_dir, "ks_match_err", "",
                  "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_match_rel_err[cent], hopts, copts, out_dir,
                  "ks_match_rel_err", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_hard_full[cent], hopts, copts, out_dir, "ks_hard_full", "",
                  "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_match_full[cent], hopts, copts, out_dir, "ks_match_full",
                  "", "R", "p_{T}^{const}", "TEXT COLZ");
  }

  // ------------------------------------------------------------------------------------
  // with all keys run, we can now start to print the grids

  // each point on the grid is its own TPad - so first we create the TCanvases
  // and split them
  double grid_x_low_margin = 0.10;
  double grid_x_high_margin = 0.10;
  double grid_y_low_margin = 0.10;
  double grid_y_high_margin = 0.10;

  LOG(INFO) << "creating pads";
  std::vector<TCanvas *> canvas_hard;
  std::vector<std::vector<std::vector<TPad *>>> hard_pads;
  std::vector<TCanvas *> canvas_match;
  std::vector<std::vector<std::vector<TPad *>>> match_pads;
  std::vector<TCanvas *> canvas_hard_full;
  std::vector<std::vector<std::vector<TPad *>>> hard_pads_full;
  std::vector<TCanvas *> canvas_match_full;
  std::vector<std::vector<std::vector<TPad *>>> match_pads_full;
  std::vector<TCanvas *> canvas_oa;
  std::vector<std::vector<std::vector<TPad *>>> oa_pads;
  std::vector<TCanvas *> canvas_hard_err_frac;
  std::vector<std::vector<std::vector<TPad *>>> hard_err_frac_pads;
  std::vector<TCanvas *> canvas_match_err_frac;
  std::vector<std::vector<std::vector<TPad *>>> match_err_frac_pads;
  std::vector<TCanvas *> canvas_hard_sys;
  std::vector<std::vector<std::vector<TPad *>>> hard_sys_pads;
  std::vector<TCanvas *> canvas_match_sys;
  std::vector<std::vector<std::vector<TPad *>>> match_sys_pads;

  for (int cent = 0; cent < refcent_string.size(); ++cent) {
    canvas_hard.push_back(
        new TCanvas(dijetcore::MakeString("hard_canvas_", cent).c_str()));
    hard_pads.push_back(dijetcore::CanvasPartition(canvas_hard[cent],
                                                   radii.size(), constpt.size(),
                                                   0.10, 0.10, 0.10, 0.10));
    canvas_match.push_back(
        new TCanvas(dijetcore::MakeString("match_canvas_", cent).c_str()));
    match_pads.push_back(
        dijetcore::CanvasPartition(canvas_match[cent], radii.size(),
                                   constpt.size(), 0.10, 0.10, 0.10, 0.10));
    canvas_hard_full.push_back(
        new TCanvas(dijetcore::MakeString("hard_canvas_full_", cent).c_str()));
    hard_pads_full.push_back(
        dijetcore::CanvasPartition(canvas_hard_full[cent], radii.size(),
                                   constpt.size(), 0.10, 0.10, 0.10, 0.10));
    canvas_match_full.push_back(
        new TCanvas(dijetcore::MakeString("match_canvas_full_", cent).c_str()));
    match_pads_full.push_back(
        dijetcore::CanvasPartition(canvas_match_full[cent], radii.size(),
                                   constpt.size(), 0.10, 0.10, 0.10, 0.10));
    canvas_oa.push_back(
        new TCanvas(dijetcore::MakeString("oa_canvas_", cent).c_str()));
    oa_pads.push_back(dijetcore::CanvasPartition(
        canvas_oa[cent], radii.size(), constpt.size(), 0.10, 0.10, 0.10, 0.10));
    canvas_hard_err_frac.push_back(new TCanvas(
        dijetcore::MakeString("hard_err_frac_canvas_", cent).c_str()));
    hard_err_frac_pads.push_back(
        dijetcore::CanvasPartition(canvas_hard_err_frac[cent], radii.size(),
                                   constpt.size(), 0.10, 0.10, 0.10, 0.10));
    canvas_match_err_frac.push_back(new TCanvas(
        dijetcore::MakeString("match_err_frac_canvas_", cent).c_str()));
    match_err_frac_pads.push_back(
        dijetcore::CanvasPartition(canvas_match_err_frac[cent], radii.size(),
                                   constpt.size(), 0.10, 0.10, 0.10, 0.10));
    canvas_hard_sys.push_back(
        new TCanvas(dijetcore::MakeString("hard_sys_canvas_", cent).c_str()));
    hard_sys_pads.push_back(
        dijetcore::CanvasPartition(canvas_hard_sys[cent], radii.size(),
                                   constpt.size(), 0.10, 0.10, 0.10, 0.10));
    canvas_match_sys.push_back(
        new TCanvas(dijetcore::MakeString("match_sys_canvas_", cent).c_str()));
    match_sys_pads.push_back(
        dijetcore::CanvasPartition(canvas_match_sys[cent], radii.size(),
                                   constpt.size(), 0.10, 0.10, 0.10, 0.10));
  }

  // ------------------------------------------------------------------------------------
  // have to build the 2D vectors of entries, and set the histogram options for
  // drawing
  LOG(INFO) << "build 2D TH1 grids";
  std::vector<std::vector<std::vector<TH1D *>>> auau_grid_hard(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> auau_grid_match(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_hard(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_match(
      refcent_string.size());

  std::vector<std::vector<std::vector<TH1D *>>> auau_grid_hard_full(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> auau_grid_match_full(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_hard_full(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_match_full(
      refcent_string.size());

  std::vector<std::vector<std::vector<TH1D *>>> auau_oa_grid(
      refcent_string.size());
  std::vector<std::vector<std::vector<TGraphErrors *>>> error_grid_hard(
      refcent_string.size());
  std::vector<std::vector<std::vector<TGraphErrors *>>> error_grid_match(
      refcent_string.size());
  std::vector<std::vector<std::vector<TGraphErrors *>>> error_grid_hard_full(
      refcent_string.size());
  std::vector<std::vector<std::vector<TGraphErrors *>>> error_grid_match_full(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> err_frac_grid_hard(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> err_frac_grid_match(
      refcent_string.size());

  // the systematic variations
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_hard_tow_p(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_hard_tow_m(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_match_tow_p(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_match_tow_m(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_hard_trk_p(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_hard_trk_m(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_match_trk_p(
      refcent_string.size());
  std::vector<std::vector<std::vector<TH1D *>>> pp_grid_match_trk_m(
      refcent_string.size());

  for (int cent = 0; cent < refcent_string.size(); ++cent) {
    auau_grid_hard[cent].resize(radii.size());
    auau_grid_match[cent].resize(radii.size());
    pp_grid_hard[cent].resize(radii.size());
    pp_grid_match[cent].resize(radii.size());
    auau_grid_hard_full[cent].resize(radii.size());
    auau_grid_match_full[cent].resize(radii.size());
    pp_grid_hard_full[cent].resize(radii.size());
    pp_grid_match_full[cent].resize(radii.size());
    auau_oa_grid[cent].resize(radii.size());
    error_grid_hard[cent].resize(radii.size());
    error_grid_match[cent].resize(radii.size());
    error_grid_hard_full[cent].resize(radii.size());
    error_grid_match_full[cent].resize(radii.size());
    err_frac_grid_hard[cent].resize(radii.size());
    err_frac_grid_match[cent].resize(radii.size());
    pp_grid_hard_tow_p[cent].resize(radii.size());
    pp_grid_hard_tow_m[cent].resize(radii.size());
    pp_grid_match_tow_p[cent].resize(radii.size());
    pp_grid_match_tow_m[cent].resize(radii.size());
    pp_grid_hard_trk_p[cent].resize(radii.size());
    pp_grid_hard_trk_m[cent].resize(radii.size());
    pp_grid_match_trk_p[cent].resize(radii.size());
    pp_grid_match_trk_m[cent].resize(radii.size());
    for (int rad = 0; rad < radii.size(); ++rad) {
      auau_grid_hard[cent][rad].resize(constpt.size());
      auau_grid_match[cent][rad].resize(constpt.size());
      pp_grid_hard[cent][rad].resize(constpt.size());
      pp_grid_match[cent][rad].resize(constpt.size());
      auau_grid_hard_full[cent][rad].resize(constpt.size());
      auau_grid_match_full[cent][rad].resize(constpt.size());
      pp_grid_hard_full[cent][rad].resize(constpt.size());
      pp_grid_match_full[cent][rad].resize(constpt.size());
      auau_oa_grid[cent][rad].resize(constpt.size());
      error_grid_hard[cent][rad].resize(constpt.size());
      error_grid_match[cent][rad].resize(constpt.size());
      error_grid_hard_full[cent][rad].resize(constpt.size());
      error_grid_match_full[cent][rad].resize(constpt.size());
      err_frac_grid_hard[cent][rad].resize(constpt.size());
      err_frac_grid_match[cent][rad].resize(constpt.size());
      pp_grid_hard_tow_p[cent][rad].resize(constpt.size());
      pp_grid_hard_tow_m[cent][rad].resize(constpt.size());
      pp_grid_match_tow_p[cent][rad].resize(constpt.size());
      pp_grid_match_tow_m[cent][rad].resize(constpt.size());
      pp_grid_hard_trk_p[cent][rad].resize(constpt.size());
      pp_grid_hard_trk_m[cent][rad].resize(constpt.size());
      pp_grid_match_trk_p[cent][rad].resize(constpt.size());
      pp_grid_match_trk_m[cent][rad].resize(constpt.size());
      for (int pt = 0; pt < constpt.size(); ++pt) {
        string key = grid_keys[rad][pt];
        auau_grid_hard[cent][rad][pt] = hard_aj_cent[AUAU][key][cent];
        auau_grid_match[cent][rad][pt] = match_aj_cent[AUAU][key][cent];
        pp_grid_hard[cent][rad][pt] = hard_aj_cent[PP][key][cent];
        pp_grid_match[cent][rad][pt] = match_aj_cent[PP][key][cent];
        auau_grid_hard_full[cent][rad][pt] = hard_aj_full_cent[AUAU][key][cent];
        auau_grid_match_full[cent][rad][pt] =
            match_aj_full_cent[AUAU][key][cent];
        pp_grid_hard_full[cent][rad][pt] = hard_aj_full_cent[PP][key][cent];
        pp_grid_match_full[cent][rad][pt] = match_aj_full_cent[PP][key][cent];
        if (off_axis_aj_cent[AUAU].count(key)) {
          auau_oa_grid[cent][rad][pt] = off_axis_aj_cent[AUAU][key][cent];
        }
        error_grid_hard[cent][rad][pt] = systematic_errors_hard[key][cent];
        error_grid_match[cent][rad][pt] = systematic_errors_match[key][cent];
        error_grid_hard_full[cent][rad][pt] =
            systematic_errors_hard_full[key][cent];
        error_grid_match_full[cent][rad][pt] =
            systematic_errors_match_full[key][cent];
        err_frac_grid_hard[cent][rad][pt] = hard_aj_cent_err_frac[key][cent];
        err_frac_grid_match[cent][rad][pt] = match_aj_cent_err_frac[key][cent];
        pp_grid_hard_tow_p[cent][rad][pt] = hard_aj_cent[TOWP][key][cent];
        pp_grid_hard_tow_m[cent][rad][pt] = hard_aj_cent[TOWM][key][cent];
        pp_grid_match_tow_p[cent][rad][pt] = match_aj_cent[TOWP][key][cent];
        pp_grid_match_tow_m[cent][rad][pt] = match_aj_cent[TOWM][key][cent];
        pp_grid_hard_trk_p[cent][rad][pt] = hard_aj_cent[TRACKP][key][cent];
        pp_grid_hard_trk_m[cent][rad][pt] = hard_aj_cent[TRACKM][key][cent];
        pp_grid_match_trk_p[cent][rad][pt] = match_aj_cent[TRACKP][key][cent];
        pp_grid_match_trk_m[cent][rad][pt] = match_aj_cent[TRACKM][key][cent];
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // format each histogram for visually - color choices, etc. Each histogram
  // type (hard auau aj, match pp aj, etc) will be formatted the same way
  LOG(INFO) << "formatting histograms";
  for (auto &container1 : auau_grid_hard)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.5);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : auau_grid_match)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.5);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : pp_grid_hard)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kRed);
        h->SetMarkerColor(kRed);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.5);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : pp_grid_match)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kRed);
        h->SetMarkerColor(kRed);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.5);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : auau_grid_hard_full)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.5);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : auau_grid_match_full)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.5);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : pp_grid_hard_full)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kRed);
        h->SetMarkerColor(kRed);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.5);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : pp_grid_match_full)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kRed);
        h->SetMarkerColor(kRed);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.5);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : auau_oa_grid)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kAzure);
        h->SetMarkerColor(kAzure);
        h->SetFillColorAlpha(kAzure, 0.65);
        h->SetMarkerSize(0);
        h->SetLineWidth(0);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : error_grid_hard)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetFillStyle(1001);
        h->SetLineWidth(0);
        h->SetMarkerSize(0);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : error_grid_match)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetFillStyle(1001);
        h->SetLineWidth(0);
        h->SetMarkerSize(0);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : error_grid_hard_full)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetFillStyle(1001);
        h->SetLineWidth(0);
        h->SetMarkerSize(0);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : error_grid_match_full)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetFillStyle(1001);
        h->SetLineWidth(0);
        h->SetMarkerSize(0);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : err_frac_grid_hard)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.1);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : err_frac_grid_match)
    for (auto &container2 : container1)
      for (auto &h : container2) {
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
        h->SetLineWidth(1);
        h->SetMarkerSize(0.1);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->GetXaxis()->SetNdivisions(305);
        h->GetYaxis()->SetNdivisions(305);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetLabelSize(10);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetLabelSize(10);
        h->GetYaxis()->SetRangeUser(0.0001, 0.249);
      }

  for (auto &container1 : {pp_grid_hard_tow_p, pp_grid_hard_trk_p,
                           pp_grid_match_tow_p, pp_grid_match_trk_p})
    for (auto &container2 : container1)
      for (auto &container3 : container2)
        for (auto &h : container3) {
          h->SetLineColor(kAzure);
          h->SetMarkerColor(kAzure);
          h->SetLineWidth(1);
          h->SetMarkerSize(0.1);
          h->GetXaxis()->SetTitle("");
          h->GetYaxis()->SetTitle("");
          h->GetXaxis()->SetNdivisions(305);
          h->GetYaxis()->SetNdivisions(305);
          h->GetXaxis()->SetLabelFont(43);
          h->GetXaxis()->SetLabelSize(10);
          h->GetYaxis()->SetLabelFont(43);
          h->GetYaxis()->SetLabelSize(10);
          h->GetYaxis()->SetRangeUser(0.0001, 0.249);
        }

  for (auto &container1 : {pp_grid_hard_tow_m, pp_grid_hard_trk_m,
                           pp_grid_match_tow_m, pp_grid_match_trk_m})
    for (auto &container2 : container1)
      for (auto &container3 : container2)
        for (auto &h : container3) {
          h->SetLineColor(kOrange);
          h->SetMarkerColor(kOrange);
          h->SetLineWidth(1);
          h->SetMarkerSize(0.1);
          h->GetXaxis()->SetTitle("");
          h->GetYaxis()->SetTitle("");
          h->GetXaxis()->SetNdivisions(305);
          h->GetYaxis()->SetNdivisions(305);
          h->GetXaxis()->SetLabelFont(43);
          h->GetXaxis()->SetLabelSize(10);
          h->GetYaxis()->SetLabelFont(43);
          h->GetYaxis()->SetLabelSize(10);
          h->GetYaxis()->SetRangeUser(0.0001, 0.249);
        }

  // ------------------------------------------------------------------------------------
  // define the sets of text that are to be printed in each bin on the canvas
  // - the outer vector is for centrality, the inner is the list of text
  // objects for each grid
  LOG(INFO) << "creating text for grid printouts";
  std::vector<std::vector<dijetcore::GridTextObject>> hard_aj_text(
      refcent_string.size());
  std::vector<std::vector<dijetcore::GridTextObject>> match_aj_text(
      refcent_string.size());
  std::vector<std::vector<dijetcore::GridTextObject>> oa_aj_text(
      refcent_string.size());
  std::vector<std::vector<dijetcore::GridTextObject>> hard_err_text(
      refcent_string.size());
  std::vector<std::vector<dijetcore::GridTextObject>> match_err_text(
      refcent_string.size());

  // add legends - legends will be in the top left corner (0xlast)
  // the legend position will have to be adjusted when the grid changes,
  // unfortunately. Haven't come up with a better way of changing it
  // automatically
  for (int cent = 0; cent < refcent_string.size(); ++cent) {
    unsigned leg_x_index = 0;
    unsigned leg_y_index = constpt.size() - 1;
    double x_low = 0.65;
    double x_high = 0.95;
    double y_low = 0.4;
    double y_high = 0.6;

    hard_aj_text[cent].push_back(dijetcore::GridTextObject());
    hard_aj_text[cent].back().location_x = leg_x_index;
    hard_aj_text[cent].back().location_y = leg_y_index;
    hard_aj_text[cent].back().legend = std::move(
        dijetcore::make_unique<TLegend>(x_low, y_low, x_high, y_high));
    hard_aj_text[cent].back().legend->AddEntry(
        auau_grid_hard[cent][leg_x_index][leg_y_index], "Au+Au HT");
    hard_aj_text[cent].back().legend->AddEntry(
        pp_grid_hard[cent][leg_x_index][leg_y_index], "p+p #oplus Au+Au MB");

    match_aj_text[cent].push_back(dijetcore::GridTextObject());
    match_aj_text[cent].back().location_x = leg_x_index;
    match_aj_text[cent].back().location_y = leg_y_index;
    match_aj_text[cent].back().legend = std::move(
        dijetcore::make_unique<TLegend>(x_low, y_low, x_high, y_high));
    match_aj_text[cent].back().legend->AddEntry(
        auau_grid_match[cent][leg_x_index][leg_y_index], "Au+Au HT");
    match_aj_text[cent].back().legend->AddEntry(
        pp_grid_match[cent][leg_x_index][leg_y_index], "p+p #oplus Au+Au MB");

    oa_aj_text[cent].push_back(dijetcore::GridTextObject());
    oa_aj_text[cent].back().location_x = leg_x_index;
    oa_aj_text[cent].back().location_y = leg_y_index;
    oa_aj_text[cent].back().legend = std::move(
        dijetcore::make_unique<TLegend>(x_low, y_low, x_high, y_high));
    oa_aj_text[cent].back().legend->AddEntry(
        auau_grid_match[cent][leg_x_index][leg_y_index], "Au+Au HT");
    oa_aj_text[cent].back().legend->AddEntry(
        auau_oa_grid[cent][leg_x_index][leg_y_index],
        "Au+Au HT #oplus Au+Au MB");

    hard_err_text[cent].push_back(dijetcore::GridTextObject());
    hard_err_text[cent].back().location_x = leg_x_index;
    hard_err_text[cent].back().location_y = leg_y_index;
    hard_err_text[cent].back().legend = std::move(
        dijetcore::make_unique<TLegend>(x_low, y_low, x_high, y_high));
    hard_err_text[cent].back().legend->AddEntry(
        err_frac_grid_hard[cent][leg_x_index][leg_y_index],
        "p+p relative sys. err.");

    match_err_text[cent].push_back(dijetcore::GridTextObject());
    match_err_text[cent].back().location_x = leg_x_index;
    match_err_text[cent].back().location_y = leg_y_index;
    match_err_text[cent].back().legend = std::move(
        dijetcore::make_unique<TLegend>(x_low, y_low, x_high, y_high));
    match_err_text[cent].back().legend->AddEntry(
        err_frac_grid_match[cent][leg_x_index][leg_y_index],
        "p+p relative sys. err.");
  }

  // add a text box with centrality, jet pT, etc. This is almost same for all
  // grids
  for (int cent = 0; cent < refcent_string.size(); ++cent) {
    unsigned text_x_index = 1;
    unsigned text_y_index = constpt.size() - 1;
    double x_low = 0.68;
    double x_high = 0.88;
    double y_low = 0.3;
    double y_high = 0.6;

    std::string cent_string = "Au+Au, " + refcent_string[cent];
    std::string pt_lead_string = dijetcore::MakeString(
        "p_{T}^{hard,lead} > ", grid_key_params[0][0].lead_init_pt, " GeV/c");
    std::string pt_sub_string = dijetcore::MakeString(
        "p_{T}^{hard,sublead} > ", grid_key_params[0][0].sub_init_pt, " GeV/c");
    std::string match_pt_string = "p_{T}^{match const} > 0.2 GeV/c";

    // hard di-jets get this one, which is lacking match jet pTconst
    TPaveText gen_hard_info(x_low, y_low, x_high, y_high, "NB NDC");
    gen_hard_info.SetFillStyle(0);
    gen_hard_info.SetBorderSize(0);
    gen_hard_info.AddText(cent_string.c_str())->SetTextSize(0.038);
    gen_hard_info.AddText(pt_lead_string.c_str())->SetTextSize(0.038);
    gen_hard_info.AddText(pt_sub_string.c_str())->SetTextSize(0.038);

    TPaveText gen_match_info(x_low, y_low, x_high, y_high, "NB NDC");
    gen_match_info.SetFillStyle(0);
    gen_match_info.SetBorderSize(0);
    gen_match_info.AddText(cent_string.c_str())->SetTextSize(0.038);
    gen_match_info.AddText(pt_lead_string.c_str())->SetTextSize(0.038);
    gen_match_info.AddText(pt_sub_string.c_str())->SetTextSize(0.038);
    gen_match_info.AddText(match_pt_string.c_str())->SetTextSize(0.038);

    hard_aj_text[cent].push_back(dijetcore::GridTextObject());
    hard_aj_text[cent].back().location_x = text_x_index;
    hard_aj_text[cent].back().location_y = text_y_index;
    hard_aj_text[cent].back().text =
        std::move(dijetcore::make_unique<TPaveText>(gen_hard_info));

    match_aj_text[cent].push_back(dijetcore::GridTextObject());
    match_aj_text[cent].back().location_x = text_x_index;
    match_aj_text[cent].back().location_y = text_y_index;
    match_aj_text[cent].back().text =
        std::move(dijetcore::make_unique<TPaveText>(gen_match_info));

    oa_aj_text[cent].push_back(dijetcore::GridTextObject());
    oa_aj_text[cent].back().location_x = text_x_index;
    oa_aj_text[cent].back().location_y = text_y_index;
    oa_aj_text[cent].back().text =
        std::move(dijetcore::make_unique<TPaveText>(gen_match_info));

    hard_err_text[cent].push_back(dijetcore::GridTextObject());
    hard_err_text[cent].back().location_x = text_x_index;
    hard_err_text[cent].back().location_y = text_y_index;
    hard_err_text[cent].back().text =
        std::move(dijetcore::make_unique<TPaveText>(gen_hard_info));

    match_err_text[cent].push_back(dijetcore::GridTextObject());
    match_err_text[cent].back().location_x = text_x_index;
    match_err_text[cent].back().location_y = text_y_index;
    match_err_text[cent].back().text =
        std::move(dijetcore::make_unique<TPaveText>(gen_match_info));
  }

  // ------------------------------------------------------------------------------------
  // define axis labels for grids

  TText *tmp_text = nullptr;
  TPaveText y_axis_text(0.02, 0.2, 0.08, 0.8, "NB NDC");
  y_axis_text.SetFillStyle(0);
  y_axis_text.SetBorderSize(0);
  tmp_text = y_axis_text.AddText("event fraction");
  tmp_text->SetTextSize(0.04);
  tmp_text->SetTextAngle(90);
  TPaveText y_axis_text_err_frac(0.02, 0.2, 0.08, 0.8, "NB NDC");
  y_axis_text_err_frac.SetFillStyle(0);
  y_axis_text_err_frac.SetBorderSize(0);
  tmp_text = y_axis_text_err_frac.AddText("percent error");
  tmp_text->SetTextSize(0.04);
  tmp_text->SetTextAngle(90);
  TPaveText x_axis_text(0.4, 0.02, 0.6, 0.08, "NB NDC");
  x_axis_text.SetFillStyle(0);
  x_axis_text.SetBorderSize(0);
  tmp_text = x_axis_text.AddText("|A_{J}|");
  tmp_text->SetTextSize(0.04);
  TPaveText y_axis_text_right(0.94, 0.2, 0.99, 0.8, "NB NDC");
  y_axis_text_right.SetFillStyle(0);
  y_axis_text_right.SetBorderSize(0);
  tmp_text = y_axis_text_right.AddText("p_{T}^{hard const} [GeV/c]");
  tmp_text->SetTextSize(0.04);
  tmp_text->SetTextAngle(270);
  TPaveText x_axis_text_top(0.4, 0.94, 0.6, 0.99, "NB NDC");
  x_axis_text_top.SetFillStyle(0);
  x_axis_text_top.SetBorderSize(0);
  tmp_text = x_axis_text_top.AddText("resolution parameter");
  tmp_text->SetTextSize(0.04);

  std::vector<TPaveText> axis_text_aj;
  std::vector<TPaveText> axis_text_err_frac;

  axis_text_aj.push_back(y_axis_text);
  axis_text_aj.push_back(x_axis_text);
  axis_text_aj.push_back(y_axis_text_right);
  axis_text_aj.push_back(x_axis_text_top);

  axis_text_err_frac.push_back(y_axis_text_err_frac);
  axis_text_err_frac.push_back(x_axis_text);
  axis_text_err_frac.push_back(y_axis_text_right);
  axis_text_err_frac.push_back(x_axis_text_top);

  // and add the actual grid labels
  std::vector<TPaveText> grid_labels = dijetcore::GetGridLabelsForAj(
      radii.size(), constpt.size(), radii.front(), radii.back(),
      constpt.front(), constpt.back(), grid_x_low_margin, grid_x_high_margin,
      grid_y_low_margin, grid_y_high_margin);

  axis_text_aj.insert(axis_text_aj.end(), grid_labels.begin(),
                      grid_labels.end());
  axis_text_err_frac.insert(axis_text_err_frac.end(), grid_labels.begin(),
                            grid_labels.end());

  // ------------------------------------------------------------------------------------
  // print grids
  LOG(INFO) << "printing grids";
  for (int cent = 0; cent < refcent_string.size(); ++cent) {
    // create output
    boost::filesystem::path grid_dir(FLAGS_outputDir);
    grid_dir /= "grid_print";
    grid_dir /= refcent_string[cent];
    boost::filesystem::create_directories(grid_dir);

    // create dummy grid for when we only have 2 or 1 layers
    std::vector<std::vector<TH1D *>> dummy_grid;

    // create invisible full width/height pad for drawing Axis text, etc
    TPad *invis = new TPad("invis_pad", "", 0, 0, 1, 1);
    invis->SetFillStyle(4000);

    // print hard aj
    dijetcore::GridPrintOptions hard_opts;
    hard_opts.layer_2_print = "e2";
    boost::filesystem::path hard_aj_name = grid_dir;
    hard_aj_name /= "hard_aj.pdf";
    dijetcore::PrintGrid(canvas_hard[cent], hard_pads[cent], pp_grid_hard[cent],
                         error_grid_hard[cent], auau_grid_hard[cent], hard_opts,
                         hard_aj_text[cent], invis, axis_text_aj,
                         hard_aj_name.string());

    // print match aj
    dijetcore::GridPrintOptions match_opts;
    match_opts.layer_2_print = "e2";
    boost::filesystem::path match_aj_name = grid_dir;
    match_aj_name /= "match_aj.pdf";
    dijetcore::PrintGrid(canvas_match[cent], match_pads[cent],
                         pp_grid_match[cent], error_grid_match[cent],
                         auau_grid_match[cent], match_opts, match_aj_text[cent],
                         invis, axis_text_aj, match_aj_name.string());

    // print hard aj full
    dijetcore::GridPrintOptions hard_opts_full;
    hard_opts_full.layer_2_print = "e2";
    boost::filesystem::path hard_aj_full_name = grid_dir;
    hard_aj_full_name /= "hard_aj_full.pdf";
    dijetcore::PrintGrid(
        canvas_hard_full[cent], hard_pads_full[cent], pp_grid_hard_full[cent],
        error_grid_hard_full[cent], auau_grid_hard_full[cent], hard_opts_full,
        hard_aj_text[cent], invis, axis_text_aj, hard_aj_full_name.string());

    // print match aj
    dijetcore::GridPrintOptions match_opts_full;
    match_opts_full.layer_2_print = "e2";
    boost::filesystem::path match_aj_full_name = grid_dir;
    match_aj_full_name /= "match_aj_full.pdf";
    dijetcore::PrintGrid(canvas_match_full[cent], match_pads_full[cent],
                         pp_grid_match_full[cent], error_grid_match_full[cent],
                         auau_grid_match_full[cent], match_opts_full,
                         match_aj_text[cent], invis, axis_text_aj,
                         match_aj_full_name.string());

    // print off axis aj
    dijetcore::GridPrintOptions oa_opts;
    oa_opts.layer_3_active = false;
    oa_opts.layer_1_print = "e3";
    boost::filesystem::path oa_aj_name = grid_dir;
    oa_aj_name /= "oa_aj.pdf";
    dijetcore::PrintGrid(canvas_oa[cent], oa_pads[cent], auau_oa_grid[cent],
                         auau_grid_match[cent], dummy_grid, oa_opts,
                         oa_aj_text[cent], invis, axis_text_aj,
                         oa_aj_name.string());
  }

  return 0;
}
