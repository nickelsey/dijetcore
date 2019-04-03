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
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <string>
#include "dijetcore/lib/map.h"

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

enum class DATATYPE {
  AUAU = 0,
  PP = 1,
  TOWP = 2,
  TOWM = 3,
  TRACKP = 4,
  TRACKM = 5,
  SIZE = 6
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
void GetTreesFromFile(const TFile& file,
                      std::unordered_map<string, TTree*>& map) {
  TKey* key;
  TIter next(file.GetListOfKeys());

  while ((key = (TKey*)next())) {
    // check if its a TTree
    TClass* cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TTree")) continue;

    // check if its produced by the DijetWorker
    // (in a rather handwavy fashion)
    string tmp(key->GetName());
    if (tmp.find("LEAD_INIT") == string::npos ||
        tmp.find("SUB_INIT") == string::npos)
      continue;

    map.insert({tmp, (TTree*)key->ReadObj()});
  }
}

std::vector<TH1D*> SplitByCentrality(TH2D* h) {
  std::vector<TH1D*> ret;
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    string name = string(h->GetName()) + std::to_string(i);
    TH1D* tmp = h->ProjectionY(name.c_str(), i, i);
    ret.push_back(tmp);
  }
  return ret;
}

std::vector<TH2D*> SplitByCentrality3D(TH3D* h) {
  std::vector<TH2D*> ret;
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    string name = string(h->GetName()) + std::to_string(i);
    h->GetXaxis()->SetRange(i, i);
    TH2D* tmp = (TH2D*)h->Project3D("zy");
    tmp->SetName(name.c_str());
    ret.push_back(tmp);
  }
  return ret;
}

std::vector<TH1D*> SplitByCentralityNormalized(TH2D* h, int bins = 9) {
  std::vector<TH1D*> ret;
  for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
    string name = string(h->GetName()) + std::to_string(i);
    TH1D* tmp = h->ProjectionY(name.c_str(), i, i);
    if (tmp->Integral() > 0) tmp->Scale(1.0 / tmp->Integral());
    ret.push_back(tmp);
  }
  return ret;
}

TGraphErrors* GetSystematic(TH1D* nom, TH1D* var1_a, TH1D* var1_b, TH1D* var2_a,
                            TH1D* var2_b) {
  int nBins = nom->GetNbinsX();
  double x_[nBins];
  double y_[nBins];
  double x_err_[nBins];
  double y_err_[nBins];

  for (int i = 0; i < nBins; ++i) {
    x_[i] = nom->GetBinCenter(i + 1);
    y_[i] = nom->GetBinContent(i + 1);
    x_err_[i] = nom->GetXaxis()->GetBinWidth(1) / 2.0;
    double diff_var_1_a =
        fabs(nom->GetBinContent(i + 1) - var1_a->GetBinContent(i + 1));
    double diff_var_1_b =
        fabs(nom->GetBinContent(i + 1) - var1_b->GetBinContent(i + 1));
    double diff_var_2_a =
        fabs(nom->GetBinContent(i + 1) - var2_a->GetBinContent(i + 1));
    double diff_var_2_b =
        fabs(nom->GetBinContent(i + 1) - var2_b->GetBinContent(i + 1));
    double max_var_1 =
        (diff_var_1_a > diff_var_1_b ? diff_var_1_a : diff_var_1_b);
    double max_var_2 =
        (diff_var_2_a > diff_var_2_b ? diff_var_2_a : diff_var_2_b);
    y_err_[i] = sqrt(max_var_1 * max_var_1 + max_var_2 * max_var_2);
  }
  TGraphErrors* ret = new TGraphErrors(nBins, x_, y_, x_err_, y_err_);
  return ret;
}

template <class T>
TH1D* GetSystematicFractional(T* nom, TGraphErrors* errors) {
  if (nom->GetXaxis()->GetNbins() != errors->GetN()) {
    LOG(ERROR) << "bin mismatch: exiting";
    return nullptr;
  }

  string name = nom->GetName() + "errfrac";
  TH1D* ret = new TH1D(name.c_str(), "", nom->GetXaxis()->GetXbins(),
                       nom->GetXaxis()->GetXmin(), nom->GetXaxis()->GetXmax());
  for (int i = 1; i < nom->GetXaxis()->GetXbins(); ++i) {
    double bin_content = nom->GetBinContent(i);
    double bin_error = errors->GetErrorY(i - 1);
    ret->SetBinContent(i, bin_error / bin_content);
  }
  return ret;
}

template <typename H>
void AjPrintout(H* h1, H* h2, TGraphErrors* sys, double y_min, double y_max,
                double x_min, double x_max, TPaveText& text,
                std::string h1_title, std::string h2_title,
                dijetcore::histogramOpts hopts, dijetcore::canvasOpts copts,
                std::string output_loc, std::string output_name,
                std::string canvas_title, std::string x_axis_label,
                std::string y_axis_label, std::string legend_title = "") {
  // we assume the output location exists, so create
  // the final output string that will be used for pdf creation
  std::string canvas_name = output_loc + "/" + output_name + ".pdf";

  // axis labels, and title
  h1->SetTitle(canvas_title.c_str());
  h1->GetXaxis()->SetTitle(x_axis_label.c_str());
  h1->GetYaxis()->SetTitle(y_axis_label.c_str());
  h1->GetYaxis()->SetRangeUser(y_min, y_max);
  h1->GetXaxis()->SetRangeUser(x_min, x_max);

  // set draw options
  hopts.SetHistogram(h1);
  hopts.SetHistogram(h2);

  if (sys == nullptr) {
    h2->SetLineColor(kBlue);
    h2->SetMarkerColor(kBlue);
  }

  TCanvas c;
  copts.SetMargins(&c);
  copts.SetLogScale(&c);

  h1->Draw("9");
  h2->Draw("9SAME");

  if (sys != nullptr) {
    sys->SetFillColorAlpha(h2->GetLineColor(), 0.45);
    sys->SetFillStyle(1001);
    sys->SetLineWidth(0);
    sys->SetMarkerSize(0);
    sys->Draw("9e2SAME");
  }

  TLegend* leg = copts.Legend();
  if (leg != nullptr) {
    leg->SetHeader(legend_title.c_str());
    leg->SetTextFont(62);
    leg->AddEntry(h1, h1_title.c_str(), "p")->SetTextSize(.035);
    leg->AddEntry(h2, h2_title.c_str(), "p")->SetTextSize(.035);
    leg->Draw();
  }

  // and STAR preliminary message
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.065);
  latex.SetTextColor(kRed + 3);
  // latex.DrawLatex(0.19, 0.8, "STAR Preliminary");

  text.Draw();

  c.SaveAs(canvas_name.c_str());
}

TLegend* GetLegend(int x_bin, int y_bin, int x_max, int y_max, double radius,
                   double pt, double scale = 1.0) {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05 * scale);
  latex.SetTextColor(kBlack);
  TLegend* leg;
  std::stringstream in;
  in << std::setprecision(2) << radius;
  std::string radius_label = "Radius = " + in.str();
  in.clear();
  in << std::setprecision(2) << pt;
  std::string pt_label = "p_{T}^{const} > " + in.str() + " GeV/c";
  if (x_bin == 0 && y_bin == 0) {
    // latex.DrawLatex(0.8, 0.8, radius_label.c_str());
    // latex.DrawLatex(0.8, 0.73, pt_label.c_str());
    leg = new TLegend(0.6, 0.8, 0.95, 0.95);
  } else if (x_bin == x_max && y_bin == y_max) {
    // latex.DrawLatex(0.3, 0.5, radius_label.c_str());
    // latex.DrawLatex(0.3, 0.42, pt_label.c_str());
    leg = new TLegend(0.28, 0.6, 0.5, 0.76);
  } else if (x_bin == x_max && y_bin == 0) {
    // latex.DrawLatex(0.3, 0.8, radius_label.c_str());
    // latex.DrawLatex(0.3, 0.73, pt_label.c_str());
    leg = new TLegend(0.28, 0.8, 0.5, 0.95);
  } else if (x_bin == 0 && y_bin == y_max) {
    // latex.DrawLatex(0.8, 0.5, radius_label.c_str());
    // latex.DrawLatex(0.8, 0.42, pt_label.c_str());
    leg = new TLegend(0.6, 0.6, 0.95, 0.76);
  } else if (x_bin == 0) {
    // latex.DrawLatex(0.8, 0.6, radius_label.c_str());
    // latex.DrawLatex(0.8, 0.5, pt_label.c_str());
    leg = new TLegend(0.6, 0.7, 0.95, 0.95);
  } else if (x_bin == x_max) {
    // latex.DrawLatex(0.3, 0.6, radius_label.c_str());
    // latex.DrawLatex(0.3, 0.5, pt_label.c_str());
    leg = new TLegend(0.28, 0.7, 0.5, 0.95);
  } else if (y_bin == 0) {
    // latex.DrawLatex(0.7, 0.8, radius_label.c_str());
    // latex.DrawLatex(0.7, 0.73, pt_label.c_str());
    leg = new TLegend(0.48, 0.8, 0.95, 0.95);
  } else if (y_bin == y_max) {
    // latex.DrawLatex(0.7, 0.5, radius_label.c_str());
    // latex.DrawLatex(0.7, 0.42, pt_label.c_str());
    leg = new TLegend(0.48, 0.6, 0.95, 0.76);
  } else {
    // latex.DrawLatex(0.6, 0.5, radius_label.c_str());
    // latex.DrawLatex(0.6, 0.45, pt_label.c_str());
    leg = new TLegend(0.48, 0.7, 0.95, 0.95);
  }
  leg->SetTextFont(62);
  return leg;
}

TGraphErrors* MakeGraph(TH1D* h) {
  int bins = h->GetNbinsX();
  double x[bins];
  double x_err[bins];
  double y[bins];
  double y_err[bins];

  for (int i = 0; i < bins; ++i) {
    x[i] = h->GetBinCenter(i + 1);
    y[i] = h->GetBinContent(i + 1);
    x_err[i] = 0.0;
    y_err[i] = h->GetBinError(i + 1);
  }

  TGraphErrors* err = new TGraphErrors(bins, x, y, x_err, y_err);
  return err;
}

int main(int argc, char* argv[]) {
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

  // check to make sure we have valid inputs
  std::vector<string> inputs{FLAGS_auau, FLAGS_ppDir + "/tow_0_track_0.root"};
  for (auto& file : inputs) {
    if (!boost::filesystem::exists(file)) {
      std::cout << "input file " << file;
      std::cout << "doesn't exist: exiting" << std::endl;
      return 1;
    }
  }

  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty()) FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);

  // read in the file
  TFile auau_file(FLAGS_auau.c_str(), "READ");
  TFile pp_file((FLAGS_ppDir + "/tow_0_track_0.root").c_str(), "READ");
  TFile tow_p_file((FLAGS_ppDir + "/tow_1_track_0.root").c_str(), "READ");
  TFile tow_m_file((FLAGS_ppDir + "/tow_-1_track_0.root").c_str(), "READ");
  TFile track_p_file((FLAGS_ppDir + "/tow_0_track_1.root").c_str(), "READ");
  TFile track_m_file((FLAGS_ppDir + "/tow_0_track_-1.root").c_str(), "READ");

  // define centralities
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

  // and the radii & constituent pT we want to use
  std::vector<double> radii =
      dijetcore::ParseArgStringToVec<double>(FLAGS_radii);
  std::sort(radii.begin(), radii.end());
  std::vector<string> radii_string;
  for (auto& val : radii) {
    std::stringstream stream;
    stream << val;
    radii_string.push_back(stream.str());
  }
  std::vector<double> constpt =
      dijetcore::ParseArgStringToVec<double>(FLAGS_constPt);
  std::sort(constpt.begin(), constpt.end());
  std::vector<string> constpt_string;
  for (auto& val : constpt) {
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
  std::vector<std::unordered_map<string, TTree*>> trees(TYPEMAP.size());

  GetTreesFromFile(auau_file, trees[TYPEMAP[DATATYPE::AUAU]]);
  GetTreesFromFile(pp_file, trees[TYPEMAP[DATATYPE::PP]]);
  GetTreesFromFile(tow_p_file, trees[TYPEMAP[DATATYPE::TOWP]]);
  GetTreesFromFile(tow_m_file, trees[TYPEMAP[DATATYPE::TOWM]]);
  GetTreesFromFile(track_p_file, trees[TYPEMAP[DATATYPE::TRACKP]]);
  GetTreesFromFile(track_m_file, trees[TYPEMAP[DATATYPE::TRACKM]]);

  // match keys
  int auau_index = static_cast<int>(DATATYPE::AUAU);
  int pp_index = static_cast<int>(DATATYPE::PP);
  for (auto entry : trees[auau_index]) {
    if (trees[pp_index].find(entry.first) != trees[pp_index].end()) {
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
      for (auto& key_string : keys) {
        dijetcore::DijetKey& key = parsed_keys[key_string];

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

  const int num_datatypes = TYPEMAP.size();

  std::vector<std::unordered_map<string, TH2D*>> hard_lead_eta(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_lead_phi(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_lead_pt(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_lead_const(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_lead_rho(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_lead_sig(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_lead_eta(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_lead_phi(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_lead_pt(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_lead_const(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_lead_rho(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_lead_sig(num_datatypes);

  std::vector<std::unordered_map<string, TH2D*>> hard_sub_eta(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_sub_phi(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_sub_pt(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_sub_const(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_sub_rho(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_sub_sig(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_sub_eta(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_sub_phi(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_sub_pt(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_sub_const(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_sub_rho(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_sub_sig(num_datatypes);

  std::vector<std::unordered_map<string, TH2D*>> npart(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_dphi(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_dphi(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> lead_dr(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> sub_dr(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> lead_dpt(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> sub_dpt(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> lead_dpt_frac(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> sub_dpt_frac(num_datatypes);

  std::vector<std::unordered_map<string, TH2D*>> hard_aj(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_aj(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_aj_test(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_aj_test(num_datatypes);

  std::vector<std::unordered_map<string, TH2D*>> off_axis_aj(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> off_axis_aj_test(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> off_axis_rho(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> off_axis_sig(num_datatypes);

  std::vector<std::unordered_map<string, TH2D*>> hard_pp_only_aj(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_pp_only_aj(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_lead_pp_dpt(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_lead_pp_dpt(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_sub_pp_dpt(num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_sub_pp_dpt(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_lead_pp_dpt_frac(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_lead_pp_dpt_frac(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> hard_sub_pp_dpt_frac(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> match_sub_pp_dpt_frac(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> pp_lead_hard_match_fraction(
      num_datatypes);
  std::vector<std::unordered_map<string, TH2D*>> pp_sub_hard_match_fraction(
      num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_lead_eta_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_lead_phi_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_lead_pt_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_lead_const_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_lead_rho_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_lead_sig_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_lead_eta_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_lead_phi_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_lead_pt_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_lead_const_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_lead_rho_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_lead_sig_cent(num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_sub_eta_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_sub_phi_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_sub_pt_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_sub_const_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_sub_rho_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_sub_sig_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_sub_eta_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_sub_phi_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> match_sub_pt_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_sub_const_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_sub_rho_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_sub_sig_cent(num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D*>>> npart_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_dphi_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> match_dphi_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> lead_dr_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> sub_dr_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> lead_dpt_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> sub_dpt_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      lead_dpt_frac_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> sub_dpt_frac_cent(
      num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_aj_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> match_aj_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> hard_aj_test_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_aj_test_cent(num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D*>>> off_axis_aj_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      off_axis_aj_test_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> off_axis_rho_cent(
      num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>> off_axis_sig_cent(
      num_datatypes);

  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_pp_only_aj_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_pp_only_aj_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_lead_pp_dpt_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_lead_pp_dpt_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_sub_pp_dpt_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_sub_pp_dpt_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_lead_pp_dpt_frac_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_lead_pp_dpt_frac_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      hard_sub_pp_dpt_frac_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      match_sub_pp_dpt_frac_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      pp_lead_hard_match_fraction_cent(num_datatypes);
  std::vector<std::unordered_map<string, std::vector<TH1D*>>>
      pp_sub_hard_match_fraction_cent(num_datatypes);

  std::unordered_map<string, std::vector<TGraphErrors*>> systematic_errors_hard;
  std::unordered_map<string, std::vector<TGraphErrors*>>
      systematic_errors_match;

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

  // count which entry we're on
  int entry = -1;

  // loop over all keys
  for (auto& key : keys) {
    // increment entry counter
    entry++;
    LOG(INFO) << "loading: " << key;
    // make our output directory
    string key_loc = FLAGS_outputDir + "/" + key;
    boost::filesystem::path dir(key_loc.c_str());
    boost::filesystem::create_directories(dir);

    // now iterate over all our datasets
    for (auto type : TYPEMAP) {
      auto& type_enum = type.first;
      auto& data_index = type.second;

      // make our output directory
      string out_loc = key_loc + "/" + TYPENAME[type_enum];
      boost::filesystem::path dir(out_loc.c_str());
      boost::filesystem::create_directories(dir);

      // load our reader
      TTree* current_tree = trees[data_index][key];
      TTreeReader reader(current_tree);

      // create readervalues for auau first
      dijetcore::unique_ptr<TTreeReaderValue<int>> runid =
          dijetcore::make_unique<TTreeReaderValue<int>>(reader, "runid");
      dijetcore::unique_ptr<TTreeReaderValue<int>> eventid =
          dijetcore::make_unique<TTreeReaderValue<int>>(reader, "eventid");
      dijetcore::unique_ptr<TTreeReaderValue<double>> vz =
          dijetcore::make_unique<TTreeReaderValue<double>>(reader, "vz");
      dijetcore::unique_ptr<TTreeReaderValue<int>> refmult =
          dijetcore::make_unique<TTreeReaderValue<int>>(reader, "refmult");
      dijetcore::unique_ptr<TTreeReaderValue<int>> grefmult =
          dijetcore::make_unique<TTreeReaderValue<int>>(reader, "grefmult");
      dijetcore::unique_ptr<TTreeReaderValue<double>> refmultcorr =
          dijetcore::make_unique<TTreeReaderValue<double>>(reader,
                                                           "refmultcorr");
      dijetcore::unique_ptr<TTreeReaderValue<double>> grefmultcorr =
          dijetcore::make_unique<TTreeReaderValue<double>>(reader,
                                                           "grefmultcorr");
      dijetcore::unique_ptr<TTreeReaderValue<int>> cent =
          dijetcore::make_unique<TTreeReaderValue<int>>(reader, "cent");
      dijetcore::unique_ptr<TTreeReaderValue<double>> zdcrate =
          dijetcore::make_unique<TTreeReaderValue<double>>(reader, "zdcrate");
      dijetcore::unique_ptr<TTreeReaderValue<double>> rp =
          dijetcore::make_unique<TTreeReaderValue<double>>(reader, "rp");
      dijetcore::unique_ptr<TTreeReaderValue<int>> nglobal =
          dijetcore::make_unique<TTreeReaderValue<int>>(reader, "nglobal");
      dijetcore::unique_ptr<TTreeReaderValue<int>> nprt =
          dijetcore::make_unique<TTreeReaderValue<int>>(reader, "npart");
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jl =
          dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
                                                                   "jl");
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> js =
          dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
                                                                   "js");
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jlm =
          dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
                                                                   "jlm");
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jsm =
          dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
                                                                   "jsm");
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

      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jloa;
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jsoa;
      dijetcore::unique_ptr<TTreeReaderValue<int>> jloaconst;
      dijetcore::unique_ptr<TTreeReaderValue<double>> jloarho;
      dijetcore::unique_ptr<TTreeReaderValue<double>> jloasig;
      dijetcore::unique_ptr<TTreeReaderValue<int>> jsoaconst;
      dijetcore::unique_ptr<TTreeReaderValue<double>> jsoarho;
      dijetcore::unique_ptr<TTreeReaderValue<double>> jsoasig;

      // build all embedding branches (only used if pp_embedded is true)
      dijetcore::unique_ptr<TTreeReaderValue<int>> embed_eventid;
      dijetcore::unique_ptr<TTreeReaderValue<int>> embed_runid;
      dijetcore::unique_ptr<TTreeReaderValue<int>> embed_refmult;
      dijetcore::unique_ptr<TTreeReaderValue<int>> embed_grefmult;
      dijetcore::unique_ptr<TTreeReaderValue<int>> embed_nprt;
      dijetcore::unique_ptr<TTreeReaderValue<double>> embed_refmultcorr;
      dijetcore::unique_ptr<TTreeReaderValue<double>> embed_grefmultcorr;
      dijetcore::unique_ptr<TTreeReaderValue<int>> embed_cent;
      dijetcore::unique_ptr<TTreeReaderValue<double>> embed_rp;
      dijetcore::unique_ptr<TTreeReaderValue<double>> embed_zdcrate;
      dijetcore::unique_ptr<TTreeReaderValue<double>> embed_vz;

      dijetcore::unique_ptr<TTreeReaderValue<bool>> only_found_match;
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pp_only_jl;
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pp_only_js;
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pp_only_jlm;
      dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> pp_only_jsm;

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

      // initialize source dependent reader values

      if (auau_off_axis_present) {
        jloa = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
                                                                        "jloa");
        jsoa = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader,
                                                                        "jsoa");
        jloaconst =
            dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jloaconst");
        jloarho =
            dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jloarho");
        jloasig =
            dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jloasig");
        jsoaconst =
            dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jsoaconst");
        jsoarho =
            dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jsoarho");
        jsoasig =
            dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jsoasig");
      }

      if (pp_embedded_present) {
        embed_eventid = dijetcore::make_unique<TTreeReaderValue<int>>(
            reader, "embed_eventid");
        embed_runid = dijetcore::make_unique<TTreeReaderValue<int>>(
            reader, "embed_runid");
        embed_refmult = dijetcore::make_unique<TTreeReaderValue<int>>(
            reader, "embed_refmult");
        embed_grefmult = dijetcore::make_unique<TTreeReaderValue<int>>(
            reader, "embed_grefmult");
        embed_nprt = dijetcore::make_unique<TTreeReaderValue<int>>(
            reader, "embed_npart");
        embed_refmultcorr = dijetcore::make_unique<TTreeReaderValue<double>>(
            reader, "embed_refmultcorr");
        embed_grefmultcorr = dijetcore::make_unique<TTreeReaderValue<double>>(
            reader, "embed_grefmultcorr");
        embed_cent =
            dijetcore::make_unique<TTreeReaderValue<int>>(reader, "embed_cent");
        embed_rp = dijetcore::make_unique<TTreeReaderValue<double>>(reader,
                                                                    "embed_rp");
        embed_zdcrate = dijetcore::make_unique<TTreeReaderValue<double>>(
            reader, "embed_zdcrate");
        embed_vz = dijetcore::make_unique<TTreeReaderValue<double>>(reader,
                                                                    "embed_vz");
      }

      if (pp_only_present) {
        only_found_match =
            dijetcore::make_unique<TTreeReaderValue<bool>>(reader, "foundpp");
        pp_only_jl = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(
            reader, "ppjl");
        pp_only_js = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(
            reader, "ppjs");
        pp_only_jlm = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(
            reader, "ppjlm");
        pp_only_jsm = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(
            reader, "ppjsm");
      }

      // build prefix for names so histograms don't get confused
      // (root fuckin sucks)
      std::string key_prefix = "key_" + std::to_string(entry) + "_";
      std::string datatype_prefix = TYPENAME[type_enum];
      std::string hist_prefix = key_prefix + datatype_prefix;

      // create the histograms
      hard_lead_eta[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardleadeta").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, -1, 1);
      hard_lead_phi[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardleadphi").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, -TMath::Pi(), TMath::Pi());
      hard_lead_pt[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardleadpt").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, 50);
      hard_lead_const[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardleadconst").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 100);
      hard_lead_rho[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardleadrho").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 100);
      hard_lead_sig[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardleadsig").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 20);

      match_lead_eta[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchleadeta").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, -1, 1);
      match_lead_phi[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchleadphi").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, -TMath::Pi(), TMath::Pi());
      match_lead_pt[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchleadpt").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 50);
      match_lead_const[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchleadconst").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 100);
      match_lead_rho[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchleadrho").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 100);
      match_lead_sig[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchleadsig").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 20);

      hard_sub_eta[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardsubeta").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, -1, 1);
      hard_sub_phi[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardsubphi").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, -TMath::Pi(), TMath::Pi());
      hard_sub_pt[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardsubpt").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, 50);
      hard_sub_const[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardsubconst").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 100);
      hard_sub_rho[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardsubrho").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, 100);
      hard_sub_sig[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "hardsubsig").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, 20);

      match_sub_eta[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchsubeta").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, -1, 1);
      match_sub_phi[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchsubphi").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, -TMath::Pi(), TMath::Pi());
      match_sub_pt[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchsubpt").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, 50);
      match_sub_const[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchsubconst").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 100);
      match_sub_rho[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchsubrho").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 100);
      match_sub_sig[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchsubsig").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, 0, 20);

      npart[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "npart").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   100, 0, 1200);
      hard_dphi[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "harddphi").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, TMath::Pi());
      match_dphi[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "matchdphi").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, TMath::Pi());
      lead_dr[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "leaddr").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, 0.5);
      sub_dr[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "subdr").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, 0, 0.5);
      lead_dpt[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "leaddpt").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, -20, 20);
      sub_dpt[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "subdptfrac").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, -20, 20);
      lead_dpt_frac[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "leaddptfrac").c_str(),
                   "", cent_boundaries.size(), -0.5,
                   cent_boundaries.size() - 0.5, 50, -2.0, 2.0);
      sub_dpt_frac[data_index][key] =
          new TH2D(dijetcore::MakeString(hist_prefix, "subdpt").c_str(), "",
                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                   50, -2.0, 2.0);

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

      if (auau_off_axis_present) {
        off_axis_aj[data_index][key] =
            new TH2D(dijetcore::MakeString(hist_prefix, "offaxisaj").c_str(),
                     "", cent_boundaries.size(), -0.5,
                     cent_boundaries.size() - 0.5, 15, 0.0001, 0.9);
        off_axis_aj_test[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "offaxisajtest").c_str(), "",
            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000,
            0.0001, 0.9);
        off_axis_rho[data_index][key] =
            new TH2D(dijetcore::MakeString(hist_prefix, "offaxisrho").c_str(),
                     "", cent_boundaries.size(), -0.5,
                     cent_boundaries.size() - 0.5, 50, 0, 100);
        off_axis_sig[data_index][key] =
            new TH2D(dijetcore::MakeString(hist_prefix, "offaxissig").c_str(),
                     "", cent_boundaries.size(), -0.5,
                     cent_boundaries.size() - 0.5, 50, 0, 20);
      }
      if (pp_only_present) {
        hard_pp_only_aj[data_index][key] =
            new TH2D(dijetcore::MakeString(hist_prefix, "pphardaj").c_str(), "",
                     cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                     15, 0.0001, 0.9);
        match_pp_only_aj[data_index][key] =
            new TH2D(dijetcore::MakeString(hist_prefix, "ppmatchaj").c_str(),
                     "", cent_boundaries.size(), -0.5,
                     cent_boundaries.size() - 0.5, 15, 0.0001, 0.9);
        hard_lead_pp_dpt[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "pphardleaddpt").c_str(), "",
            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 50, 0,
            30);
        match_lead_pp_dpt[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "ppmatchleaddpt").c_str(), "",
            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 50, 0,
            30);
        hard_sub_pp_dpt[data_index][key] =
            new TH2D(dijetcore::MakeString(hist_prefix, "pphardsubdpt").c_str(),
                     "", cent_boundaries.size(), -0.5,
                     cent_boundaries.size() - 0.5, 50, 0, 30);
        match_sub_pp_dpt[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "ppmatchsubdpt").c_str(), "",
            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 50, 0,
            30);
        hard_lead_pp_dpt_frac[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "pphardleaddptfrac").c_str(), "",
            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 50, 0,
            30);
        match_lead_pp_dpt_frac[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "ppmatchleaddptfrac").c_str(),
            "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 50,
            0, 30);
        hard_sub_pp_dpt_frac[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "pphardsubdptfrac").c_str(), "",
            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 50, 0,
            30);
        match_sub_pp_dpt_frac[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "ppmatchsubdptfrac").c_str(), "",
            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 50, 0,
            30);
        pp_lead_hard_match_fraction[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "ppleadhardmatchfrac").c_str(),
            "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 2,
            -0.5, 1.5);
        pp_sub_hard_match_fraction[data_index][key] = new TH2D(
            dijetcore::MakeString(hist_prefix, "ppsubhardmatchfrac").c_str(),
            "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 2,
            -0.5, 1.5);
      }

      // event loop

      while (reader.Next()) {
        int cent_bin = -1;
        for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
          auto& pair = cent_bin_boundaries[i];
          if (**cent >= pair.first && **cent <= pair.second) {
            cent_bin = i;
            break;
          }
        }

        if (cent_bin == -1) continue;

        // check jet eta
        if (FLAGS_setScanning &&
            (fabs((*jl)->Eta()) > max_eta || fabs((*js)->Eta()) > max_eta))
          continue;

        hard_lead_eta[data_index][key]->Fill(cent_bin, (*jl)->Eta());
        hard_lead_phi[data_index][key]->Fill(cent_bin, (*jl)->Phi());
        hard_lead_pt[data_index][key]->Fill(cent_bin, (*jl)->Pt());
        hard_lead_const[data_index][key]->Fill(cent_bin, **jlconst);
        hard_lead_rho[data_index][key]->Fill(cent_bin, **jlrho);
        hard_lead_sig[data_index][key]->Fill(cent_bin, **jlsig);
        match_lead_eta[data_index][key]->Fill(cent_bin, (*jlm)->Eta());
        match_lead_phi[data_index][key]->Fill(cent_bin, (*jlm)->Phi());
        match_lead_pt[data_index][key]->Fill(cent_bin, (*jlm)->Pt());
        match_lead_const[data_index][key]->Fill(cent_bin, **jlmconst);
        match_lead_rho[data_index][key]->Fill(cent_bin, **jlmrho);
        match_lead_sig[data_index][key]->Fill(cent_bin, **jlmsig);

        hard_sub_eta[data_index][key]->Fill(cent_bin, (*js)->Eta());
        hard_sub_phi[data_index][key]->Fill(cent_bin, (*js)->Phi());
        hard_sub_pt[data_index][key]->Fill(cent_bin, (*js)->Pt());
        hard_sub_const[data_index][key]->Fill(cent_bin, **jsconst);
        hard_sub_rho[data_index][key]->Fill(cent_bin, **jsrho);
        hard_sub_sig[data_index][key]->Fill(cent_bin, **jssig);
        match_sub_eta[data_index][key]->Fill(cent_bin, (*jsm)->Eta());
        match_sub_phi[data_index][key]->Fill(cent_bin, (*jsm)->Phi());
        match_sub_pt[data_index][key]->Fill(cent_bin, (*jsm)->Pt());
        match_sub_const[data_index][key]->Fill(cent_bin, **jsmconst);
        match_sub_rho[data_index][key]->Fill(cent_bin, **jsmrho);
        match_sub_sig[data_index][key]->Fill(cent_bin, **jsmsig);

        npart[data_index][key]->Fill(cent_bin, **nprt);
        hard_dphi[data_index][key]->Fill(cent_bin, fabs((*jl)->DeltaPhi(**js)));
        match_dphi[data_index][key]->Fill(cent_bin,
                                          fabs((*jlm)->DeltaPhi(**jsm)));
        lead_dr[data_index][key]->Fill(cent_bin, (*jl)->DeltaR(**jlm));
        sub_dr[data_index][key]->Fill(cent_bin, (*js)->DeltaR(**jsm));
        lead_dpt[data_index][key]->Fill(cent_bin, (*jl)->Pt() - (*jlm)->Pt());
        sub_dpt[data_index][key]->Fill(cent_bin, (*js)->Pt() - (*jsm)->Pt());
        lead_dpt_frac[data_index][key]->Fill(
            cent_bin, ((*jl)->Pt() - (*jlm)->Pt()) / (*jl)->Pt());
        sub_dpt_frac[data_index][key]->Fill(
            cent_bin, ((*js)->Pt() - (*jsm)->Pt()) / (*js)->Pt());

        hard_aj[data_index][key]->Fill(
            cent_bin,
            fabs((*jl)->Pt() - (*js)->Pt()) / ((*jl)->Pt() + (*js)->Pt()));
        match_aj[data_index][key]->Fill(
            cent_bin,
            fabs((*jlm)->Pt() - (*jsm)->Pt()) / ((*jlm)->Pt() + (*jsm)->Pt()));
        hard_aj_test[data_index][key]->Fill(
            cent_bin,
            fabs((*jl)->Pt() - (*js)->Pt()) / ((*jl)->Pt() + (*js)->Pt()));
        match_aj_test[data_index][key]->Fill(
            cent_bin,
            fabs((*jlm)->Pt() - (*jsm)->Pt()) / ((*jlm)->Pt() + (*jsm)->Pt()));

        if (auau_off_axis_present) {
          off_axis_aj[data_index][key]->Fill(
              cent_bin, fabs((*jloa)->Pt() - (*jsoa)->Pt()) /
                            ((*jloa)->Pt() + (*jsoa)->Pt()));
          off_axis_aj_test[data_index][key]->Fill(
              cent_bin, fabs((*jloa)->Pt() - (*jsoa)->Pt()) /
                            ((*jloa)->Pt() + (*jsoa)->Pt()));
          off_axis_rho[data_index][key]->Fill(cent_bin, **jloarho);
          off_axis_sig[data_index][key]->Fill(cent_bin, **jloasig);
        }

        if (pp_only_present) {
          if ((*pp_only_jl)->Pt() > 0) {
            pp_lead_hard_match_fraction[data_index][key]->Fill(cent_bin, 1);
            if ((*pp_only_js)->Pt() > 0) {
              pp_sub_hard_match_fraction[data_index][key]->Fill(cent_bin, 1);
              hard_pp_only_aj[data_index][key]->Fill(
                  cent_bin, fabs((*pp_only_jl)->Pt() - (*pp_only_js)->Pt()) /
                                ((*pp_only_jl)->Pt() + (*pp_only_js)->Pt()));
              match_pp_only_aj[data_index][key]->Fill(
                  cent_bin, fabs((*pp_only_jlm)->Pt() - (*pp_only_jsm)->Pt()) /
                                ((*pp_only_jlm)->Pt() + (*pp_only_jsm)->Pt()));
              hard_lead_pp_dpt[data_index][key]->Fill(
                  cent_bin, (*jl)->Pt() - (*pp_only_jl)->Pt());
              match_lead_pp_dpt[data_index][key]->Fill(
                  cent_bin, (*jlm)->Pt() - (*pp_only_jlm)->Pt());
              hard_sub_pp_dpt[data_index][key]->Fill(
                  cent_bin, (*js)->Pt() - (*pp_only_js)->Pt());
              match_sub_pp_dpt[data_index][key]->Fill(
                  cent_bin, (*jsm)->Pt() - (*pp_only_jsm)->Pt());
              hard_lead_pp_dpt_frac[data_index][key]->Fill(
                  cent_bin, ((*jl)->Pt() - (*pp_only_jl)->Pt()) / (*jl)->Pt());
              match_lead_pp_dpt_frac[data_index][key]->Fill(
                  cent_bin,
                  ((*jlm)->Pt() - (*pp_only_jlm)->Pt()) / (*jlm)->Pt());
              hard_sub_pp_dpt_frac[data_index][key]->Fill(
                  cent_bin, ((*js)->Pt() - (*pp_only_js)->Pt()) / (*js)->Pt());
              match_sub_pp_dpt_frac[data_index][key]->Fill(
                  cent_bin,
                  ((*jsm)->Pt() - (*pp_only_jsm)->Pt()) / (*jsm)->Pt());
            } else {
              pp_sub_hard_match_fraction[data_index][key]->Fill(cent_bin, 0);
            }
          } else
            pp_lead_hard_match_fraction[data_index][key]->Fill(cent_bin, 0);
        }
      }

      hard_lead_eta_cent[data_index][key] =
          SplitByCentralityNormalized(hard_lead_eta[data_index][key]);
      hard_lead_phi_cent[data_index][key] =
          SplitByCentralityNormalized(hard_lead_phi[data_index][key]);
      hard_lead_pt_cent[data_index][key] =
          SplitByCentralityNormalized(hard_lead_pt[data_index][key]);
      hard_lead_const_cent[data_index][key] =
          SplitByCentralityNormalized(hard_lead_const[data_index][key]);
      hard_lead_rho_cent[data_index][key] =
          SplitByCentralityNormalized(hard_lead_rho[data_index][key]);
      hard_lead_sig_cent[data_index][key] =
          SplitByCentralityNormalized(hard_lead_sig[data_index][key]);
      match_lead_eta_cent[data_index][key] =
          SplitByCentralityNormalized(match_lead_eta[data_index][key]);
      match_lead_phi_cent[data_index][key] =
          SplitByCentralityNormalized(match_lead_phi[data_index][key]);
      match_lead_pt_cent[data_index][key] =
          SplitByCentralityNormalized(match_lead_pt[data_index][key]);
      match_lead_const_cent[data_index][key] =
          SplitByCentralityNormalized(match_lead_const[data_index][key]);
      match_lead_rho_cent[data_index][key] =
          SplitByCentralityNormalized(match_lead_rho[data_index][key]);
      match_lead_sig_cent[data_index][key] =
          SplitByCentralityNormalized(match_lead_sig[data_index][key]);

      hard_sub_eta_cent[data_index][key] =
          SplitByCentralityNormalized(hard_sub_eta[data_index][key]);
      hard_sub_phi_cent[data_index][key] =
          SplitByCentralityNormalized(hard_sub_phi[data_index][key]);
      hard_sub_pt_cent[data_index][key] =
          SplitByCentralityNormalized(hard_sub_pt[data_index][key]);
      hard_sub_const_cent[data_index][key] =
          SplitByCentralityNormalized(hard_sub_const[data_index][key]);
      hard_sub_rho_cent[data_index][key] =
          SplitByCentralityNormalized(hard_sub_rho[data_index][key]);
      hard_sub_sig_cent[data_index][key] =
          SplitByCentralityNormalized(hard_sub_sig[data_index][key]);
      match_sub_eta_cent[data_index][key] =
          SplitByCentralityNormalized(match_sub_eta[data_index][key]);
      match_sub_phi_cent[data_index][key] =
          SplitByCentralityNormalized(match_sub_phi[data_index][key]);
      match_sub_pt_cent[data_index][key] =
          SplitByCentralityNormalized(match_sub_pt[data_index][key]);
      match_sub_const_cent[data_index][key] =
          SplitByCentralityNormalized(match_sub_const[data_index][key]);
      match_sub_rho_cent[data_index][key] =
          SplitByCentralityNormalized(match_sub_rho[data_index][key]);
      match_sub_sig_cent[data_index][key] =
          SplitByCentralityNormalized(match_sub_sig[data_index][key]);

      npart_cent[data_index][key] =
          SplitByCentralityNormalized(npart[data_index][key]);
      hard_dphi_cent[data_index][key] =
          SplitByCentralityNormalized(hard_dphi[data_index][key]);
      match_dphi_cent[data_index][key] =
          SplitByCentralityNormalized(match_dphi[data_index][key]);
      lead_dr_cent[data_index][key] =
          SplitByCentralityNormalized(lead_dr[data_index][key]);
      sub_dr_cent[data_index][key] =
          SplitByCentralityNormalized(sub_dr[data_index][key]);
      lead_dpt_cent[data_index][key] =
          SplitByCentralityNormalized(lead_dpt[data_index][key]);
      sub_dpt_cent[data_index][key] =
          SplitByCentralityNormalized(sub_dpt[data_index][key]);
      lead_dpt_frac_cent[data_index][key] =
          SplitByCentralityNormalized(lead_dpt_frac[data_index][key]);
      sub_dpt_frac_cent[data_index][key] =
          SplitByCentralityNormalized(sub_dpt_frac[data_index][key]);

      hard_aj_cent[data_index][key] =
          SplitByCentralityNormalized(hard_aj[data_index][key]);
      match_aj_cent[data_index][key] =
          SplitByCentralityNormalized(match_aj[data_index][key]);
      hard_aj_test_cent[data_index][key] =
          SplitByCentralityNormalized(hard_aj_test[data_index][key]);
      match_aj_test_cent[data_index][key] =
          SplitByCentralityNormalized(match_aj_test[data_index][key]);

      if (auau_off_axis_present) {
        off_axis_aj_cent[data_index][key] =
            SplitByCentralityNormalized(off_axis_aj[data_index][key]);
        off_axis_aj_test_cent[data_index][key] =
            SplitByCentralityNormalized(off_axis_aj_test[data_index][key]);
        off_axis_rho_cent[data_index][key] =
            SplitByCentralityNormalized(off_axis_rho[data_index][key]);
        off_axis_sig_cent[data_index][key] =
            SplitByCentralityNormalized(off_axis_sig[data_index][key]);
      }

      if (pp_only_present) {
        hard_pp_only_aj_cent[data_index][key] =
            SplitByCentralityNormalized(hard_pp_only_aj[data_index][key]);
        match_pp_only_aj_cent[data_index][key] =
            SplitByCentralityNormalized(match_pp_only_aj[data_index][key]);
        hard_lead_pp_dpt_cent[data_index][key] =
            SplitByCentralityNormalized(hard_lead_pp_dpt[data_index][key]);
        match_lead_pp_dpt_cent[data_index][key] =
            SplitByCentralityNormalized(match_lead_pp_dpt[data_index][key]);
        hard_sub_pp_dpt_cent[data_index][key] =
            SplitByCentralityNormalized(hard_sub_pp_dpt[data_index][key]);
        match_sub_pp_dpt_cent[data_index][key] =
            SplitByCentralityNormalized(match_sub_pp_dpt[data_index][key]);
        hard_lead_pp_dpt_frac_cent[data_index][key] =
            SplitByCentralityNormalized(hard_lead_pp_dpt_frac[data_index][key]);
        match_lead_pp_dpt_frac_cent[data_index][key] =
            SplitByCentralityNormalized(
                match_lead_pp_dpt_frac[data_index][key]);
        hard_sub_pp_dpt_frac_cent[data_index][key] =
            SplitByCentralityNormalized(hard_sub_pp_dpt_frac[data_index][key]);
        match_sub_pp_dpt_frac_cent[data_index][key] =
            SplitByCentralityNormalized(match_sub_pp_dpt_frac[data_index][key]);
        pp_lead_hard_match_fraction_cent[data_index][key] =
            SplitByCentralityNormalized(
                pp_lead_hard_match_fraction[data_index][key]);
        pp_sub_hard_match_fraction_cent[data_index][key] =
            SplitByCentralityNormalized(
                pp_sub_hard_match_fraction[data_index][key]);
      }

      if (auau_off_axis_present) {
        Overlay1D(off_axis_aj_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "off_axis_aj", "", "|A_{J}|",
                  "event fraction", "Centrality");
        Overlay1D(off_axis_rho_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "off_axis_rho", "", "rho", "event fraction",
                  "Centrality");
        Overlay1D(off_axis_sig_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "off_axis_sig", "", "sig", "event fraction",
                  "Centrality");
      }

      if (pp_only_present) {
        Overlay1D(hard_pp_only_aj_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "pp_only_hard_aj", "", "|A_{J}|",
                  "event fraction", "Centrality");
        Overlay1D(match_pp_only_aj_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "pp_only_match_aj", "", "|A_{J}|",
                  "event fraction", "Centrality");
        Overlay1D(hard_lead_pp_dpt_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "pp_only_hard_lead_dpt", "", "dp_{T}",
                  "event fraction", "Centrality");
        Overlay1D(hard_sub_pp_dpt_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "pp_only_hard_sub_dpt", "", "dp_{T}",
                  "event fraction", "Centrality");
        Overlay1D(match_lead_pp_dpt_cent[data_index][key], refcent_string,
                  hopts, copts, out_loc, "pp_only_match_lead_dpt", "", "dp_{T}",
                  "event fraction", "Centrality");
        Overlay1D(match_sub_pp_dpt_cent[data_index][key], refcent_string, hopts,
                  copts, out_loc, "pp_only_match_sub_dpt", "", "dp_{T}",
                  "event fraction", "Centrality");
        Overlay1D(hard_lead_pp_dpt_frac_cent[data_index][key], refcent_string,
                  hopts, copts, out_loc, "pp_only_hard_lead_dpt_frac", "",
                  "dp_{T}/p_{T}", "event fraction", "Centrality");
        Overlay1D(hard_sub_pp_dpt_frac_cent[data_index][key], refcent_string,
                  hopts, copts, out_loc, "pp_only_hard_sub_dpt_frac", "",
                  "dp_{T}/p_{T}", "event fraction", "Centrality");
        Overlay1D(match_lead_pp_dpt_frac_cent[data_index][key], refcent_string,
                  hopts, copts, out_loc, "pp_only_match_lead_dpt_frac", "",
                  "dp_{T}/p_{T}", "event fraction", "Centrality");
        Overlay1D(match_sub_pp_dpt_frac_cent[data_index][key], refcent_string,
                  hopts, copts, out_loc, "pp_only_match_sub_dpt_frac", "",
                  "dp_{T}/p_{T}", "event fraction", "Centrality");
        Overlay1D(pp_lead_hard_match_fraction_cent[data_index][key],
                  refcent_string, hopts, copts, out_loc,
                  "pp_lead_hard_match_fraction", "", "", "event fraction",
                  "Centrality");
        Overlay1D(pp_sub_hard_match_fraction_cent[data_index][key],
                  refcent_string, hopts, copts, out_loc,
                  "pp_sub_hard_match_fraction", "", "", "event fraction",
                  "Centrality");
      }
    }  // datatype

    // now we will get systematics and print out comparisons - do everything in
    // bins of centrality
    for (int i = 0; i < refcent_string.size(); ++i) {
      // make our output directory
      string out_loc = key_loc + "/cent_" + refcent_string_fs[i];
      boost::filesystem::path dir(out_loc.c_str());
      boost::filesystem::create_directories(dir);

      // get systematic errors
      systematic_errors_hard[key].push_back(
          GetSystematic(hard_aj_cent[pp_index][key][i],
                        hard_aj_cent[TYPEMAP[DATATYPE::TOWP]][key][i],
                        hard_aj_cent[TYPEMAP[DATATYPE::TOWM]][key][i],
                        hard_aj_cent[TYPEMAP[DATATYPE::TRACKP]][key][i],
                        hard_aj_cent[TYPEMAP[DATATYPE::TRACKM]][key][i]));
      systematic_errors_match[key].push_back(
          GetSystematic(match_aj_cent[pp_index][key][i],
                        match_aj_cent[TYPEMAP[DATATYPE::TOWP]][key][i],
                        match_aj_cent[TYPEMAP[DATATYPE::TOWM]][key][i],
                        match_aj_cent[TYPEMAP[DATATYPE::TRACKP]][key][i],
                        match_aj_cent[TYPEMAP[DATATYPE::TRACKM]][key][i]));

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
      AjPrintout(hard_aj_cent[auau_index][key][i],
                 hard_aj_cent[pp_index][key][i], systematic_errors_hard[key][i],
                 0.0, 0.25, 0.0001, 0.9, hardPave, "Au+Au HT",
                 "p+p HT #oplus Au+Au MB", hopts, copts, out_loc, "aj_hard", "",
                 "|A_{J}|", "event fraction");
      AjPrintout(match_aj_cent[auau_index][key][i],
                 match_aj_cent[pp_index][key][i],
                 systematic_errors_match[key][i], 0.0, 0.3, 0.0001, 0.9,
                 matchPave, "Au+Au HT", "p+p HT #oplus Au+Au MB", hopts, copts,
                 out_loc, "aj_match", "", "|A_{J}|", "event fraction");
      AjPrintout(match_aj_cent[auau_index][key][i],
                 off_axis_aj_cent[auau_index][key][i], nullptr, 0.0, 0.3, 0.0,
                 0.9, matchPave, "Au+Au HT", "Au+Au HC embedded", hopts, copts,
                 out_loc, "aj_embed", "", "|A_{J}|", "event fraction");

      // now print off-axis AJ with the matched
      std::vector<TH1D*> off_axis_match_compare{
          match_aj_cent[auau_index][key][i], match_aj_cent[pp_index][key][i],
          off_axis_aj_cent[auau_index][key][i]};
      std::vector<string> off_axis_match_string{"Au+Au", "embedded p+p",
                                                "embedded HC Au+Au"};

      // now we will compare all the other observables
      Overlay1D(npart_cent[auau_index][key][i], npart_cent[pp_index][key][i],
                "Au+Au", "p+p embedded", hopts, copts, out_loc, "npart", "",
                "N_{part}", "event fraction");
      Overlay1D(hard_dphi_cent[auau_index][key][i],
                hard_dphi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_dphi", "", "d#phi",
                "event fraction");
      Overlay1D(match_dphi_cent[auau_index][key][i],
                match_dphi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_dphi", "", "d#phi",
                "event fraction");
      Overlay1D(lead_dr_cent[auau_index][key][i],
                lead_dr_cent[pp_index][key][i], "Au+Au", "p+p embedded", hopts,
                copts, out_loc, "lead_dr", "", "d#phi", "event fraction");
      Overlay1D(sub_dr_cent[auau_index][key][i], sub_dr_cent[pp_index][key][i],
                "Au+Au", "p+p embedded", hopts, copts, out_loc, "sub_dr", "",
                "d#phi", "event fraction");
      Overlay1D(lead_dpt_cent[auau_index][key][i],
                lead_dpt_cent[pp_index][key][i], "Au+Au", "p+p embedded", hopts,
                copts, out_loc, "lead_dpt", "", "dp_{T}", "event fraction");
      Overlay1D(sub_dpt_cent[auau_index][key][i],
                sub_dpt_cent[pp_index][key][i], "Au+Au", "p+p embedded", hopts,
                copts, out_loc, "sub_dpt", "", "dp_{T}", "event fraction");
      Overlay1D(lead_dpt_frac_cent[auau_index][key][i],
                lead_dpt_frac_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "lead_dpt_frac", "", "dp_{T}",
                "event fraction");
      Overlay1D(sub_dpt_frac_cent[auau_index][key][i],
                sub_dpt_frac_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "sub_dpt_frac", "", "dp_{T}",
                "event fraction");
      Overlay1D(hard_lead_eta_cent[auau_index][key][i],
                hard_lead_eta_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_lead_eta", "", "#eta",
                "event fraction");
      Overlay1D(match_lead_eta_cent[auau_index][key][i],
                match_lead_eta_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_lead_eta", "", "#eta",
                "event fraction");
      Overlay1D(hard_sub_eta_cent[auau_index][key][i],
                hard_sub_eta_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_sub_eta", "", "#eta",
                "event fraction");
      Overlay1D(match_sub_eta_cent[auau_index][key][i],
                match_sub_eta_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_sub_eta", "", "#eta",
                "event fraction");
      Overlay1D(hard_lead_phi_cent[auau_index][key][i],
                hard_lead_phi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_lead_phi", "", "#phi",
                "event fraction");
      Overlay1D(match_lead_phi_cent[auau_index][key][i],
                match_lead_phi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_lead_phi", "", "#phi",
                "event fraction");
      Overlay1D(hard_sub_phi_cent[auau_index][key][i],
                hard_sub_phi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_sub_phi", "", "#phi",
                "event fraction");
      Overlay1D(match_sub_phi_cent[auau_index][key][i],
                match_sub_phi_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_sub_phi", "", "#phi",
                "event fraction");
      Overlay1D(hard_lead_pt_cent[auau_index][key][i],
                hard_lead_pt_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_lead_pt", "", "p_{T}",
                "event fraction");
      Overlay1D(match_lead_pt_cent[auau_index][key][i],
                match_lead_pt_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_lead_pt", "", "p_{T}",
                "event fraction");
      Overlay1D(hard_sub_pt_cent[auau_index][key][i],
                hard_sub_pt_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_sub_pt", "", "p_{T}",
                "event fraction");
      Overlay1D(match_sub_pt_cent[auau_index][key][i],
                match_sub_pt_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_sub_pt", "", "p_{T}",
                "event fraction");
      Overlay1D(hard_lead_const_cent[auau_index][key][i],
                hard_lead_const_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_lead_const", "", "N_{part}",
                "event fraction");
      Overlay1D(match_lead_const_cent[auau_index][key][i],
                match_lead_const_cent[pp_index][key][i], "Au+Au",
                "p+p embedded", hopts, copts, out_loc, "match_lead_const", "",
                "N_{part}", "event fraction");
      Overlay1D(hard_sub_const_cent[auau_index][key][i],
                hard_sub_const_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_sub_const", "", "N_{part}",
                "event fraction");
      Overlay1D(match_sub_const_cent[auau_index][key][i],
                match_sub_const_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_sub_const", "", "N_{part}",
                "event fraction");
      Overlay1D(hard_lead_rho_cent[auau_index][key][i],
                hard_lead_rho_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_lead_rho", "", "rho",
                "event fraction");
      Overlay1D(match_lead_rho_cent[auau_index][key][i],
                match_lead_rho_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_lead_rho", "", "rho",
                "event fraction");
      Overlay1D(hard_sub_rho_cent[auau_index][key][i],
                hard_sub_rho_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_sub_rho", "", "rho",
                "event fraction");
      Overlay1D(match_sub_rho_cent[auau_index][key][i],
                match_sub_rho_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_sub_rho", "", "rho",
                "event fraction");
      Overlay1D(hard_lead_sig_cent[auau_index][key][i],
                hard_lead_sig_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_lead_sig", "", "sigma",
                "event fraction");
      Overlay1D(match_lead_sig_cent[auau_index][key][i],
                match_lead_sig_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_lead_sig", "", "sigma",
                "event fraction");
      Overlay1D(hard_sub_sig_cent[auau_index][key][i],
                hard_sub_sig_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "hard_sub_sig", "", "sigma",
                "event fraction");
      Overlay1D(match_sub_sig_cent[auau_index][key][i],
                match_sub_sig_cent[pp_index][key][i], "Au+Au", "p+p embedded",
                hopts, copts, out_loc, "match_sub_sig", "", "sigma",
                "event fraction");

    }  // centrality

  }  // key

  for (int cent = 0; cent < cent_boundaries.size(); ++cent) {
    // make a directory for our grid outputs
    string out_loc_grid =
        FLAGS_outputDir + "/grid_cent_" + refcent_string_fs[cent];
    LOG(INFO) << "output directory: " << out_loc_grid;
    // build output directory if it doesn't exist, using boost::filesystem
    boost::filesystem::path dir_cent(out_loc_grid.c_str());
    boost::filesystem::create_directories(dir_cent);

    // now build histograms
    TH2D* p_values_hard =
        new TH2D(dijetcore::MakeString("p_values_hard_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_hard =
        new TH2D(dijetcore::MakeString("ks_values_hard_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* p_values_match =
        new TH2D(dijetcore::MakeString("p_values_match_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_match =
        new TH2D(dijetcore::MakeString("ks_values_match_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* p_values_bkg =
        new TH2D(dijetcore::MakeString("p_values_bkg_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_bkg =
        new TH2D(dijetcore::MakeString("ks_values_bkg_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);

    TH2D* p_values_hard_error =
        new TH2D(dijetcore::MakeString("p_values_hard_error_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_hard_error =
        new TH2D(dijetcore::MakeString("ks_values_hard_error_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* p_values_match_error =
        new TH2D(dijetcore::MakeString("p_values_match_error_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_match_error =
        new TH2D(dijetcore::MakeString("ks_values_match_error_", cent).c_str(),
                 ";R;p_{T}^{const}", radii.size(), -0.5, radii.size() - 0.5,
                 constpt.size(), -0.5, constpt.size() - 0.5);

    for (int j = 0; j < radii_string.size(); ++j) {
      p_values_hard->GetXaxis()->SetBinLabel(j + 1, radii_string[j].c_str());
      ks_values_hard->GetXaxis()->SetBinLabel(j + 1, radii_string[j].c_str());
      p_values_match->GetXaxis()->SetBinLabel(j + 1, radii_string[j].c_str());
      ks_values_match->GetXaxis()->SetBinLabel(j + 1, radii_string[j].c_str());
      p_values_bkg->GetXaxis()->SetBinLabel(j + 1, radii_string[j].c_str());
      ks_values_bkg->GetXaxis()->SetBinLabel(j + 1, radii_string[j].c_str());
      p_values_hard_error->GetXaxis()->SetBinLabel(j + 1,
                                                   radii_string[j].c_str());
      ks_values_hard_error->GetXaxis()->SetBinLabel(j + 1,
                                                    radii_string[j].c_str());
      p_values_match_error->GetXaxis()->SetBinLabel(j + 1,
                                                    radii_string[j].c_str());
      ks_values_match_error->GetXaxis()->SetBinLabel(j + 1,
                                                     radii_string[j].c_str());
    }
    for (int j = 0; j < constpt_string.size(); ++j) {
      p_values_hard->GetYaxis()->SetBinLabel(j + 1, constpt_string[j].c_str());
      ks_values_hard->GetYaxis()->SetBinLabel(j + 1, constpt_string[j].c_str());
      p_values_match->GetYaxis()->SetBinLabel(j + 1, constpt_string[j].c_str());
      ks_values_match->GetYaxis()->SetBinLabel(j + 1,
                                               constpt_string[j].c_str());
      p_values_bkg->GetYaxis()->SetBinLabel(j + 1, constpt_string[j].c_str());
      ks_values_bkg->GetYaxis()->SetBinLabel(j + 1, constpt_string[j].c_str());
      p_values_hard_error->GetYaxis()->SetBinLabel(j + 1,
                                                   constpt_string[j].c_str());
      ks_values_hard_error->GetYaxis()->SetBinLabel(j + 1,
                                                    constpt_string[j].c_str());
      p_values_match_error->GetYaxis()->SetBinLabel(j + 1,
                                                    constpt_string[j].c_str());
      ks_values_match_error->GetYaxis()->SetBinLabel(j + 1,
                                                     constpt_string[j].c_str());
    }

    TCanvas* c_hard =
        new TCanvas(dijetcore::MakeString("hard_canvas_", cent).c_str());
    std::vector<std::vector<TPad*>> hard_pads = dijetcore::CanvasPartition(
        c_hard, radii.size(), constpt.size(), 0.10, 0.10, 0.12, 0.05);

    TCanvas* c_match =
        new TCanvas(dijetcore::MakeString("match_canvas_", cent).c_str());
    std::vector<std::vector<TPad*>> matched_pads = dijetcore::CanvasPartition(
        c_match, radii.size(), constpt.size(), 0.10, 0.10, 0.12, 0.05);

    TCanvas* c_match_oa =
        new TCanvas(dijetcore::MakeString("match_canvas_oa_", cent).c_str());
    std::vector<std::vector<TPad*>> matched_oa_pads =
        dijetcore::CanvasPartition(c_match_oa, radii.size(), constpt.size(),
                                   0.10, 0.10, 0.12, 0.05);

    // for drawing x/y axis labels
    TPad* invis = new TPad("invis_pad", "", 0, 0, 1, 1);
    invis->SetFillStyle(4000);
    TLatex* x_label_text = new TLatex(0.47, 0.03, "|A_{J}|");
    x_label_text->SetTextSize(0.05);
    TLatex* y_label_text = new TLatex(0.05, 0.4, "event fraction");
    y_label_text->SetTextSize(0.05);
    y_label_text->SetTextAngle(90);

    for (int rad = 0; rad < radii.size(); ++rad) {
      for (int pt = 0; pt < constpt.size(); ++pt) {
        int x_bin = p_values_hard->GetXaxis()->FindBin(rad);
        int y_bin = p_values_hard->GetYaxis()->FindBin(pt);
        string key = grid_keys[rad][pt];
        dijetcore::DijetKey key_params = grid_key_params[rad][pt];

        TH1D* aj_hard_test = hard_aj_test_cent[auau_index][key][cent];
        TH1D* aj_match_test = match_aj_test_cent[auau_index][key][cent];

        TH1D* aj_hard_pp_test = hard_aj_test_cent[pp_index][key][cent];
        TH1D* aj_match_pp_test = match_aj_test_cent[pp_index][key][cent];

        TH1D* aj_hard = hard_aj_cent[auau_index][key][cent];
        TH1D* aj_match = match_aj_cent[auau_index][key][cent];

        TH1D* aj_hard_pp = hard_aj_cent[pp_index][key][cent];
        TH1D* aj_match_pp = match_aj_cent[pp_index][key][cent];

        TGraphErrors* hard_sys_error = systematic_errors_hard[key][cent];
        TGraphErrors* match_sys_error = systematic_errors_match[key][cent];

        TH1D* aj_off_axis = off_axis_aj_cent[auau_index][key][cent];
        TH1D* aj_off_axis_test = off_axis_aj_test_cent[auau_index][key][cent];

        // print

        c_hard->cd(0);
        aj_hard->GetXaxis()->SetTitle(
            dijetcore::MakeString("R=", radii_string[rad]).c_str());
        aj_hard->GetXaxis()->SetTitleOffset(1.15);
        aj_hard->GetXaxis()->CenterTitle();
        aj_hard->GetYaxis()->SetTitle(
            dijetcore::MakeString("p_{T}^{const}=", constpt_string[pt])
                .c_str());
        aj_hard->GetYaxis()->SetTitleOffset(1.07);
        if (pt != 0 && pt != constpt.size() -1)
            aj_hard->GetYaxis()->SetTitleSize(aj_hard->GetYaxis()->GetTitleSize() * 1.05);
        if (pt == 0)
            aj_hard->GetYaxis()->SetTitleSize(aj_hard->GetYaxis()->GetTitleSize() * 0.95);
        aj_hard->GetYaxis()->CenterTitle();
        aj_hard->GetYaxis()->SetRangeUser(0.00001, 0.2499);
        aj_hard->SetMarkerSize(0.2);
        aj_hard_pp->SetMarkerSize(0.2);
        hard_sys_error->SetFillColorAlpha(
            hard_aj_cent[pp_index][key][cent]->GetLineColor(), 0.4);
        hard_sys_error->SetFillStyle(1001);
        hard_sys_error->SetLineWidth(0);
        hard_sys_error->SetMarkerSize(0);
        TPad* hard_pad = hard_pads[rad][pt];
        double scale_factor =
            hard_pads[0][0]->GetAbsHNDC() / hard_pad->GetAbsHNDC();
        hard_pad->Draw();
        hard_pad->SetFillStyle(4000);
        hard_pad->SetFrameFillStyle(4000);
        hard_pad->cd();
        aj_hard->GetXaxis()->SetNdivisions(305);
        aj_hard->GetXaxis()->SetLabelSize(0.08);
        aj_hard->GetYaxis()->SetNdivisions(305);
        aj_hard->GetYaxis()->SetLabelSize(0.08);
        if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
          aj_hard->GetYaxis()->SetLabelSize(0.11);

        aj_hard->Draw();
        aj_hard_pp->Draw("SAME");
        hard_sys_error->Draw("9e2SAME");

        // draw the STAR Prelim
        TLegend* leg1 = GetLegend(rad, pt, radii_string.size() - 1,
                                  constpt_string.size() - 1, radii[rad],
                                  constpt[pt], scale_factor);
        leg1->AddEntry(aj_hard, "Au+Au HT", "lep")
            ->SetTextSize(.04 * scale_factor);
        leg1->AddEntry(aj_hard_pp, "p+p HT #oplus Au+Au MB", "lep")
            ->SetTextSize(.04 * scale_factor);
        if (pt == constpt.size() - 1 && rad == 0) leg1->Draw();
        c_hard->cd(0);

        c_match->cd(0);
        aj_match->GetXaxis()->SetTitle(
            dijetcore::MakeString("R=", radii_string[rad]).c_str());
        aj_match->GetXaxis()->SetTitleOffset(1.15);
        aj_match->GetXaxis()->CenterTitle();
        aj_match->GetYaxis()->SetTitle(
            dijetcore::MakeString("p_{T}^{const}=", constpt_string[pt])
                .c_str());
        aj_match->GetYaxis()->SetTitleOffset(1.07);
        if (pt != 0 && pt != constpt.size() -1)
            aj_match->GetYaxis()->SetTitleSize(aj_match->GetYaxis()->GetTitleSize() * 1.05);
        if (pt == 0)
            aj_match->GetYaxis()->SetTitleSize(aj_match->GetYaxis()->GetTitleSize() * 0.95);
        aj_match->GetYaxis()->CenterTitle();

        aj_match->GetYaxis()->SetRangeUser(0.00001, 0.2499);
        aj_match->SetMarkerSize(0.2);
        aj_match_pp->SetMarkerSize(0.2);
        match_sys_error->SetFillColorAlpha(
            match_aj_cent[pp_index][key][cent]->GetLineColor(), 0.4);
        match_sys_error->SetFillStyle(1001);
        match_sys_error->SetLineWidth(0);
        match_sys_error->SetMarkerSize(0);
        TPad* match_pad = matched_pads[rad][pt];
        match_pad->SetFillStyle(4000);
        match_pad->SetFrameFillStyle(4000);
        match_pad->cd();
        aj_match->GetXaxis()->SetNdivisions(305);
        aj_match->GetXaxis()->SetLabelSize(0.08);
        aj_match->GetYaxis()->SetNdivisions(305);
        aj_match->GetYaxis()->SetLabelSize(0.08);
        if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
          aj_match->GetYaxis()->SetLabelSize(0.11);

        aj_match->Draw();
        aj_match_pp->Draw("SAME");
        match_sys_error->Draw("9e2SAME");

        // draw STAR prelim
        TLegend* leg2 = GetLegend(rad, pt, radii_string.size() - 1,
                                  constpt_string.size() - 1, radii[rad],
                                  constpt[pt], scale_factor);
        leg2->AddEntry(aj_match, "Au+Au HT", "lep")
            ->SetTextSize(.04 * scale_factor);
        leg2->AddEntry(aj_match_pp, "p+p HT #oplus Au+Au MB", "lep")
            ->SetTextSize(.04 * scale_factor);
        if (pt == constpt.size() - 1 && rad == 0) leg2->Draw();
        c_match->cd(0);

        c_match_oa->cd(0);
        aj_off_axis->GetXaxis()->SetTitle(
            dijetcore::MakeString("R=", radii_string[rad]).c_str());
        aj_off_axis->GetXaxis()->CenterTitle();
        aj_off_axis->GetXaxis()->SetTitleOffset(1.15);
        aj_off_axis->GetYaxis()->SetTitle(
            dijetcore::MakeString("p_{T}^{const}=", constpt_string[pt])
                .c_str());
        aj_off_axis->GetYaxis()->SetTitleOffset(1.07);
        if (pt != 0 && pt != constpt.size() -1)
            aj_off_axis->GetYaxis()->SetTitleSize(aj_off_axis->GetYaxis()->GetTitleSize() * 1.05);
        if (pt == 0)
            aj_off_axis->GetYaxis()->SetTitleSize(aj_off_axis->GetYaxis()->GetTitleSize() * 0.95);
        
        aj_off_axis->GetYaxis()->CenterTitle();

        aj_off_axis->GetYaxis()->SetRangeUser(0.00001, 0.2499);
        aj_off_axis->SetMarkerSize(0.2);
        aj_off_axis->SetMarkerStyle(23);
        aj_off_axis->SetMarkerColor(kAzure);
        aj_off_axis->SetLineColor(kAzure);
        aj_off_axis->SetLineWidth(2);
        aj_off_axis->SetFillColorAlpha(kAzure, 0.65);
        aj_off_axis->SetLineWidth(0);
        aj_off_axis->SetMarkerSize(0);
        TPad* match_oa_pad = matched_oa_pads[rad][pt];
        match_oa_pad->SetFillStyle(4000);
        match_oa_pad->SetFrameFillStyle(4000);
        match_oa_pad->cd();
        aj_off_axis->GetXaxis()->SetNdivisions(305);
        aj_off_axis->GetXaxis()->SetLabelSize(0.08);
        aj_off_axis->GetYaxis()->SetNdivisions(305);
        aj_off_axis->GetYaxis()->SetLabelSize(0.08);
        if (rad == 0 && !(pt == 0 || pt == constpt.size() - 1))
          aj_off_axis->GetYaxis()->SetLabelSize(0.11);

        aj_off_axis->Draw("9e3");
        aj_match->Draw("9SAME");

        // draw STAR prelim
        TLegend* leg3 = GetLegend(rad, pt, radii_string.size() - 1,
                                  constpt_string.size() - 1, radii[rad],
                                  constpt[pt], scale_factor);
        leg3->AddEntry(aj_match, "Au+Au HT", "lep")
            ->SetTextSize(.04 * scale_factor);
        leg3->AddEntry(aj_off_axis, "Au+Au HT #oplus Au+Au MB", "f")
            ->SetTextSize(.04 * scale_factor);
        if (pt == constpt.size() - 1 && rad == 0) leg3->Draw();

        c_match_oa->cd(0);

        c_match_oa->cd(0);

        std::vector<TH1D*> aj_hard_pp_var{
            hard_aj_cent[TYPEMAP[DATATYPE::TOWP]][key][cent],
            hard_aj_cent[TYPEMAP[DATATYPE::TOWM]][key][cent],
            hard_aj_cent[TYPEMAP[DATATYPE::TRACKP]][key][cent],
            hard_aj_cent[TYPEMAP[DATATYPE::TRACKM]][key][cent]};
        std::vector<TH1D*> aj_match_pp_var{
            match_aj_cent[TYPEMAP[DATATYPE::TOWP]][key][cent],
            match_aj_cent[TYPEMAP[DATATYPE::TOWM]][key][cent],
            match_aj_cent[TYPEMAP[DATATYPE::TRACKP]][key][cent],
            match_aj_cent[TYPEMAP[DATATYPE::TRACKM]][key][cent]};

        std::vector<TH1D*> aj_hard_pp_var_test{
            hard_aj_test_cent[TYPEMAP[DATATYPE::TOWP]][key][cent],
            hard_aj_test_cent[TYPEMAP[DATATYPE::TOWM]][key][cent],
            hard_aj_test_cent[TYPEMAP[DATATYPE::TRACKP]][key][cent],
            hard_aj_test_cent[TYPEMAP[DATATYPE::TRACKM]][key][cent]};
        std::vector<TH1D*> aj_match_pp_var_test{
            match_aj_test_cent[TYPEMAP[DATATYPE::TOWP]][key][cent],
            match_aj_test_cent[TYPEMAP[DATATYPE::TOWM]][key][cent],
            match_aj_test_cent[TYPEMAP[DATATYPE::TRACKP]][key][cent],
            match_aj_test_cent[TYPEMAP[DATATYPE::TRACKM]][key][cent]};

        double p_value_hard = aj_hard->Chi2Test(aj_hard_pp, "UU NORM");
        double p_value_match = aj_match->Chi2Test(aj_match_pp, "UU NORM");
        double p_value_bkg = aj_match->Chi2Test(aj_off_axis, "UU NORM");
        double ks_value_hard = aj_hard_test->KolmogorovTest(aj_hard_pp_test);
        double ks_value_match = aj_match_test->KolmogorovTest(aj_match_pp_test);
        double ks_value_bkg = aj_match_test->KolmogorovTest(aj_off_axis_test);

        double max_hard_p_deviation = 0;
        double max_match_p_deviation = 0;
        double max_hard_ks_deviation = 0;
        double max_match_ks_deviation = 0;

        for (int i = 0; i < aj_hard_pp_var.size(); ++i) {
          double hard_p_deviation =
              aj_hard->Chi2Test(aj_hard_pp_var[i], "UU NORM");
          double hard_ks_deviation =
              aj_hard_test->KolmogorovTest(aj_hard_pp_var_test[i]);
          double match_p_deviation =
              aj_match->Chi2Test(aj_match_pp_var[i], "UU NORM");
          double match_ks_deviation =
              aj_match_test->KolmogorovTest(aj_match_pp_var_test[i]);

          if (fabs(p_value_hard - hard_p_deviation) > max_hard_p_deviation)
            max_hard_p_deviation = fabs(hard_p_deviation);
          if (fabs(ks_value_hard - hard_ks_deviation) > max_hard_ks_deviation)
            max_hard_ks_deviation = fabs(hard_ks_deviation);
          if (fabs(p_value_match - match_p_deviation) > max_match_p_deviation)
            max_match_p_deviation = fabs(match_p_deviation);
          if (fabs(ks_value_match - match_ks_deviation) >
              max_match_ks_deviation)
            max_match_ks_deviation = fabs(match_ks_deviation);
        }

        // set bins
        p_values_hard->SetBinContent(x_bin, y_bin, p_value_hard);
        p_values_match->SetBinContent(x_bin, y_bin, p_value_match);
        p_values_bkg->SetBinContent(x_bin, y_bin, p_value_bkg);
        ks_values_hard->SetBinContent(x_bin, y_bin, ks_value_hard);
        ks_values_match->SetBinContent(x_bin, y_bin, ks_value_match);
        ks_values_bkg->SetBinContent(x_bin, y_bin, ks_value_bkg);
        p_values_hard_error->SetBinContent(x_bin, y_bin, max_hard_p_deviation);
        p_values_match_error->SetBinContent(x_bin, y_bin,
                                            max_match_p_deviation);
        ks_values_hard_error->SetBinContent(x_bin, y_bin,
                                            max_hard_ks_deviation);
        ks_values_match_error->SetBinContent(x_bin, y_bin,
                                             max_match_ks_deviation);
      }
    }

    Print2DSimple(p_values_hard, hopts, copts, out_loc_grid, "ajphard", "", "R",
                  "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(p_values_match, hopts, copts, out_loc_grid, "ajpmatch", "",
                  "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(p_values_bkg, hopts, copts, out_loc_grid, "ajpbkg", "", "R",
                  "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_hard, hopts, copts, out_loc_grid, "ajkshard", "",
                  "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_match, hopts, copts, out_loc_grid, "ajksmatch", "",
                  "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_bkg, hopts, copts, out_loc_grid, "ajksbkg", "", "R",
                  "p_{T}^{const}", "TEXT COLZ");

    Print2DSimple(p_values_hard_error, hopts, copts, out_loc_grid,
                  "ajpharderror", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(p_values_match_error, hopts, copts, out_loc_grid,
                  "ajpmatcherror", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_hard_error, hopts, copts, out_loc_grid,
                  "ajksharderror", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_match_error, hopts, copts, out_loc_grid,
                  "ajksmatcherror", "", "R", "p_{T}^{const}", "TEXT COLZ");

    // draw axis labels
    c_hard->cd();
    invis->Draw("clone");
    invis->cd();
    x_label_text->Draw();
    y_label_text->Draw();
    c_hard->SaveAs(
        dijetcore::MakeString(out_loc_grid, "/ajhardgrid.pdf").c_str());
    c_match->cd();
    // draw axis labels
    invis->Draw("clone");
    invis->cd();
    x_label_text->Draw();
    y_label_text->Draw();
    c_match->SaveAs(
        dijetcore::MakeString(out_loc_grid, "/ajmatchgrid.pdf").c_str());
    c_match_oa->cd();
    // draw axis labels
    invis->Draw("clone");
    invis->cd();
    x_label_text->Draw();
    y_label_text->Draw();
    c_match_oa->SaveAs(
        dijetcore::MakeString(out_loc_grid, "/ajmatchgridwithoa.pdf").c_str());
  }  // centrality

  return 0;
}
