#include "dijetcore/util/root/root_print_routines.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/flags.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"

#include "dijetcore/lib/map.h"
#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(auau, "", "input file for Au+Au");
DIJETCORE_DEFINE_string(ppDir, "", "input directory for p+p");
DIJETCORE_DEFINE_string(outputDir, "results", "directory for output");
DIJETCORE_DEFINE_string(radii, "0.2,0.25,0.3,0.35,0.4", "radii to put in grid");
DIJETCORE_DEFINE_string(constPt, "1.0,1.5,2.0,2.5,3.0", "radii to put in grid");
DIJETCORE_DEFINE_bool(useSingleCentrality, true, "do things in 3 centrality bins or 1");

enum class DATATYPE {
  AUAU = 0,
  PP = 1,
  TOWP = 2,
  TOWM = 3,
  TRACKP = 4,
  TRACKM = 5,
  SIZE = 6
};

const std::set<DATATYPE> data_types
{DATATYPE::AUAU, DATATYPE::PP, DATATYPE::TOWP, DATATYPE::TOWM, DATATYPE::TRACKP, DATATYPE::TRACKM};


// adds to the map all TTrees that conform to
// the DijetWorker's naming convention
void GetTreesFromFile(const TFile& file, std::unordered_map<string, TTree*>& map) {
  TKey *key;
  TIter next(file.GetListOfKeys());
  
  while ((key = (TKey*) next())) {
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
    
    map.insert({tmp, (TTree*) key->ReadObj()});
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
    TH2D* tmp = (TH2D*) h->Project3D("zy");
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
    if (tmp->Integral() > 0)
    tmp->Scale(1.0 / tmp->Integral());
    ret.push_back(tmp);
  }
  return ret;
}

std::vector<TH1D*> AddBins(std::vector<TH1D*> container, std::vector<std::pair<int, int>> bins) {
  
  std::vector<TH1D*> ret(bins.size(), nullptr);
  for (int i = 0; i < container.size(); ++i) {
    for (int j = 0; j < bins.size(); ++j) {
      if (i >=bins[j].first && i <= bins[j].second) {
        if (ret[j] == nullptr) {
          string name = string(container[i]->GetName()) + std::to_string(j);
          ret[j] = (TH1D*) container[i]->Clone(name.c_str());
        }
        else {
          ret[j]->Add(container[i]);
        }
      }
    }
  }
  return ret;
}

TGraphErrors* GetSystematic(TH1D* nom, TH1D* var1_a, TH1D* var1_b, TH1D* var2_a, TH1D* var2_b) {
  int nBins = nom->GetNbinsX();
  double x_[nBins];
  double y_[nBins];
  double x_err_[nBins];
  double y_err_[nBins];
  
  for (int i = 0; i < nBins; ++i) {
    x_[i] = nom->GetBinCenter(i+1);
    y_[i] = nom->GetBinContent(i+1);
    x_err_[i] = nom->GetXaxis()->GetBinWidth(1) / 2.0;
    double diff_var_1_a = fabs(nom->GetBinContent(i+1)  - var1_a->GetBinContent(i+1));
    double diff_var_1_b = fabs(nom->GetBinContent(i+1)  - var1_b->GetBinContent(i+1));
    double diff_var_2_a = fabs(nom->GetBinContent(i+1)  - var2_a->GetBinContent(i+1));
    double diff_var_2_b = fabs(nom->GetBinContent(i+1)  - var2_b->GetBinContent(i+1));
    double max_var_1 = (diff_var_1_a > diff_var_1_b ? diff_var_1_a : diff_var_1_b);
    double max_var_2 = (diff_var_2_a > diff_var_2_b ? diff_var_2_a : diff_var_2_b);
    y_err_[i] = sqrt(max_var_1 * max_var_1 + max_var_2 * max_var_2);
  }
  TGraphErrors* ret = new TGraphErrors(nBins, x_, y_, x_err_, y_err_);
  return ret;
}

template<class T>
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
    double bin_error = errors->GetErrorY(i-1);
    ret->SetBinContent(i, bin_error / bin_content);
  }
  return ret;
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
  gErrorIgnoreLevel = kInfo+1;
  
  // check to make sure we have valid inputs
  std::vector<string> inputs{FLAGS_auau, FLAGS_ppDir + "/nom.root"};
  for (auto& file : inputs) {
    if (!boost::filesystem::exists(file)) {
      std::cout << "input file " << file;
      std::cout << "doesn't exist: exiting" << std::endl;
      return 1;
    }
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);
  
  // read in the file
  TFile auau_file(FLAGS_auau.c_str(), "READ");
  TFile pp_file((FLAGS_ppDir + "/nom.root").c_str(), "READ");
  TFile tow_p_file((FLAGS_ppDir + "/tow_p.root").c_str(), "READ");
  TFile tow_m_file((FLAGS_ppDir + "/tow_m.root").c_str(), "READ");
  TFile track_p_file((FLAGS_ppDir + "/track_p.root").c_str(), "READ");
  TFile track_m_file((FLAGS_ppDir + "/track_m.root").c_str(), "READ");

  // define centralities
  std::vector<unsigned> cent_boundaries;
  std::vector<std::pair<int, int>> cent_bin_boundaries;
  std::vector<string> refcent_string;
  
  if (FLAGS_useSingleCentrality) {
    cent_boundaries = {269};
    cent_bin_boundaries = {{0, 2}};
    refcent_string = {"0-20%"};
  }
  else {
    cent_boundaries = {485, 399, 269};
    cent_bin_boundaries = {{0, 0}, {1, 1}, {2, 2}};
    refcent_string = {"0-5%", "5-10%", "10-20%"};
  }
  
  // and the radii we want to use
  std::vector<double> radii = dijetcore::ParseArgStringToVec<double>(FLAGS_radii);
  std::sort(radii.begin(), radii.end());
  std::vector<string> radii_string;
  for (auto& val : radii)
    radii_string.push_back(std::to_string(val));
  std::vector<double> constpt = dijetcore::ParseArgStringToVec<double>(FLAGS_constPt);
  std::sort(constpt.begin(), constpt.end());
  std::vector<string> constpt_string;
  for (auto& val : constpt)
    constpt_string.push_back(std::to_string(val));
  LOG(INFO) << "R x const Pt matrix selected";
  LOG(INFO) << "R: " << radii_string;
  LOG(INFO) << "constPt: " << constpt_string;
  
  // now we'll get the trees from the files, ignoring any objects
  // in the file that don't conform to the naming conventions from
  // the DijetWorker. There are also coincidence histograms to save
  std::vector<string> keys;
  std::unordered_map<std::string, dijetcore::DijetKey> parsed_keys;
  std::unordered_map<std::string, TTree*> auau_trees;
  std::unordered_map<std::string, TTree*> pp_trees;
  std::unordered_map<std::string, TTree*> tow_p_trees;
  std::unordered_map<std::string, TTree*> tow_m_trees;
  std::unordered_map<std::string, TTree*> track_p_trees;
  std::unordered_map<std::string, TTree*> track_m_trees;
  
  
  GetTreesFromFile(auau_file, auau_trees);
  GetTreesFromFile(pp_file, pp_trees);
  GetTreesFromFile(tow_p_file, tow_p_trees);
  GetTreesFromFile(tow_m_file, tow_m_trees);
  GetTreesFromFile(track_p_file, track_p_trees);
  GetTreesFromFile(track_m_file, track_m_trees);
  
  // match keys
  for (auto entry : auau_trees) {
    if (pp_trees.find(entry.first) != pp_trees.end()) {
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
      for (int key_index = 0; key_index < keys.size(); ++key_index) {
        string key_string = keys[key_index];
        dijetcore::DijetKey& key = parsed_keys[key_string];
        if (key.lead_init_r == radii[rad] &&
            key.sub_init_r == radii[rad] &&
            key.lead_match_r == radii[rad] &&
            key.sub_match_r == radii[rad] &&
            key.lead_init_const_pt == constpt[pt] &&
            key.sub_init_const_pt == constpt[pt]) {
          grid_keys[rad].push_back(keys[key_index]);
          grid_key_params[rad].push_back(key);
        }
      }
      if (grid_keys[rad].size() < pt + 1) {
        LOG(INFO) << "could not find a key for R=" << radii[rad] << " const Pt=" << constpt[pt];
        return 1;
      }
      if (grid_keys[rad].size() > pt + 1) {
        LOG(INFO) << "found multiple keys for R=" << radii[rad] << " const Pt=" << constpt[pt];
        return 1;
      }
    }
  }

  // save all histograms so we can do comparisons
  // between different keys if we want
  std::unordered_map<string, TH2D*> auau_hard_lead_pt;
  std::unordered_map<string, TH2D*> auau_hard_sub_pt;
  std::unordered_map<string, TH2D*> auau_match_lead_pt;
  std::unordered_map<string, TH2D*> auau_match_sub_pt;
  std::unordered_map<string, TH2D*> auau_hard_aj;
  std::unordered_map<string, TH2D*> auau_match_aj;
  std::unordered_map<string, TH2D*> auau_hard_aj_test;
  std::unordered_map<string, TH2D*> auau_match_aj_test;
  std::unordered_map<string, TH2D*> auau_dphi;
  std::unordered_map<string, TH2D*> auau_hard_lead_rp;
  std::unordered_map<string, TH2D*> auau_off_axis_lead_pt;
  std::unordered_map<string, TH2D*> auau_off_axis_sub_pt;
  std::unordered_map<string, TH2D*> auau_off_axis_aj;
  
  std::unordered_map<string, TH2D*> pp_hard_lead_pt;
  std::unordered_map<string, TH2D*> pp_hard_sub_pt;
  std::unordered_map<string, TH2D*> pp_match_lead_pt;
  std::unordered_map<string, TH2D*> pp_match_sub_pt;
  std::unordered_map<string, TH2D*> pp_hard_aj;
  std::unordered_map<string, TH2D*> pp_match_aj;
  std::unordered_map<string, TH2D*> pp_hard_aj_test;
  std::unordered_map<string, TH2D*> pp_match_aj_test;
  std::unordered_map<string, TH2D*> pp_dphi;
  std::unordered_map<string, TH2D*> pp_hard_lead_rp;
  
  std::unordered_map<string, TH2D*> pp_only_hard_lead_pt;
  std::unordered_map<string, TH2D*> pp_only_hard_sub_pt;
  std::unordered_map<string, TH2D*> pp_only_match_lead_pt;
  std::unordered_map<string, TH2D*> pp_only_match_sub_pt;
  std::unordered_map<string, TH2D*> pp_only_hard_aj;
  std::unordered_map<string, TH2D*> pp_only_match_aj;
  
  std::unordered_map<string, TH2D*> pp_hard_lead_dpt;
  std::unordered_map<string, TH2D*> pp_hard_sub_dpt;
  std::unordered_map<string, TH2D*> pp_match_lead_dpt;
  std::unordered_map<string, TH2D*> pp_match_sub_dpt;
  
  std::unordered_map<string, TH2D*> pp_lead_hard_match_fraction;
  std::unordered_map<string, TH2D*> pp_sub_hard_match_fraction;
  
  std::unordered_map<string, TH2D*> pp_hard_lead_dpt_frac;
  std::unordered_map<string, TH2D*> pp_hard_sub_dpt_frac;
  std::unordered_map<string, TH2D*> pp_match_lead_dpt_frac;
  std::unordered_map<string, TH2D*> pp_match_sub_dpt_frac;
  
  std::unordered_map<string, TH3D*> pp_hard_lead_pp_3d;
  std::unordered_map<string, TH3D*> pp_hard_sub_pp_3d;
  std::unordered_map<string, TH3D*> pp_match_lead_pp_3d;
  std::unordered_map<string, TH3D*> pp_match_sub_pp_3d;
  
  std::unordered_map<string, TH2D*> auau_hard_lead_const;
  std::unordered_map<string, TH2D*> auau_hard_sub_const;
  std::unordered_map<string, TH2D*> auau_match_lead_const;
  std::unordered_map<string, TH2D*> auau_match_sub_const;
  std::unordered_map<string, TH2D*> auau_off_axis_lead_const;
  std::unordered_map<string, TH2D*> auau_off_axis_sub_const;
  
  std::unordered_map<string, TH2D*> pp_hard_lead_const;
  std::unordered_map<string, TH2D*> pp_hard_sub_const;
  std::unordered_map<string, TH2D*> pp_match_lead_const;
  std::unordered_map<string, TH2D*> pp_match_sub_const;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_aj_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_aj_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_aj_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_test_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_only_hard_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_only_hard_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_only_match_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_only_match_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_only_hard_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_only_match_aj_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_dpt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_dpt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_dpt_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_dpt_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_lead_hard_match_fraction_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_sub_hard_match_fraction_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_dpt_frac_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_dpt_frac_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_dpt_frac_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_dpt_frac_cent;
  
  std::unordered_map<string, std::vector<TH2D*>> pp_hard_lead_pp_3d_cent;
  std::unordered_map<string, std::vector<TH2D*>> pp_hard_sub_pp_3d_cent;
  std::unordered_map<string, std::vector<TH2D*>> pp_match_lead_pp_3d_cent;
  std::unordered_map<string, std::vector<TH2D*>> pp_match_sub_pp_3d_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_lead_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_sub_const_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_const_cent;
  
  std::unordered_map<string, TH2D*> auau_hard_lead_rho;
  std::unordered_map<string, TH2D*> auau_hard_sub_rho;
  std::unordered_map<string, TH2D*> auau_match_lead_rho;
  std::unordered_map<string, TH2D*> auau_match_sub_rho;
  
  std::unordered_map<string, TH2D*> pp_hard_lead_rho;
  std::unordered_map<string, TH2D*> pp_hard_sub_rho;
  std::unordered_map<string, TH2D*> pp_match_lead_rho;
  std::unordered_map<string, TH2D*> pp_match_sub_rho;
  
  std::unordered_map<string, TH2D*> auau_hard_lead_sig;
  std::unordered_map<string, TH2D*> auau_hard_sub_sig;
  std::unordered_map<string, TH2D*> auau_match_lead_sig;
  std::unordered_map<string, TH2D*> auau_match_sub_sig;
  
  std::unordered_map<string, TH2D*> pp_hard_lead_sig;
  std::unordered_map<string, TH2D*> pp_hard_sub_sig;
  std::unordered_map<string, TH2D*> pp_match_lead_sig;
  std::unordered_map<string, TH2D*> pp_match_sub_sig;
  
  std::unordered_map<string, TH2D*> auau_off_axis_lead_rho;
  std::unordered_map<string, TH2D*> auau_off_axis_sub_rho;
  std::unordered_map<string, TH2D*> auau_off_axis_lead_sig;
  std::unordered_map<string, TH2D*> auau_off_axis_sub_sig;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_rho_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_rho_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_sig_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_sub_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_sub_sig_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_sub_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_sub_sig_cent;
  
  std::unordered_map<string, TH2D*> auau_npart;
  std::unordered_map<string, TH2D*> pp_npart;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_npart_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_npart_cent;
  
  // for systematics
  std::unordered_map<string, TH2D*> pp_hard_aj_tow_p;
  std::unordered_map<string, TH2D*> pp_match_aj_tow_p;
  std::unordered_map<string, TH2D*> pp_hard_aj_tow_m;
  std::unordered_map<string, TH2D*> pp_match_aj_tow_m;
  std::unordered_map<string, TH2D*> pp_hard_aj_track_p;
  std::unordered_map<string, TH2D*> pp_match_aj_track_p;
  std::unordered_map<string, TH2D*> pp_hard_aj_track_m;
  std::unordered_map<string, TH2D*> pp_match_aj_track_m;
  
  std::unordered_map<string, TH2D*> pp_hard_aj_tow_p_test;
  std::unordered_map<string, TH2D*> pp_match_aj_tow_p_test;
  std::unordered_map<string, TH2D*> pp_hard_aj_tow_m_test;
  std::unordered_map<string, TH2D*> pp_match_aj_tow_m_test;
  std::unordered_map<string, TH2D*> pp_hard_aj_track_p_test;
  std::unordered_map<string, TH2D*> pp_match_aj_track_p_test;
  std::unordered_map<string, TH2D*> pp_hard_aj_track_m_test;
  std::unordered_map<string, TH2D*> pp_match_aj_track_m_test;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_tow_p_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_tow_p_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_tow_m_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_tow_m_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_track_p_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_track_p_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_track_m_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_track_m_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_tow_p_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_tow_p_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_tow_m_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_tow_m_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_track_p_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_track_p_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_hard_aj_track_m_test_cent;
  std::unordered_map<string, std::vector<TH1D*>> pp_match_aj_track_m_test_cent;
  
  std::unordered_map<string, std::vector<TGraphErrors*>> systematic_errors_hard;
  std::unordered_map<string, std::vector<TGraphErrors*>> systematic_errors_match;
  
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
  
  // count which key we are on
  int entry = -1;
  // loop over all matched trees
  for (auto key : keys) {
    
    //count which key we're on
    entry++;
    
    // make output directory
    string out_loc = FLAGS_outputDir + "/" + key;
    boost::filesystem::path dir(out_loc.c_str());
    boost::filesystem::create_directories(dir);
    
    // key prefix for names so histograms don't get confused
    // (root fuckin sucks)
    std::string key_prefix = "key_" + std::to_string(entry) + "_";
    
    // and create the file name prefix
    string file_prefix = out_loc + "/";
    
    // process the trees
    TTree* auau_tree = auau_trees[key];
    TTreeReader auau_reader(auau_tree);
    
    TTree* pp_tree = pp_trees[key];
    TTreeReader pp_reader(pp_tree);
    
    TTree* tow_p_tree = tow_p_trees[key];
    TTree* tow_m_tree = tow_m_trees[key];
    TTree* track_p_tree = track_p_trees[key];
    TTree* track_m_tree = track_m_trees[key];
    TTreeReader tow_p_reader(tow_p_tree);
    TTreeReader tow_m_reader(tow_m_tree);
    TTreeReader track_p_reader(track_p_tree);
    TTreeReader track_m_reader(track_m_tree);
    
    // create readervalues for auau first
    TTreeReaderValue<int> auau_runid(auau_reader, "runid");
    TTreeReaderValue<int> auau_eventid(auau_reader, "eventid");
    TTreeReaderValue<double> auau_vz(auau_reader, "vz");
    TTreeReaderValue<int> auau_refmult(auau_reader, "refmult");
    TTreeReaderValue<int> auau_grefmult(auau_reader, "grefmult");
    TTreeReaderValue<double> auau_refmultcorr(auau_reader, "refmultcorr");
    TTreeReaderValue<double> auau_grefmultcorr(auau_reader, "grefmultcorr");
    TTreeReaderValue<int> auau_cent(auau_reader, "cent");
    TTreeReaderValue<double> auau_zdcrate(auau_reader, "zdcrate");
    TTreeReaderValue<double> auau_rp(auau_reader, "rp");
    TTreeReaderValue<int> auau_nglobal(auau_reader, "nglobal");
    TTreeReaderValue<int> auau_nprt(auau_reader, "npart");
    TTreeReaderValue<TLorentzVector> auau_jl(auau_reader, "jl");
    TTreeReaderValue<TLorentzVector> auau_js(auau_reader, "js");
    TTreeReaderValue<TLorentzVector> auau_jlm(auau_reader, "jlm");
    TTreeReaderValue<TLorentzVector> auau_jsm(auau_reader, "jsm");
    TTreeReaderValue<TLorentzVector> auau_jloa(auau_reader, "jloa");
    TTreeReaderValue<TLorentzVector> auau_jsoa(auau_reader, "jsoa");
    TTreeReaderValue<int> auau_jlconst(auau_reader, "jlconst");
    TTreeReaderValue<double> auau_jlrho(auau_reader, "jlrho");
    TTreeReaderValue<double> auau_jlsig(auau_reader, "jlsig");
    TTreeReaderValue<int> auau_jlmconst(auau_reader, "jlmconst");
    TTreeReaderValue<double> auau_jlmrho(auau_reader, "jlmrho");
    TTreeReaderValue<double> auau_jlmsig(auau_reader, "jlmsig");
    TTreeReaderValue<int> auau_jsconst(auau_reader, "jsconst");
    TTreeReaderValue<double> auau_jsrho(auau_reader, "jsrho");
    TTreeReaderValue<double> auau_jssig(auau_reader, "jssig");
    TTreeReaderValue<int> auau_jsmconst(auau_reader, "jsmconst");
    TTreeReaderValue<double> auau_jsmrho(auau_reader, "jsmrho");
    TTreeReaderValue<double> auau_jsmsig(auau_reader, "jsmsig");
    TTreeReaderValue<int> auau_jloaconst(auau_reader, "jloaconst");
    TTreeReaderValue<double> auau_jloarho(auau_reader, "jloarho");
    TTreeReaderValue<double> auau_jloasig(auau_reader, "jloasig");
    TTreeReaderValue<int> auau_jsoaconst(auau_reader, "jsoaconst");
    TTreeReaderValue<double> auau_jsoarho(auau_reader, "jsoarho");
    TTreeReaderValue<double> auau_jsoasig(auau_reader, "jsoasig");
    
    // create readervalues for pp
    TTreeReaderValue<int> pp_runid(pp_reader, "runid");
    TTreeReaderValue<int> pp_eventid(pp_reader, "eventid");
    TTreeReaderValue<double> pp_vz(pp_reader, "vz");
    TTreeReaderValue<int> pp_refmult(pp_reader, "refmult");
    TTreeReaderValue<int> pp_grefmult(pp_reader, "grefmult");
    TTreeReaderValue<double> pp_refmultcorr(pp_reader, "refmultcorr");
    TTreeReaderValue<double> pp_grefmultcorr(pp_reader, "grefmultcorr");
    TTreeReaderValue<int> pp_cent(pp_reader, "cent");
    TTreeReaderValue<double> pp_zdcrate(pp_reader, "zdcrate");
    TTreeReaderValue<double> pp_rp(pp_reader, "rp");
    TTreeReaderValue<int> pp_nglobal(pp_reader, "nglobal");
    TTreeReaderValue<int> pp_nprt(pp_reader, "npart");
    TTreeReaderValue<TLorentzVector> pp_jl(pp_reader, "jl");
    TTreeReaderValue<TLorentzVector> pp_js(pp_reader, "js");
    TTreeReaderValue<TLorentzVector> pp_jlm(pp_reader, "jlm");
    TTreeReaderValue<TLorentzVector> pp_jsm(pp_reader, "jsm");
    TTreeReaderValue<bool> pp_only_found_match(pp_reader, "foundpp");
    TTreeReaderValue<TLorentzVector> pp_only_jl(pp_reader, "ppjl");
    TTreeReaderValue<TLorentzVector> pp_only_js(pp_reader, "ppjs");
    TTreeReaderValue<TLorentzVector> pp_only_jlm(pp_reader, "ppjlm");
    TTreeReaderValue<TLorentzVector> pp_only_jsm(pp_reader, "ppjsm");
    TTreeReaderValue<int> pp_jlconst(pp_reader, "jlconst");
    TTreeReaderValue<double> pp_jlrho(pp_reader, "jlrho");
    TTreeReaderValue<double> pp_jlsig(pp_reader, "jlsig");
    TTreeReaderValue<int> pp_jlmconst(pp_reader, "jlmconst");
    TTreeReaderValue<double> pp_jlmrho(pp_reader, "jlmrho");
    TTreeReaderValue<double> pp_jlmsig(pp_reader, "jlmsig");
    TTreeReaderValue<int> pp_jsconst(pp_reader, "jsconst");
    TTreeReaderValue<double> pp_jsrho(pp_reader, "jsrho");
    TTreeReaderValue<double> pp_jssig(pp_reader, "jssig");
    TTreeReaderValue<int> pp_jsmconst(pp_reader, "jsmconst");
    TTreeReaderValue<double> pp_jsmrho(pp_reader, "jsmrho");
    TTreeReaderValue<double> pp_jsmsig(pp_reader, "jsmsig");
    
    // create the readers for the 4 systematic variations
    TTreeReaderValue<TLorentzVector> tow_p_jl(tow_p_reader, "jl");
    TTreeReaderValue<TLorentzVector> tow_p_js(tow_p_reader, "js");
    TTreeReaderValue<TLorentzVector> tow_p_jlm(tow_p_reader, "jlm");
    TTreeReaderValue<TLorentzVector> tow_p_jsm(tow_p_reader, "jsm");
    TTreeReaderValue<int> tow_p_cent(tow_p_reader, "cent");
    
    TTreeReaderValue<TLorentzVector> tow_m_jl(tow_m_reader, "jl");
    TTreeReaderValue<TLorentzVector> tow_m_js(tow_m_reader, "js");
    TTreeReaderValue<TLorentzVector> tow_m_jlm(tow_m_reader, "jlm");
    TTreeReaderValue<TLorentzVector> tow_m_jsm(tow_m_reader, "jsm");
    TTreeReaderValue<int> tow_m_cent(tow_m_reader, "cent");
    
    TTreeReaderValue<TLorentzVector> track_p_jl(track_p_reader, "jl");
    TTreeReaderValue<TLorentzVector> track_p_js(track_p_reader, "js");
    TTreeReaderValue<TLorentzVector> track_p_jlm(track_p_reader, "jlm");
    TTreeReaderValue<TLorentzVector> track_p_jsm(track_p_reader, "jsm");
    TTreeReaderValue<int> track_p_cent(track_p_reader, "cent");
    
    TTreeReaderValue<TLorentzVector> track_m_jl(track_m_reader, "jl");
    TTreeReaderValue<TLorentzVector> track_m_js(track_m_reader, "js");
    TTreeReaderValue<TLorentzVector> track_m_jlm(track_m_reader, "jlm");
    TTreeReaderValue<TLorentzVector> track_m_jsm(track_m_reader, "jsm");
    TTreeReaderValue<int> track_m_cent(track_m_reader, "cent");
    
    // and check if we have embedding in the pp tree
    bool pp_embedded = false;
    if (pp_tree->GetBranch("embed_eventid"))
      pp_embedded = true;
    // build all embedding branches (only used if pp_embedded is true)
    TTreeReaderValue<int> embed_eventid(pp_reader, "embed_eventid");
    TTreeReaderValue<int> embed_runid(pp_reader, "embed_runid");
    TTreeReaderValue<int> embed_refmult(pp_reader, "embed_refmult");
    TTreeReaderValue<int> embed_grefmult(pp_reader, "embed_grefmult");
    TTreeReaderValue<int> embed_nprt(pp_reader, "embed_npart");
    TTreeReaderValue<double> embed_refmultcorr(pp_reader, "embed_refmultcorr");
    TTreeReaderValue<double> embed_grefmultcorr(pp_reader, "embed_grefmultcorr");
    TTreeReaderValue<int> embed_cent(pp_reader, "embed_cent");
    TTreeReaderValue<double> embed_rp(pp_reader, "embed_rp");
    TTreeReaderValue<double> embed_zdcrate(pp_reader, "embed_zdcrate");
    TTreeReaderValue<double> embed_vz(pp_reader, "embed_vz");
    
    // insert into the dictionaries
    auau_hard_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardleadpt").c_str(),
                                      "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_hard_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardsubpt").c_str(),
                                     "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_match_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchleadpt").c_str(),
                                       "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_match_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchsubpt").c_str(),
                                      "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_hard_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardaj").c_str(),
                                 "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    auau_match_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchaj").c_str(),
                                  "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    auau_hard_aj_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardajtest").c_str(),
                                 "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    auau_match_aj_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchajtest").c_str(),
                                  "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    auau_dphi[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaudphi").c_str(),
                              "d#phi", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 2*TMath::Pi());
    
    pp_hard_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardleadpt").c_str(),
                                      "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_hard_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardsubpt").c_str(),
                                     "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_match_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchleadpt").c_str(),
                                       "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_match_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchsubpt").c_str(),
                                      "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_hard_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardaj").c_str(),
                                 "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_match_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchaj").c_str(),
                                  "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_hard_aj_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardajtest").c_str(),
                               "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_match_aj_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchajtest").c_str(),
                                "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_dphi[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdphi").c_str(),
                              "d#phi", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 2*TMath::Pi());
    
    pp_only_hard_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "pponlyhardleadpt").c_str(), "",
                                    cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_only_hard_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "pponlyhardsubpt").c_str(), "",
                                   cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_only_match_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "pponlymatchleadpt").c_str(), "",
                                     cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_only_match_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "pponlymatchsubpt").c_str(), "",
                                    cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_only_hard_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "pponlyhardaj").c_str(),
                               "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_only_match_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "pponlymatchaj").c_str(),
                                "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    
    pp_hard_lead_dpt[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdptleadhard").c_str(),
                                   "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 30);
    pp_hard_sub_dpt[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdptsubhard").c_str(),
                               "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 30);
    pp_match_lead_dpt[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdptleadmatch").c_str(),
                                 "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 30);
    pp_match_sub_dpt[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdptsubmatch").c_str(),
                                "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 30);
    
    pp_hard_lead_pp_3d[key] = new TH3D(dijetcore::MakeString(key_prefix, "ppleadhardpt3d").c_str(), "",
                                       cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                                       50, 0, 50, 50, 0, 50);
    pp_hard_sub_pp_3d[key] = new TH3D(dijetcore::MakeString(key_prefix, "ppsubhardpt3d").c_str(), "",
                                      cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                                      50, 0, 50, 50, 0, 50);
    pp_match_lead_pp_3d[key] = new TH3D(dijetcore::MakeString(key_prefix, "ppleadmatchpt3d").c_str(), "",
                                        cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                                        50, 0, 50, 50, 0, 50);
    pp_match_sub_pp_3d[key] = new TH3D(dijetcore::MakeString(key_prefix, "ppsubmatchpt3d").c_str(), "",
                                       cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                                       50, 0, 50, 50, 0, 50);
    
    pp_lead_hard_match_fraction[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppleadhardmatchfrac").c_str(),
                                                "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                                                2, -0.5, 1.5);
    pp_sub_hard_match_fraction[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppsubhardmatchfrac").c_str(),
                                               "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5,
                                               2, -0.5, 1.5);
    
    pp_hard_lead_dpt_frac[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdptleadhardfrac").c_str(),
                                     "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 20, 0, 1.01);
    pp_hard_sub_dpt_frac[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdptsubhardfrac").c_str(),
                                    "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 20, 0, 1.01);
    pp_match_lead_dpt_frac[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdptleadmatchfrac").c_str(),
                                      "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 20, 0, 1.01);
    pp_match_sub_dpt_frac[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppdptsubmatchfrac").c_str(),
                                     "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 20, 0, 1.01);
    
    
    
    auau_off_axis_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisleadpt").c_str(),
                                          "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_off_axis_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxissubpt").c_str(),
                                         "p_{T}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_off_axis_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisaj").c_str(),
                                     "A_{J}", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    
    auau_hard_lead_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardleadconst").c_str(),
                                         "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_hard_sub_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardsubconst").c_str(),
                                        "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_match_lead_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchleadconst").c_str(),
                                          "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_match_sub_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchsubconst").c_str(),
                                         "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_off_axis_lead_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisleadconst").c_str(),
                                             "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_off_axis_sub_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxissubconst").c_str(),
                                            "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    
    pp_hard_lead_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardleadconst").c_str(),
                                         "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_hard_sub_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardsubconst").c_str(),
                                        "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_match_lead_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchleadconst").c_str(),
                                          "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_match_sub_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchsubconst").c_str(),
                                         "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    
    auau_hard_lead_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardleadrho").c_str(),
                                       "auau hard lead rho", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_hard_lead_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardleadsig").c_str(),
                                       "auau hard lead sig", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    auau_hard_sub_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardsubrho").c_str(),
                                      "auau hard sub rho", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_hard_sub_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardsubsig").c_str(),
                                      "auau hard sub sig", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    auau_match_lead_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchleadrho").c_str(),
                                        "auau match lead rho", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_match_lead_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchleadsig").c_str(),
                                        "auau match lead sig", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    auau_match_sub_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchsubrho").c_str(),
                                       "auau match sub rho", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_match_sub_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchsubsig").c_str(),
                                       "auau match sub sig", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    
    pp_hard_lead_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardleadrho").c_str(),
                                       "pp hard lead rho", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_hard_lead_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardleadsig").c_str(),
                                       "pp hard lead sig", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    pp_hard_sub_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardsubrho").c_str(),
                                      "pp hard sub rho", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_hard_sub_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "pphardsubsig").c_str(),
                                      "pp hard sub sig", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    pp_match_lead_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchleadrho").c_str(),
                                        "pp match lead rho", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_match_lead_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchleadsig").c_str(),
                                        "pp match lead sig", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    pp_match_sub_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchsubrho").c_str(),
                                       "pp match sub rho", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    pp_match_sub_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppmatchsubsig").c_str(),
                                       "pp match sub sig", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    
    auau_off_axis_lead_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisleadrho").c_str(),
                                           "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_off_axis_sub_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxissubrho").c_str(),
                                          "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 100);
    auau_off_axis_lead_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisleadsig").c_str(),
                                           "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    auau_off_axis_sub_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxissubsig").c_str(),
                                          "", cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0, 20);
    
    pp_hard_aj_tow_p[key] = new TH2D(dijetcore::MakeString(key_prefix, "towpajhard").c_str(), "",
                                     cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_match_aj_tow_p[key] = new TH2D(dijetcore::MakeString(key_prefix, "towpajmatch").c_str(), "",
                                      cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_hard_aj_tow_m[key] = new TH2D(dijetcore::MakeString(key_prefix, "towmajhard").c_str(), "",
                                     cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_match_aj_tow_m[key] = new TH2D(dijetcore::MakeString(key_prefix, "towmajmatch").c_str(), "",
                                      cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_hard_aj_track_p[key] = new TH2D(dijetcore::MakeString(key_prefix, "trackpajhard").c_str(), "",
                                       cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_match_aj_track_p[key] = new TH2D(dijetcore::MakeString(key_prefix, "trackpajmatch").c_str(), "",
                                        cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_hard_aj_track_m[key] = new TH2D(dijetcore::MakeString(key_prefix, "trackmajhard").c_str(), "",
                                       cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    pp_match_aj_track_m[key] = new TH2D(dijetcore::MakeString(key_prefix, "trackmajmatch").c_str(), "",
                                        cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 30, 0, 0.9);
    
    pp_hard_aj_tow_p_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "towpajhardtest").c_str(), "",
                                          cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_match_aj_tow_p_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "towpajmatchtest").c_str(), "",
                                           cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_hard_aj_tow_m_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "towmajhardtest").c_str(), "",
                                          cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_match_aj_tow_m_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "towmajmatchtest").c_str(), "",
                                           cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_hard_aj_track_p_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "trackpajhardtest").c_str(), "",
                                            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_match_aj_track_p_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "trackpajmatchtest").c_str(), "",
                                             cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_hard_aj_track_m_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "trackmajhardtest").c_str(), "",
                                            cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    pp_match_aj_track_m_test[key] = new TH2D(dijetcore::MakeString(key_prefix, "trackmajmatchtest").c_str(), "",
                                             cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 10000, 0, 0.9);
    
    auau_npart[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaunpart").c_str(), ";refmult;nPart",
                               cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0.5, 2500.5);
    pp_npart[key] = new TH2D(dijetcore::MakeString(key_prefix, "ppnpart").c_str(), ";refmult;nPart",
                               cent_boundaries.size(), -0.5, cent_boundaries.size() - 0.5, 100, 0.5, 2500.5);
    
    
    
    // loop over the data & fill histograms
    while (auau_reader.Next()) {
      
      int cent_bin = -1;
      for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
        auto& pair = cent_bin_boundaries[i];
        if (*auau_cent >= pair.first &&
            *auau_cent <= pair.second) {
          cent_bin = i;
          break;
        }
      }
      
      if (cent_bin == -1)
        continue;
      
      // auau jet pt
      auau_hard_lead_pt[key]->Fill(cent_bin, (*auau_jl).Pt());
      auau_hard_lead_rho[key]->Fill(cent_bin, *auau_jlrho);
      auau_hard_lead_sig[key]->Fill(cent_bin, *auau_jlsig);
      auau_hard_sub_pt[key]->Fill(cent_bin, (*auau_js).Pt());
      auau_hard_sub_rho[key]->Fill(cent_bin, *auau_jsrho);
      auau_hard_sub_sig[key]->Fill(cent_bin, *auau_jssig);
      auau_match_lead_pt[key]->Fill(cent_bin,(*auau_jlm).Pt());
      auau_match_lead_rho[key]->Fill(cent_bin, *auau_jlmrho);
      auau_match_lead_sig[key]->Fill(cent_bin, *auau_jlmsig);
      auau_match_sub_pt[key]->Fill(cent_bin, (*auau_jsm).Pt());
      auau_match_sub_rho[key]->Fill(cent_bin, *auau_jsmrho);
      auau_match_sub_sig[key]->Fill(cent_bin, *auau_jsmsig);
      auau_hard_lead_const[key]->Fill(cent_bin, *auau_jlconst);
      auau_hard_sub_const[key]->Fill(cent_bin, *auau_jsconst);
      auau_match_lead_const[key]->Fill(cent_bin, *auau_jlmconst);
      auau_match_sub_const[key]->Fill(cent_bin, *auau_jsmconst);
      
      if (*auau_jloaconst != 0 && *auau_jsoaconst != 0) {
        auau_off_axis_lead_pt[key]->Fill(cent_bin, (*auau_jloa).Pt());
        auau_off_axis_lead_rho[key]->Fill(cent_bin, *auau_jloarho);
        auau_off_axis_lead_sig[key]->Fill(cent_bin, *auau_jloasig);
        auau_off_axis_sub_pt[key]->Fill(cent_bin, (*auau_jsoa).Pt());
        auau_off_axis_sub_rho[key]->Fill(cent_bin, *auau_jsoarho);
        auau_off_axis_sub_sig[key]->Fill(cent_bin, *auau_jsoasig);
        auau_off_axis_lead_const[key]->Fill(cent_bin, *auau_jloaconst);
        auau_off_axis_sub_const[key]->Fill(cent_bin, *auau_jsoaconst);
        auau_off_axis_aj[key]->Fill(cent_bin,
                                    fabs((*auau_jloa).Pt() - (*auau_jsoa).Pt())/((*auau_jloa).Pt() + (*auau_jsoa).Pt()));
      }
      // auau Aj
      auau_hard_aj[key]->Fill(cent_bin,
                              fabs((*auau_jl).Pt() - (*auau_js).Pt())/((*auau_jl).Pt() + (*auau_js).Pt()));
      auau_match_aj[key]->Fill(cent_bin,
                               fabs((*auau_jlm).Pt() - (*auau_jsm).Pt())/((*auau_jlm).Pt() + (*auau_jsm).Pt()));
      auau_hard_aj_test[key]->Fill(cent_bin,
                                   fabs((*auau_jl).Pt() - (*auau_js).Pt())/((*auau_jl).Pt() + (*auau_js).Pt()));
      auau_match_aj_test[key]->Fill(cent_bin,
                                    fabs((*auau_jlm).Pt() - (*auau_jsm).Pt())/((*auau_jlm).Pt() + (*auau_jsm).Pt()));
      
      // auau dphi
      double dphi = (*auau_jl).Phi() - (*auau_js).Phi();
      
      //rotate dphi to be in [0,2*pi]
      while (dphi < 0)
      dphi += 2 * TMath::Pi();
      while (dphi > 2.0 * TMath::Pi())
      dphi -= 2.0 * TMath::Pi();
      
      auau_dphi[key]->Fill(cent_bin, dphi);
      
      // and grefmult/npart
      auau_npart[key]->Fill(cent_bin, *auau_nprt);
    }
    
    // loop over the data & fill histograms
    while (pp_reader.Next()) {
      
      int cent_bin = -1;
      for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
        auto& pair = cent_bin_boundaries[i];
        if (*pp_cent >= pair.first &&
            *pp_cent <= pair.second) {
          cent_bin = i;
          break;
        }
      }
      
      if (cent_bin == -1)
        continue;
      
      // pp jet pt
      pp_hard_lead_pt[key]->Fill(cent_bin, (*pp_jl).Pt());
      pp_hard_lead_rho[key]->Fill(cent_bin, *pp_jlrho);
      pp_hard_lead_sig[key]->Fill(cent_bin, *pp_jlsig);
      pp_hard_sub_pt[key]->Fill(cent_bin, (*pp_js).Pt());
      pp_hard_sub_rho[key]->Fill(cent_bin, *pp_jsrho);
      pp_hard_sub_sig[key]->Fill(cent_bin, *pp_jssig);
      pp_match_lead_pt[key]->Fill(cent_bin,(*pp_jlm).Pt());
      pp_match_lead_rho[key]->Fill(cent_bin, *pp_jlmrho);
      pp_match_lead_sig[key]->Fill(cent_bin, *pp_jlmsig);
      pp_match_sub_pt[key]->Fill(cent_bin, (*pp_jsm).Pt());
      pp_match_sub_rho[key]->Fill(cent_bin, *pp_jsmrho);
      pp_match_sub_sig[key]->Fill(cent_bin, *pp_jsmsig);
      pp_hard_lead_const[key]->Fill(cent_bin, *pp_jlconst);
      pp_hard_sub_const[key]->Fill(cent_bin, *pp_jsconst);
      pp_match_lead_const[key]->Fill(cent_bin, *pp_jlmconst);
      pp_match_sub_const[key]->Fill(cent_bin, *pp_jsmconst);
      
      
      // pp Aj
      pp_hard_aj[key]->Fill(cent_bin,
                              fabs((*pp_jl).Pt() - (*pp_js).Pt())/((*pp_jl).Pt() + (*pp_js).Pt()));
      pp_match_aj[key]->Fill(cent_bin,
                               fabs((*pp_jlm).Pt() - (*pp_jsm).Pt())/((*pp_jlm).Pt() + (*pp_jsm).Pt()));
      pp_hard_aj_test[key]->Fill(cent_bin,
                                 fabs((*pp_jl).Pt() - (*pp_js).Pt())/((*pp_jl).Pt() + (*pp_js).Pt()));
      pp_match_aj_test[key]->Fill(cent_bin,
                                  fabs((*pp_jlm).Pt() - (*pp_jsm).Pt())/((*pp_jlm).Pt() + (*pp_jsm).Pt()));
      
      // pp dphi
      double dphi = (*pp_jl).Phi() - (*pp_js).Phi();
      
      //rotate dphi to be in [0,2*pi]
      while (dphi < 0)
      dphi += 2 * TMath::Pi();
      while (dphi > 2.0 * TMath::Pi())
      dphi -= 2.0 * TMath::Pi();
      
      pp_dphi[key]->Fill(cent_bin, dphi);
      
      // and grefmult/npart
      pp_npart[key]->Fill(cent_bin, *pp_nprt);
      
      // pp only now
      if ((*pp_only_jl).Pt() > 0) {
        pp_lead_hard_match_fraction[key]->Fill(cent_bin, 1);
        if ((*pp_only_js).Pt() > 0) {
          pp_sub_hard_match_fraction[key]->Fill(cent_bin, 1);
          pp_only_hard_lead_pt[key]->Fill(cent_bin, (*pp_only_jl).Pt());
          pp_only_hard_sub_pt[key]->Fill(cent_bin, (*pp_only_js).Pt());
          pp_only_match_lead_pt[key]->Fill(cent_bin, (*pp_only_jlm).Pt());
          pp_only_match_sub_pt[key]->Fill(cent_bin, (*pp_only_jsm).Pt());
          pp_only_hard_aj[key]->Fill(cent_bin,
                                     fabs((*pp_only_jl).Pt() - (*pp_only_js).Pt()) / ((*pp_only_jl).Pt() + (*pp_only_js).Pt()));
          pp_only_match_aj[key]->Fill(cent_bin,
                                      fabs((*pp_only_jlm).Pt() - (*pp_only_jsm).Pt()) / ((*pp_only_jlm).Pt() + (*pp_only_jsm).Pt()));
          pp_hard_lead_dpt[key]->Fill(cent_bin, (*pp_jl).Pt() - (*pp_only_jl).Pt());
          pp_hard_sub_dpt[key]->Fill(cent_bin, (*pp_js).Pt() - (*pp_only_js).Pt());
          pp_match_lead_dpt[key]->Fill(cent_bin, (*pp_jlm).Pt() - (*pp_only_jlm).Pt());
          pp_match_sub_dpt[key]->Fill(cent_bin, (*pp_jsm).Pt() - (*pp_only_jsm).Pt());
          pp_hard_lead_dpt_frac[key]->Fill(cent_bin, ((*pp_jl).Pt() - (*pp_only_jl).Pt()) / (*pp_jl).Pt());
          pp_hard_sub_dpt_frac[key]->Fill(cent_bin, ((*pp_js).Pt() - (*pp_only_js).Pt()) / (*pp_js).Pt());
          pp_match_lead_dpt_frac[key]->Fill(cent_bin, ((*pp_jlm).Pt() - (*pp_only_jlm).Pt()) / (*pp_jlm).Pt());
          pp_match_sub_dpt_frac[key]->Fill(cent_bin, ((*pp_jsm).Pt() - (*pp_only_jsm).Pt()) / (*pp_jsm).Pt());
          
          pp_hard_lead_pp_3d[key]->Fill(cent_bin, (*pp_only_jl).Pt(), (*pp_jl).Pt());
          pp_hard_sub_pp_3d[key]->Fill(cent_bin, (*pp_only_js).Pt(), (*pp_js).Pt());
          pp_match_lead_pp_3d[key]->Fill(cent_bin, (*pp_only_jlm).Pt(), (*pp_jlm).Pt());
          pp_match_sub_pp_3d[key]->Fill(cent_bin, (*pp_only_jsm).Pt(), (*pp_jsm).Pt());
        }
        else {
          pp_sub_hard_match_fraction[key]->Fill(cent_bin, 0);
        }
      }
      else {
        pp_lead_hard_match_fraction[key]->Fill(cent_bin, 0);
      }
    }
    
    // now for the systematics
    while (tow_p_reader.Next()) {
      
      int cent_bin = -1;
      for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
        auto& pair = cent_bin_boundaries[i];
        if (*tow_p_cent >= pair.first &&
            *tow_p_cent <= pair.second) {
          cent_bin = i;
          break;
        }
      }
      
      if (cent_bin == -1)
        continue;
      
      pp_hard_aj_tow_p[key]->Fill(cent_bin,
                            fabs((*tow_p_jl).Pt() - (*tow_p_js).Pt()) / ((*tow_p_jl).Pt() + (*tow_p_js).Pt()));
      pp_match_aj_tow_p[key]->Fill(cent_bin,
                             fabs((*tow_p_jlm).Pt() - (*tow_p_jsm).Pt()) / ((*tow_p_jlm).Pt() + (*tow_p_jsm).Pt()));
      pp_hard_aj_tow_p_test[key]->Fill(cent_bin,
                                  fabs((*tow_p_jl).Pt() - (*tow_p_js).Pt()) / ((*tow_p_jl).Pt() + (*tow_p_js).Pt()));
      pp_match_aj_tow_p_test[key]->Fill(cent_bin,
                                   fabs((*tow_p_jlm).Pt() - (*tow_p_jsm).Pt()) / ((*tow_p_jlm).Pt() + (*tow_p_jsm).Pt()));
    }
    
    while (tow_m_reader.Next()) {
      
      int cent_bin = -1;
      for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
        auto& pair = cent_bin_boundaries[i];
        if (*tow_m_cent >= pair.first &&
            *tow_m_cent <= pair.second) {
          cent_bin = i;
          break;
        }
      }
      
      if (cent_bin == -1)
        continue;
      
      pp_hard_aj_tow_m[key]->Fill(cent_bin,
                            fabs((*tow_m_jl).Pt() - (*tow_m_js).Pt()) / ((*tow_m_jl).Pt() + (*tow_m_js).Pt()));
      pp_match_aj_tow_m[key]->Fill(cent_bin,
                             fabs((*tow_m_jlm).Pt() - (*tow_m_jsm).Pt()) / ((*tow_m_jlm).Pt() + (*tow_m_jsm).Pt()));
      pp_hard_aj_tow_m_test[key]->Fill(cent_bin,
                                       fabs((*tow_m_jl).Pt() - (*tow_m_js).Pt()) / ((*tow_m_jl).Pt() + (*tow_m_js).Pt()));
      pp_match_aj_tow_m_test[key]->Fill(cent_bin,
                                        fabs((*tow_m_jlm).Pt() - (*tow_m_jsm).Pt()) / ((*tow_m_jlm).Pt() + (*tow_m_jsm).Pt()));
    }
    
    while (track_p_reader.Next()) {
      
      int cent_bin = -1;
      for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
        auto& pair = cent_bin_boundaries[i];
        if (*track_p_cent >= pair.first &&
            *track_p_cent <= pair.second) {
          cent_bin = i;
          break;
        }
      }
      
      if (cent_bin == -1)
        continue;
      
      pp_hard_aj_track_p[key]->Fill(cent_bin,
                            fabs((*track_p_jl).Pt() - (*track_p_js).Pt()) / ((*track_p_jl).Pt() + (*track_p_js).Pt()));
      pp_match_aj_track_p[key]->Fill(cent_bin,
                             fabs((*track_p_jlm).Pt() - (*track_p_jsm).Pt()) / ((*track_p_jlm).Pt() + (*track_p_jsm).Pt()));
      pp_hard_aj_track_p_test[key]->Fill(cent_bin,
                                         fabs((*track_p_jl).Pt() - (*track_p_js).Pt()) / ((*track_p_jl).Pt() + (*track_p_js).Pt()));
      pp_match_aj_track_p_test[key]->Fill(cent_bin,
                                          fabs((*track_p_jlm).Pt() - (*track_p_jsm).Pt()) / ((*track_p_jlm).Pt() + (*track_p_jsm).Pt()));
    }
    
    while (track_m_reader.Next()) {
      
      int cent_bin = -1;
      for (int i = 0; i < cent_bin_boundaries.size(); ++i) {
        auto& pair = cent_bin_boundaries[i];
        if (*track_m_cent >= pair.first &&
            *track_m_cent <= pair.second) {
          cent_bin = i;
          break;
        }
      }
      
      if (cent_bin == -1)
        continue;

      
      pp_hard_aj_track_m[key]->Fill(cent_bin,
                            fabs((*track_m_jl).Pt() - (*track_m_js).Pt()) / ((*track_m_jl).Pt() + (*track_m_js).Pt()));
      pp_match_aj_track_m[key]->Fill(cent_bin,
                             fabs((*track_m_jlm).Pt() - (*track_m_jsm).Pt()) / ((*track_m_jlm).Pt() + (*track_m_jsm).Pt()));
      pp_hard_aj_track_m_test[key]->Fill(cent_bin,
                                    fabs((*track_m_jl).Pt() - (*track_m_js).Pt()) / ((*track_m_jl).Pt() + (*track_m_js).Pt()));
      pp_match_aj_track_m_test[key]->Fill(cent_bin,
                                     fabs((*track_m_jlm).Pt() - (*track_m_jsm).Pt()) / ((*track_m_jlm).Pt() + (*track_m_jsm).Pt()));
    }
    
    
    
    // split by centrality
    auau_npart_cent[key] = SplitByCentralityNormalized(auau_npart[key]);
    pp_npart_cent[key] = SplitByCentralityNormalized(pp_npart[key]);
    
    auau_hard_lead_pt_cent[key] = SplitByCentralityNormalized(auau_hard_lead_pt[key]);
    auau_hard_sub_pt_cent[key] = SplitByCentralityNormalized(auau_hard_sub_pt[key]);
    auau_match_lead_pt_cent[key] = SplitByCentralityNormalized(auau_match_lead_pt[key]);
    auau_match_sub_pt_cent[key] = SplitByCentralityNormalized(auau_match_sub_pt[key]);
    auau_hard_aj_cent[key] = SplitByCentralityNormalized(auau_hard_aj[key]);
    auau_match_aj_cent[key] = SplitByCentralityNormalized(auau_match_aj[key]);
    auau_hard_aj_test_cent[key] = SplitByCentralityNormalized(auau_hard_aj_test[key]);
    auau_match_aj_test_cent[key] = SplitByCentralityNormalized(auau_match_aj_test[key]);
    auau_off_axis_lead_pt_cent[key] = SplitByCentralityNormalized(auau_off_axis_lead_pt[key]);
    auau_off_axis_sub_pt_cent[key] = SplitByCentralityNormalized(auau_off_axis_sub_pt[key]);
    auau_off_axis_aj_cent[key] = SplitByCentralityNormalized(auau_off_axis_aj[key]);
    
    pp_hard_lead_pt_cent[key] = SplitByCentralityNormalized(pp_hard_lead_pt[key]);
    pp_hard_sub_pt_cent[key] = SplitByCentralityNormalized(pp_hard_sub_pt[key]);
    pp_match_lead_pt_cent[key] = SplitByCentralityNormalized(pp_match_lead_pt[key]);
    pp_match_sub_pt_cent[key] = SplitByCentralityNormalized(pp_match_sub_pt[key]);
    pp_hard_aj_cent[key] = SplitByCentralityNormalized(pp_hard_aj[key]);
    pp_match_aj_cent[key] = SplitByCentralityNormalized(pp_match_aj[key]);
    pp_hard_aj_test_cent[key] = SplitByCentralityNormalized(pp_hard_aj_test[key]);
    pp_match_aj_test_cent[key] = SplitByCentralityNormalized(pp_match_aj_test[key]);
    
    pp_only_hard_lead_pt_cent[key] = SplitByCentralityNormalized(pp_only_hard_lead_pt[key]);
    pp_only_hard_sub_pt_cent[key] = SplitByCentralityNormalized(pp_only_hard_sub_pt[key]);
    pp_only_match_lead_pt_cent[key] = SplitByCentralityNormalized(pp_only_match_lead_pt[key]);
    pp_only_match_sub_pt_cent[key] = SplitByCentralityNormalized(pp_only_match_sub_pt[key]);
    pp_only_hard_aj_cent[key] = SplitByCentralityNormalized(pp_only_hard_aj[key]);
    pp_only_match_aj_cent[key] = SplitByCentralityNormalized(pp_only_match_aj[key]);
    
    pp_hard_lead_dpt_cent[key] = SplitByCentralityNormalized(pp_hard_lead_dpt[key]);
    pp_hard_sub_dpt_cent[key] = SplitByCentralityNormalized(pp_hard_sub_dpt[key]);
    pp_match_lead_dpt_cent[key] = SplitByCentralityNormalized(pp_match_lead_dpt[key]);
    pp_match_sub_dpt_cent[key] = SplitByCentralityNormalized(pp_match_sub_dpt[key]);
    
    pp_hard_lead_pp_3d_cent[key] = SplitByCentrality3D(pp_hard_lead_pp_3d[key]);
    pp_hard_sub_pp_3d_cent[key] = SplitByCentrality3D(pp_hard_sub_pp_3d[key]);
    pp_match_lead_pp_3d_cent[key] = SplitByCentrality3D(pp_match_lead_pp_3d[key]);
    pp_match_sub_pp_3d_cent[key] = SplitByCentrality3D(pp_match_sub_pp_3d[key]);
    
    pp_lead_hard_match_fraction_cent[key] = SplitByCentralityNormalized(pp_lead_hard_match_fraction[key]);
    pp_sub_hard_match_fraction_cent[key] = SplitByCentralityNormalized(pp_sub_hard_match_fraction[key]);
    
    pp_hard_lead_dpt_frac_cent[key] = SplitByCentralityNormalized(pp_hard_lead_dpt_frac[key]);
    pp_hard_sub_dpt_frac_cent[key] = SplitByCentralityNormalized(pp_hard_sub_dpt_frac[key]);
    pp_match_lead_dpt_frac_cent[key] = SplitByCentralityNormalized(pp_match_lead_dpt_frac[key]);
    pp_match_sub_dpt_frac_cent[key] = SplitByCentralityNormalized(pp_match_sub_dpt_frac[key]);
    
    auau_hard_lead_const_cent[key] = SplitByCentralityNormalized(auau_hard_lead_const[key]);
    auau_hard_sub_const_cent[key] = SplitByCentralityNormalized(auau_hard_sub_const[key]);
    auau_match_lead_const_cent[key] = SplitByCentralityNormalized(auau_match_lead_const[key]);
    auau_match_sub_const_cent[key] = SplitByCentralityNormalized(auau_match_sub_const[key]);
    auau_off_axis_lead_const_cent[key] = SplitByCentralityNormalized(auau_off_axis_lead_const[key]);
    auau_off_axis_sub_const_cent[key] = SplitByCentralityNormalized(auau_off_axis_sub_const[key]);
    
    pp_hard_lead_const_cent[key] = SplitByCentralityNormalized(pp_hard_lead_const[key]);
    pp_hard_sub_const_cent[key] = SplitByCentralityNormalized(pp_hard_sub_const[key]);
    pp_match_lead_const_cent[key] = SplitByCentralityNormalized(pp_match_lead_const[key]);
    pp_match_sub_const_cent[key] = SplitByCentralityNormalized(pp_match_sub_const[key]);
    
    auau_hard_lead_rho_cent[key] = SplitByCentralityNormalized(auau_hard_lead_rho[key]);
    auau_hard_sub_rho_cent[key] = SplitByCentralityNormalized(auau_hard_sub_rho[key]);
    auau_match_lead_rho_cent[key] = SplitByCentralityNormalized(auau_match_lead_rho[key]);
    auau_match_sub_rho_cent[key] = SplitByCentralityNormalized(auau_match_sub_rho[key]);
    
    pp_hard_lead_rho_cent[key] = SplitByCentralityNormalized(pp_hard_lead_rho[key]);
    pp_hard_sub_rho_cent[key] = SplitByCentralityNormalized(pp_hard_sub_rho[key]);
    pp_match_lead_rho_cent[key] = SplitByCentralityNormalized(pp_match_lead_rho[key]);
    pp_match_sub_rho_cent[key] = SplitByCentralityNormalized(pp_match_sub_rho[key]);
    
    auau_hard_lead_sig_cent[key] = SplitByCentralityNormalized(auau_hard_lead_sig[key]);
    auau_hard_sub_sig_cent[key] = SplitByCentralityNormalized(auau_hard_sub_sig[key]);
    auau_match_lead_sig_cent[key] = SplitByCentralityNormalized(auau_match_lead_sig[key]);
    auau_match_sub_sig_cent[key] = SplitByCentralityNormalized(auau_match_sub_sig[key]);
    
    pp_hard_lead_sig_cent[key] = SplitByCentralityNormalized(pp_hard_lead_sig[key]);
    pp_hard_sub_sig_cent[key] = SplitByCentralityNormalized(pp_hard_sub_sig[key]);
    pp_match_lead_sig_cent[key] = SplitByCentralityNormalized(pp_match_lead_sig[key]);
    pp_match_sub_sig_cent[key] = SplitByCentralityNormalized(pp_match_sub_sig[key]);
    
    auau_off_axis_lead_rho_cent[key] = SplitByCentralityNormalized(auau_off_axis_lead_rho[key]);
    auau_off_axis_sub_rho_cent[key] = SplitByCentralityNormalized(auau_off_axis_sub_rho[key]);
    auau_off_axis_lead_sig_cent[key] = SplitByCentralityNormalized(auau_off_axis_lead_sig[key]);
    auau_off_axis_sub_sig_cent[key] = SplitByCentralityNormalized(auau_off_axis_sub_sig[key]);
    
    // now for the systematics
    
    pp_hard_aj_tow_p_cent[key] = SplitByCentralityNormalized(pp_hard_aj_tow_p[key]);
    pp_match_aj_tow_p_cent[key] = SplitByCentralityNormalized(pp_match_aj_tow_p[key]);
    pp_hard_aj_tow_m_cent[key] = SplitByCentralityNormalized(pp_hard_aj_tow_m[key]);
    pp_match_aj_tow_m_cent[key] = SplitByCentralityNormalized(pp_match_aj_tow_m[key]);
    
    pp_hard_aj_tow_p_test_cent[key] = SplitByCentralityNormalized(pp_hard_aj_tow_p_test[key]);
    pp_match_aj_tow_p_test_cent[key] = SplitByCentralityNormalized(pp_match_aj_tow_p_test[key]);
    pp_hard_aj_tow_m_test_cent[key] = SplitByCentralityNormalized(pp_hard_aj_tow_m_test[key]);
    pp_match_aj_tow_m_test_cent[key] = SplitByCentralityNormalized(pp_match_aj_tow_m_test[key]);
    
    pp_hard_aj_track_p_cent[key] = SplitByCentralityNormalized(pp_hard_aj_track_p[key]);
    pp_match_aj_track_p_cent[key] = SplitByCentralityNormalized(pp_match_aj_track_p[key]);
    pp_hard_aj_track_m_cent[key] = SplitByCentralityNormalized(pp_hard_aj_track_m[key]);
    pp_match_aj_track_m_cent[key] = SplitByCentralityNormalized(pp_match_aj_track_m[key]);
    
    pp_hard_aj_track_p_test_cent[key] = SplitByCentralityNormalized(pp_hard_aj_track_p_test[key]);
    pp_match_aj_track_p_test_cent[key] = SplitByCentralityNormalized(pp_match_aj_track_p_test[key]);
    pp_hard_aj_track_m_test_cent[key] = SplitByCentralityNormalized(pp_hard_aj_track_m_test[key]);
    pp_match_aj_track_m_test_cent[key] = SplitByCentralityNormalized(pp_match_aj_track_m_test[key]);
    
    for (int i = 0; i < auau_hard_aj_cent[key].size(); ++i) {
      auau_hard_aj_cent[key][i]->RebinX(2);
      auau_match_aj_cent[key][i]->RebinX(2);
      auau_off_axis_aj_cent[key][i]->RebinX(2);
      pp_hard_aj_cent[key][i]->RebinX(2);
      pp_match_aj_cent[key][i]->RebinX(2);
      pp_hard_aj_tow_p_cent[key][i]->RebinX(2);
      pp_match_aj_tow_p_cent[key][i]->RebinX(2);
      pp_hard_aj_tow_m_cent[key][i]->RebinX(2);
      pp_match_aj_tow_m_cent[key][i]->RebinX(2);
      pp_hard_aj_track_p_cent[key][i]->RebinX(2);
      pp_match_aj_track_p_cent[key][i]->RebinX(2);
      pp_hard_aj_track_m_cent[key][i]->RebinX(2);
      pp_match_aj_track_m_cent[key][i]->RebinX(2);
    }
    
    Overlay1D(auau_npart_cent[key], refcent_string, hopts, copts, out_loc, "auau_npart_cent",
              "", "N_{part}", "fraction", "Centrality");
    Overlay1D(auau_hard_aj_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_aj",
              "", "A_{J}", "fraction", "Centrality");
    Overlay1D(auau_match_aj_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_aj",
              "", "A_{J}", "fraction", "Centrality");
    Overlay1D(auau_off_axis_aj_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_aj",
              "", "A_{J}", "fraction", "Centrality");
    Overlay1D(auau_hard_lead_pt_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_lead_pt",
              "", "p_{T}", "fraction", "Centrality");
    Overlay1D(auau_match_lead_pt_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_lead_pt",
              "", "p_{T}", "fraction", "Centrality");
    Overlay1D(auau_hard_sub_pt_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_sub_pt",
              "", "p_{T}", "fraction", "Centrality");
    Overlay1D(auau_match_sub_pt_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_sub_pt",
              "", "p_{T}", "fraction", "Centrality");
    Overlay1D(auau_off_axis_lead_pt_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_lead_pt",
              "", "p_{T}", "fraction", "Centrality");
    Overlay1D(auau_off_axis_sub_pt_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_sub_pt",
              "", "p_{T}", "fraction", "Centrality");
    
    Overlay1D(pp_npart_cent[key], refcent_string, hopts, copts, out_loc, "pp_npart_cent",
              "", "N_{part}", "fraction", "Centrality");
    Overlay1D(pp_hard_aj_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_aj",
              "", "A_{J}", "fraction", "Centrality");
    Overlay1D(pp_match_aj_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_aj",
              "", "A_{J}", "fraction", "Centrality");
    Overlay1D(pp_hard_lead_pt_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_lead_pt",
              "", "p_{T}", "fraction", "Centrality");
    Overlay1D(pp_match_lead_pt_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_lead_pt",
              "", "p_{T}", "fraction", "Centrality");
    Overlay1D(pp_hard_sub_pt_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_sub_pt",
              "", "p_{T}", "fraction", "Centrality");
    Overlay1D(pp_match_sub_pt_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_sub_pt",
              "", "p_{T}", "fraction", "Centrality");
    
    Overlay1D(auau_hard_lead_const_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_lead_const",
              "", "N_{constituents}", "fraction");
    Overlay1D(auau_hard_sub_const_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_sub_const",
              "", "N_{constituents}", "fraction");
    Overlay1D(auau_match_lead_const_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_lead_const",
              "", "N_{constituents}", "fraction");
    Overlay1D(auau_match_sub_const_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_sub_const",
              "", "N_{constituents}", "fraction");
    Overlay1D(auau_off_axis_lead_const_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_lead_const",
              "", "N_{constituents}", "fraction");
    Overlay1D(auau_off_axis_sub_const_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_sub_const",
              "", "N_{constituents}", "fraction");
    
    Overlay1D(pp_hard_lead_const_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_lead_const",
              "", "N_{constituents}", "fraction");
    Overlay1D(pp_hard_sub_const_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_sub_const",
              "", "N_{constituents}", "fraction");
    Overlay1D(pp_match_lead_const_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_lead_const",
              "", "N_{constituents}", "fraction");
    Overlay1D(pp_match_sub_const_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_sub_const",
              "", "N_{constituents}", "fraction");
    
    Overlay1D(auau_hard_lead_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_lead_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_hard_sub_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_sub_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_match_lead_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_lead_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_match_sub_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_sub_rho",
              "", "#rho", "fraction");
    
    Overlay1D(pp_hard_lead_rho_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_lead_rho",
              "", "#rho", "fraction");
    Overlay1D(pp_hard_sub_rho_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_sub_rho",
              "", "#rho", "fraction");
    Overlay1D(pp_match_lead_rho_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_lead_rho",
              "", "#rho", "fraction");
    Overlay1D(pp_match_sub_rho_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_sub_rho",
              "", "#rho", "fraction");
    
    Overlay1D(auau_hard_lead_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_lead_sig",
              "", "#sigma", "fraction");
    Overlay1D(auau_hard_sub_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_sub_sig",
              "", "#sigma", "fraction");
    Overlay1D(auau_match_lead_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_lead_sig",
              "", "#sigma", "fraction");
    Overlay1D(auau_match_sub_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_sub_sig",
              "", "#sigma", "fraction");
    
    Overlay1D(pp_hard_lead_sig_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_lead_sig",
              "", "#sigma", "fraction");
    Overlay1D(pp_hard_sub_sig_cent[key], refcent_string, hopts, copts, out_loc, "pp_hard_sub_sig",
              "", "#sigma", "fraction");
    Overlay1D(pp_match_lead_sig_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_lead_sig",
              "", "#sigma", "fraction");
    Overlay1D(pp_match_sub_sig_cent[key], refcent_string, hopts, copts, out_loc, "pp_match_sub_sig",
              "", "#sigma", "fraction");
    
    Overlay1D(auau_off_axis_lead_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_lead_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_off_axis_sub_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_sub_sub_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_off_axis_lead_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_lead_sig",
              "", "#sigma", "fraction");
    Overlay1D(auau_off_axis_sub_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_sub_sig",
              "", "#sigma", "fraction");
    
    for (int i = 0; i < auau_hard_aj_cent[key].size(); ++i) {
      std::vector<string> cent_name{"AuAu", "pp"};
      std::vector<string> cent_name_match{"AuAu", "pp", "hard-core embed"};
      std::vector<TH1D*> cent_list{auau_hard_aj_cent[key][i], pp_hard_aj_cent[key][i]};
      std::vector<TH1D*> cent_list_match{auau_match_aj_cent[key][i], pp_match_aj_cent[key][i], auau_off_axis_aj_cent[key][i]};
      
      string out_loc_cent = FLAGS_outputDir + "/" + key + "/" + "cent_" + std::to_string(i);
      boost::filesystem::path dir(out_loc_cent.c_str());
      boost::filesystem::create_directories(dir);
      
      Print2DSimple(pp_hard_lead_pp_3d_cent[key][i], hopts, coptslogz, out_loc_cent, "ppleadhardpt2d", "", "pp p_{T}", "embedded pp p_{T}");
      Print2DSimple(pp_hard_sub_pp_3d_cent[key][i], hopts, coptslogz, out_loc_cent, "ppsubhardpt2d", "", "pp p_{T}", "embedded pp p_{T}");
      Print2DSimple(pp_match_lead_pp_3d_cent[key][i], hopts, coptslogz, out_loc_cent, "ppleadmatchpt2d", "", "pp p_{T}", "embedded pp p_{T}");
      Print2DSimple(pp_match_sub_pp_3d_cent[key][i], hopts, coptslogz, out_loc_cent, "ppsubmatchpt2d", "", "pp p_{T}", "embedded pp p_{T}");
      
      Overlay1D(cent_list_match, cent_name_match, hopts, copts, out_loc_cent, "match_aj_with_off_axis", "", "A_{J}", "fraction");
      
      Overlay1D(pp_hard_aj_tow_p_cent[key][i], pp_hard_aj_tow_m_cent[key][i], "pp tow_p", "pp tow_m", hopts, copts, out_loc_cent, "tow",
                "", "AJ", "fraction");
      Overlay1D(pp_hard_aj_track_p_cent[key][i], pp_hard_aj_track_m_cent[key][i], "pp track_p", "pp track_m", hopts, copts, out_loc_cent, "track",
                "", "AJ", "fraction");
      Overlay1D(pp_match_aj_tow_p_cent[key][i], pp_match_aj_tow_m_cent[key][i], "pp tow_p", "pp tow_m", hopts, copts, out_loc_cent, "towmatch",
                "", "AJ", "fraction");
      Overlay1D(pp_match_aj_track_p_cent[key][i], pp_match_aj_track_m_cent[key][i], "pp track_p", "pp track_m", hopts, copts, out_loc_cent, "trackmatch",
                "", "AJ", "fraction");
      
      systematic_errors_hard[key].push_back( GetSystematic(pp_hard_aj_cent[key][i], pp_hard_aj_tow_p_cent[key][i], pp_hard_aj_tow_m_cent[key][i],
                                             pp_hard_aj_track_p_cent[key][i], pp_hard_aj_track_m_cent[key][i]));
      systematic_errors_match[key].push_back( GetSystematic(pp_match_aj_cent[key][i], pp_match_aj_tow_p_cent[key][i], pp_match_aj_tow_m_cent[key][i],
                                              pp_match_aj_track_p_cent[key][i], pp_match_aj_track_m_cent[key][i]));
      
      Overlay1D(auau_hard_aj_cent[key][i], pp_hard_aj_cent[key][i], systematic_errors_hard[key][i], 0.0, 0.25, 0.0, 0.9, "AuAu hard A_{J}", "PP hard A_{J}",
                "", hopts, copts, out_loc_cent, "aj_hard", "", "A_{J}", "fraction");
      Overlay1D(auau_match_aj_cent[key][i], pp_match_aj_cent[key][i], systematic_errors_match[key][i], 0.0, 0.3, 0.0, 0.9, "AuAu matched A_{J}", "PP matched A_{J}",
                "systematics", hopts, copts, out_loc_cent, "aj_match", "", "A_{J}", "fraction");
      
      Overlay1D(pp_hard_lead_dpt_cent[key][i], pp_hard_sub_dpt_cent[key][i], "leading jet", "subleading jet", hopts, copts, out_loc_cent, "pp_dpt",
                "", "dp_{T}", "fraction");
      Overlay1D(pp_hard_lead_dpt_frac_cent[key][i], pp_hard_sub_dpt_frac_cent[key][i], "leading jet", "subleading jet", hopts, copts, out_loc_cent, "pp_dpt_frac",
                "", "dp_{T}/p_{T}", "fraction");
      
      // print the fractions of successfully matched jets from embedded pp -> pp
      PrettyPrint1D(pp_lead_hard_match_fraction_cent[key][i], hopts, copts, "", out_loc_cent, "leadhardmatchingfraction", "no match / match",
                    "fraction", "");
      PrettyPrint1D(pp_sub_hard_match_fraction_cent[key][i], hopts, copts, "", out_loc_cent, "leadsubmatchingfraction", "no match / match",
                    "fraction", "");
    }
  }
  
  for (int cent = 0; cent < cent_boundaries.size(); ++cent) {
    // make a directory for our grid outputs
    string out_loc_grid = FLAGS_outputDir + "/grid_cent_" + std::to_string(cent);
    
    // build output directory if it doesn't exist, using boost::filesystem
    boost::filesystem::path dir_cent(out_loc_grid.c_str());
    boost::filesystem::create_directories(dir_cent);
    
    // now build histograms
    TH2D* p_values_hard = new TH2D(dijetcore::MakeString("p_values_hard_", cent).c_str(), ";R;p_{T}^{const}",
                                   radii.size(), - 0.5, radii.size() - 0.5,
                                   constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_hard = new TH2D(dijetcore::MakeString("ks_values_hard_", cent).c_str(), ";R;p_{T}^{const}",
                                    radii.size(), - 0.5, radii.size() - 0.5,
                                    constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* p_values_match = new TH2D(dijetcore::MakeString("p_values_match_", cent).c_str(), ";R;p_{T}^{const}",
                                    radii.size(), - 0.5, radii.size() - 0.5,
                                    constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_match = new TH2D(dijetcore::MakeString("ks_values_match_", cent).c_str(), ";R;p_{T}^{const}",
                                     radii.size(), - 0.5, radii.size() - 0.5,
                                     constpt.size(), -0.5, constpt.size() - 0.5);
    
    TH2D* p_values_hard_error = new TH2D(dijetcore::MakeString("p_values_hard_error_", cent).c_str(), ";R;p_{T}^{const}",
                                         radii.size(), - 0.5, radii.size() - 0.5,
                                         constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_hard_error = new TH2D(dijetcore::MakeString("ks_values_hard_error_", cent).c_str(), ";R;p_{T}^{const}",
                                          radii.size(), - 0.5, radii.size() - 0.5,
                                          constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* p_values_match_error = new TH2D(dijetcore::MakeString("p_values_match_error_", cent).c_str(), ";R;p_{T}^{const}",
                                          radii.size(), - 0.5, radii.size() - 0.5,
                                          constpt.size(), -0.5, constpt.size() - 0.5);
    TH2D* ks_values_match_error = new TH2D(dijetcore::MakeString("ks_values_match_error_", cent).c_str(), ";R;p_{T}^{const}",
                                           radii.size(), - 0.5, radii.size() - 0.5,
                                           constpt.size(), -0.5, constpt.size() - 0.5);
    
    for (int j = 0; j < radii.size(); ++j) {
      p_values_hard->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
      ks_values_hard->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
      p_values_match->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
      ks_values_match->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
      p_values_hard_error->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
      ks_values_hard_error->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
      p_values_match_error->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
      ks_values_match_error->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
    }
    for (int j = 0; j < radii.size(); ++j) {
      p_values_hard->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
      ks_values_hard->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
      p_values_match->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
      ks_values_match->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
      p_values_hard_error->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
      ks_values_hard_error->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
      p_values_match_error->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
      ks_values_match_error->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
    }
    
    TCanvas* c_hard = new TCanvas(dijetcore::MakeString("hard_canvas_", cent).c_str());
    std::vector<std::vector<TPad*>> hard_pads = dijetcore::CanvasPartition(c_hard, radii.size(), constpt.size(), 0.10, 0.10, 0.05, 0.05);
    
    TCanvas* c_match = new TCanvas(dijetcore::MakeString("match_canvas_", cent).c_str());
    std::vector<std::vector<TPad*>> matched_pads = dijetcore::CanvasPartition(c_match, radii.size(), constpt.size(), 0.10, 0.10, 0.05, 0.05);
    for (int rad = 0; rad < radii.size(); ++rad) {
      for (int pt = 0; pt < constpt.size(); ++pt) {
    
        int x_bin = p_values_hard->GetXaxis()->FindBin(rad);
        int y_bin = p_values_hard->GetYaxis()->FindBin(pt);
        string key = grid_keys[rad][pt];
        dijetcore::DijetKey key_params = grid_key_params[rad][pt];
        
        TH1D* aj_hard_test = auau_hard_aj_test_cent[key][cent];
        TH1D* aj_match_test = auau_match_aj_test_cent[key][cent];
        
        TH1D* aj_hard_pp_test = pp_hard_aj_test_cent[key][cent];
        TH1D* aj_match_pp_test = pp_match_aj_test_cent[key][cent];
        
        TH1D* aj_hard = auau_hard_aj_cent[key][cent];
        TH1D* aj_match = auau_match_aj_cent[key][cent];
        
        TH1D* aj_hard_pp = pp_hard_aj_cent[key][cent];
        TH1D* aj_match_pp = pp_match_aj_cent[key][cent];
        
        // print
        
        c_hard->cd(0);
        auau_hard_aj_cent[key][cent]->GetXaxis()->SetTitle(dijetcore::MakeString("R=", radii_string[rad]).c_str());
        auau_hard_aj_cent[key][cent]->GetYaxis()->SetTitle(dijetcore::MakeString("p_{T}^{const}=", constpt_string[pt]).c_str());
        auau_hard_aj_cent[key][cent]->SetMarkerSize(0.2);
        pp_hard_aj_cent[key][cent]->SetMarkerSize(0.2);
        TPad* hard_pad = hard_pads[rad][pt];
        hard_pad->Draw();
        hard_pad->SetFillStyle(4000);
        hard_pad->SetFrameFillStyle(4000);
        hard_pad->cd();
        auau_hard_aj_cent[key][cent]->Draw();
        pp_hard_aj_cent[key][cent]->Draw("SAME");
        c_hard->cd(0);
        
        c_match->cd(0);
        auau_match_aj_cent[key][cent]->GetXaxis()->SetTitle(dijetcore::MakeString("R=", radii_string[rad]).c_str());
        auau_match_aj_cent[key][cent]->GetYaxis()->SetTitle(dijetcore::MakeString("p_{T}^{const}=", constpt_string[rad]).c_str());
        auau_match_aj_cent[key][cent]->SetMarkerSize(0.2);
        pp_match_aj_cent[key][cent]->SetMarkerSize(0.2);
        TPad* match_pad = matched_pads[rad][pt];
        match_pad->SetFillStyle(4000);
        match_pad->SetFrameFillStyle(4000);
        match_pad->cd();
        auau_match_aj_cent[key][cent]->Draw();
        pp_match_aj_cent[key][cent]->Draw("SAME");
        c_match->cd(0);
        
      
        std::vector<TH1D*> aj_hard_pp_var{pp_hard_aj_tow_p_cent[key][cent], pp_hard_aj_tow_m_cent[key][cent], pp_hard_aj_track_p_cent[key][cent], pp_hard_aj_track_m_cent[key][cent]};
        std::vector<TH1D*> aj_match_pp_var{pp_match_aj_tow_p_cent[key][cent], pp_match_aj_tow_m_cent[key][cent], pp_match_aj_track_p_cent[key][cent], pp_match_aj_track_m_cent[key][cent]};
        
        std::vector<TH1D*> aj_hard_pp_var_test{pp_hard_aj_tow_p_test_cent[key][cent], pp_hard_aj_tow_m_test_cent[key][cent], pp_hard_aj_track_p_test_cent[key][cent], pp_hard_aj_track_m_test_cent[key][cent]};
        std::vector<TH1D*> aj_match_pp_var_test{pp_match_aj_tow_p_test_cent[key][cent], pp_match_aj_tow_m_test_cent[key][cent], pp_match_aj_track_p_test_cent[key][cent], pp_match_aj_track_m_test_cent[key][cent]};
        
        double p_value_hard = aj_hard->Chi2Test(aj_hard_pp, "UU NORM");
        double p_value_match = aj_match->Chi2Test(aj_match_pp, "UU NORM");
        double ks_value_hard = aj_hard_test->KolmogorovTest(aj_hard_pp_test);
        double ks_value_match = aj_match_test->KolmogorovTest(aj_match_pp_test);
        
        double max_hard_p_deviation = 0;
        double max_match_p_deviation = 0;
        double max_hard_ks_deviation = 0;
        double max_match_ks_deviation = 0;
        
        for (int i = 0; i < aj_hard_pp_var.size(); ++i) {
          double hard_p_deviation = aj_hard->Chi2Test(aj_hard_pp_var[i], "UU NORM");
          double hard_ks_deviation = aj_hard_test->KolmogorovTest(aj_hard_pp_var_test[i]);
          double match_p_deviation = aj_match->Chi2Test(aj_match_pp_var[i], "UU NORM");
          double match_ks_deviation = aj_match_test->KolmogorovTest(aj_match_pp_var_test[i]);
          
          if (fabs(p_value_hard - hard_p_deviation) > max_hard_p_deviation)
            max_hard_p_deviation = fabs(hard_p_deviation);
          if (fabs(ks_value_hard - hard_ks_deviation) > max_hard_ks_deviation)
            max_hard_ks_deviation = fabs(hard_ks_deviation);
          if (fabs(p_value_match - match_p_deviation) > max_match_p_deviation)
            max_match_p_deviation = fabs(match_p_deviation);
          if (fabs(ks_value_match - match_ks_deviation) > max_match_ks_deviation)
            max_match_ks_deviation = fabs(match_ks_deviation);
          
        }
        
        // set bins
        p_values_hard->SetBinContent(x_bin, y_bin, p_value_hard);
        p_values_match->SetBinContent(x_bin, y_bin, p_value_match);
        ks_values_hard->SetBinContent(x_bin, y_bin, ks_value_hard);
        ks_values_match->SetBinContent(x_bin, y_bin, ks_value_match);
        p_values_hard_error->SetBinContent(x_bin, y_bin, max_hard_p_deviation);
        p_values_match_error->SetBinContent(x_bin, y_bin, max_match_p_deviation);
        ks_values_hard_error->SetBinContent(x_bin, y_bin, max_hard_ks_deviation);
        ks_values_match_error->SetBinContent(x_bin, y_bin, max_match_ks_deviation);
        
      }
    }
    
    Print2DSimple(p_values_hard, hopts, copts, out_loc_grid, "ajphard", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(p_values_match, hopts, copts, out_loc_grid, "ajpmatch", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_hard, hopts, copts, out_loc_grid, "ajkshard", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_match, hopts, copts, out_loc_grid, "ajksmatch", "", "R", "p_{T}^{const}", "TEXT COLZ");
    
    Print2DSimple(p_values_hard_error, hopts, copts, out_loc_grid, "ajpharderror", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(p_values_match_error, hopts, copts, out_loc_grid, "ajpmatcherror", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_hard_error, hopts, copts, out_loc_grid, "ajksharderror", "", "R", "p_{T}^{const}", "TEXT COLZ");
    Print2DSimple(ks_values_match_error, hopts, copts, out_loc_grid, "ajksmatcherror", "", "R", "p_{T}^{const}", "TEXT COLZ");
    
    c_hard->SaveAs(dijetcore::MakeString(out_loc_grid, "/ajhardgrid.pdf").c_str());
    c_match->SaveAs(dijetcore::MakeString(out_loc_grid, "/ajmatchgrid.pdf").c_str());
    
  }
  
  return 0;
}
