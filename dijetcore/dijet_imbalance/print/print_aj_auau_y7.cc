#include "dijetcore/util/root/root_print_utils.h"
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

#include <unordered_map>
#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(outputDir, "tmp", "directory for output");

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

std::vector<TH1D*> SplitByCentrality(TH2D* h, int bins = 9) {
  
  std::vector<TH1D*> ret;
  for (int i = 1; i <= bins; ++i) {
    string name = string(h->GetName()) + std::to_string(i);
    TH1D* tmp = h->ProjectionY(name.c_str(), i, i);
    ret.push_back(tmp);
  }
  return ret;
}

std::vector<TH1D*> SplitByCentralityNormalized(TH2D* h, int bins = 9) {
  
  std::vector<TH1D*> ret;
  for (int i = 1; i <= bins; ++i) {
    string name = string(h->GetName()) + std::to_string(i);
    TH1D* tmp = h->ProjectionY(name.c_str(), i, i);
    if (tmp->Integral() > 0)
      tmp->Scale(1.0 / tmp->Integral());
    ret.push_back(tmp);
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
  
  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(FLAGS_input)) {
    std::cerr << "input file does not exist: " << FLAGS_input << std::endl;;
    return 1;
  }
  
  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
  FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);
  
  // read in the file
  TFile auau_file(FLAGS_input.c_str(), "READ");
  
  // define our centrality bins
  std::vector<unsigned> cent_boundaries = {485, 399, 269, 178, 114, 69, 39, 21, 10};
  std::vector<string> refcent_string{"0-5%", "5-10%", "10-20%", "20-30%",
    "30-40%", "40-50%", "50-60%", "60-70%", "70-80%"};
  
  
  // now we'll get the trees from the files, ignoring any objects
  // in the file that don't conform to the naming conventions from
  // the DijetWorker. There are also coincidence histograms to save
  std::vector<string> keys;
  std::unordered_map<std::string, TTree*> auau_trees;
  
  GetTreesFromFile(auau_file, auau_trees);
  
  for (auto& entry : auau_trees) {
    keys.push_back(entry.first);
  }
  
  std::unordered_map<std::string, TH1D*> auau_jl_pt;
  std::unordered_map<std::string, TH1D*> auau_jlm_pt;
  std::unordered_map<std::string, TH1D*> auau_js_pt;
  std::unordered_map<std::string, TH1D*> auau_jsm_pt;
  
  std::unordered_map<std::string, TH1D*> auau_jl_dr;
  std::unordered_map<std::string, TH1D*> auau_jlm_dr;
  std::unordered_map<std::string, TH1D*> auau_js_dr;
  std::unordered_map<std::string, TH1D*> auau_jsm_dr;
  
  for (auto key : keys) {
    auau_jl_pt.insert({key, (TH1D*) auau_file.Get(dijetcore::MakeString(key, "jlconstpt").c_str())});
    auau_jlm_pt.insert({key, (TH1D*) auau_file.Get(dijetcore::MakeString(key, "jlmconstpt").c_str())});
    auau_js_pt.insert({key, (TH1D*) auau_file.Get(dijetcore::MakeString(key, "jsconstpt").c_str())});
    auau_jsm_pt.insert({key, (TH1D*) auau_file.Get(dijetcore::MakeString(key, "jsmconstpt").c_str())});
    
    auau_jl_dr.insert({key, (TH1D*) auau_file.Get(dijetcore::MakeString(key, "jlconstdr").c_str())});
    auau_jlm_dr.insert({key, (TH1D*) auau_file.Get(dijetcore::MakeString(key, "jlmconstdr").c_str())});
    auau_js_dr.insert({key, (TH1D*) auau_file.Get(dijetcore::MakeString(key, "jsconstdr").c_str())});
    auau_jsm_dr.insert({key, (TH1D*) auau_file.Get(dijetcore::MakeString(key, "jsmconstdr").c_str())});
  }
  
  // save all histograms so we can do comparisons
  // between different keys
  std::unordered_map<string, TH2D*> auau_hard_lead_pt;
  std::unordered_map<string, TH2D*> auau_hard_sub_pt;
  std::unordered_map<string, TH2D*> auau_match_lead_pt;
  std::unordered_map<string, TH2D*> auau_match_sub_pt;
  std::unordered_map<string, TH2D*> auau_hard_aj;
  std::unordered_map<string, TH2D*> auau_match_aj;
  std::unordered_map<string, TH2D*> auau_dphi;
  std::unordered_map<string, TH2D*> auau_hard_lead_rp;
  std::unordered_map<string, TH2D*> auau_off_axis_lead_pt;
  std::unordered_map<string, TH2D*> auau_off_axis_sub_pt;
  std::unordered_map<string, TH2D*> auau_off_axis_aj;
  
  std::unordered_map<string, TH2D*> auau_hard_lead_const;
  std::unordered_map<string, TH2D*> auau_hard_sub_const;
  std::unordered_map<string, TH2D*> auau_match_lead_const;
  std::unordered_map<string, TH2D*> auau_match_sub_const;
  std::unordered_map<string, TH2D*> auau_off_axis_lead_const;
  std::unordered_map<string, TH2D*> auau_off_axis_sub_const;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_aj_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_lead_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_sub_pt_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_aj_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_lead_const_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_sub_const_cent;
  
  std::unordered_map<string, TH2D*> auau_hard_lead_rho;
  std::unordered_map<string, TH2D*> auau_hard_sub_rho;
  std::unordered_map<string, TH2D*> auau_match_lead_rho;
  std::unordered_map<string, TH2D*> auau_match_sub_rho;
  
  std::unordered_map<string, TH2D*> auau_hard_lead_sig;
  std::unordered_map<string, TH2D*> auau_hard_sub_sig;
  std::unordered_map<string, TH2D*> auau_match_lead_sig;
  std::unordered_map<string, TH2D*> auau_match_sub_sig;
  
  std::unordered_map<string, TH2D*> auau_off_axis_lead_rho;
  std::unordered_map<string, TH2D*> auau_off_axis_sub_rho;
  std::unordered_map<string, TH2D*> auau_off_axis_lead_sig;
  std::unordered_map<string, TH2D*> auau_off_axis_sub_sig;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_rho_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_hard_sub_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_match_sub_sig_cent;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_lead_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_sub_rho_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_lead_sig_cent;
  std::unordered_map<string, std::vector<TH1D*>> auau_off_axis_sub_sig_cent;
  
  std::unordered_map<string, TH2D*> auau_npart;
  
  std::unordered_map<string, std::vector<TH1D*>> auau_npart_cent;
  
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
    
    // first constituent histograms
    // pt & dr
    TH1D* auau_jl_pt_single = auau_jl_pt[key];
    TH1D* auau_jlm_pt_single = auau_jlm_pt[key];
    TH1D* auau_js_pt_single = auau_js_pt[key];
    TH1D* auau_jsm_pt_single = auau_jsm_pt[key];
    
    TH1D* auau_jl_dr_single = auau_jl_dr[key];
    TH1D* auau_jlm_dr_single = auau_jlm_dr[key];
    TH1D* auau_js_dr_single = auau_js_dr[key];
    TH1D* auau_jsm_dr_single = auau_jsm_dr[key];
    
    auau_jl_pt_single->Scale(1.0 / auau_jl_pt_single->Integral());
    auau_jlm_pt_single->Scale(1.0 / auau_jlm_pt_single->Integral());
    auau_js_pt_single->Scale(1.0 / auau_js_pt_single->Integral());
    auau_jsm_pt_single->Scale(1.0 / auau_jsm_pt_single->Integral());
    
    auau_jl_dr_single->Scale(1.0 / auau_jl_dr_single->Integral());
    auau_jlm_dr_single->Scale(1.0 / auau_jlm_dr_single->Integral());
    auau_js_dr_single->Scale(1.0 / auau_js_dr_single->Integral());
    auau_jsm_dr_single->Scale(1.0 / auau_jsm_dr_single->Integral());
    
    
    
    dijetcore::Overlay1D(auau_jl_pt_single, auau_js_pt_single, "leading jet", "subleading jet",
                         hopts, coptslogy, file_prefix, "const_pt_hard", "", "p_{T}", "fraction", "p_{T}^{const} > 2.0 GeV");
    dijetcore::Overlay1D(auau_jlm_pt_single, auau_jsm_pt_single, "leading jet", "subleading jet",
                         hopts, coptslogy, file_prefix, "const_pt_matched", "", "p_{T}", "fraction", "p_{T}^{const} > 0.2 GeV");
    dijetcore::Overlay1D(auau_jl_dr_single, auau_js_dr_single, "leading jet", "subleading jet",
                         hopts, copts, file_prefix, "const_dr_hard", "", "d_{R}", "fraction", "p_{T}^{const} > 2.0 GeV");
    dijetcore::Overlay1D(auau_jlm_dr_single, auau_jsm_dr_single, "leading jet", "subleading jet",
                         hopts, copts, file_prefix, "const_dr_matched", "", "d_{R}", "fraction", "p_{T}^{const} > 0.2 GeV");
    
    // process the trees
    TTree* auau_tree = auau_trees[key];
    TTreeReader auau_reader(auau_tree);
    
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
    
    // insert into the dictionaries
    auau_hard_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardleadpt").c_str(),
                                      "p_{T}", 9, -0.5, 8.5, 100, 0, 100);
    auau_hard_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardsubpt").c_str(),
                                    "p_{T}", 9, -0.5, 8.5, 100, 0, 100);
    auau_match_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchleadpt").c_str(),
                                       "p_{T}", 9, -0.5, 8.5, 100, 0, 100);
    auau_match_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchsubpt").c_str(),
                                      "p_{T}", 9, -0.5, 8.5, 100, 0, 100);
    auau_hard_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardaj").c_str(),
                                 "A_{J}", 9, -0.5, 8.5, 30, 0, 0.9);
    auau_match_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchaj").c_str(),
                                  "A_{J}", 9, -0.5, 8.5, 30, 0, 0.9);
    auau_dphi[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaudphi").c_str(),
                              "d#phi", 9, -0.5, 8.5, 100, 0, 2*TMath::Pi());
    
    auau_off_axis_lead_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisleadpt").c_str(),
                                          "p_{T}", 9, -0.5, 8.5, 100, 0, 100);
    auau_off_axis_sub_pt[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxissubpt").c_str(),
                                         "p_{T}", 9, -0.5, 8.5, 100, 0, 100);
    auau_off_axis_aj[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisaj").c_str(),
                                     "A_{J}", 9, -0.5, 8.5, 30, 0, 0.9);
    
    auau_hard_lead_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardleadconst").c_str(),
                                         "", 9, -0.5, 8.5, 100, 0, 100);
    auau_hard_sub_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardsubconst").c_str(),
                                        "", 9, -0.5, 8.5, 100, 0, 100);
    auau_match_lead_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchleadconst").c_str(),
                                          "", 9, -0.5, 8.5, 100, 0, 100);
    auau_match_sub_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchsubconst").c_str(),
                                         "", 9, -0.5, 8.5, 100, 0, 100);
    auau_off_axis_lead_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisleadconst").c_str(),
                                             "", 9, -0.5, 8.5, 100, 0, 100);
    auau_off_axis_sub_const[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxissubconst").c_str(),
                                            "", 9, -0.5, 8.5, 100, 0, 100);
    
    auau_hard_lead_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardleadrho").c_str(),
                                       "auau hard lead rho", 9, -0.5, 8.5, 100, 0, 100);
    auau_hard_lead_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardleadsig").c_str(),
                                       "auau hard lead sig", 9, -0.5, 8.5, 100, 0, 20);
    auau_hard_sub_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardsubrho").c_str(),
                                      "auau hard sub rho", 9, -0.5, 8.5, 100, 0, 100);
    auau_hard_sub_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauhardsubsig").c_str(),
                                      "auau hard sub sig", 9, -0.5, 8.5, 100, 0, 20);
    auau_match_lead_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchleadrho").c_str(),
                                        "auau match lead rho", 9, -0.5, 8.5, 100, 0, 100);
    auau_match_lead_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchleadsig").c_str(),
                                        "auau match lead sig", 9, -0.5, 8.5, 100, 0, 20);
    auau_match_sub_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchsubrho").c_str(),
                                       "auau match sub rho", 9, -0.5, 8.5, 100, 0, 100);
    auau_match_sub_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaumatchsubsig").c_str(),
                                       "auau match sub sig", 9, -0.5, 8.5, 100, 0, 20);
    
    auau_off_axis_lead_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisleadrho").c_str(),
                                           "", 9, -0.5, 8.5, 100, 0, 100);
    auau_off_axis_sub_rho[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxissubrho").c_str(),
                                          "", 9, -0.5, 8.5, 100, 0, 100);
    auau_off_axis_lead_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxisleadsig").c_str(),
                                           "", 9, -0.5, 8.5, 100, 0, 20);
    auau_off_axis_sub_sig[key] = new TH2D(dijetcore::MakeString(key_prefix, "auauoffaxissubsig").c_str(),
                                          "", 9, -0.5, 8.5, 100, 0, 20);
    
    auau_npart[key] = new TH2D(dijetcore::MakeString(key_prefix, "auaunpart").c_str(), ";refmult;nPart",
                               9, -0.5, 8.5, 100, 0.5, 2500.5);
    
    
    
    // loop over the data & fill histograms
    while (auau_reader.Next()) {
      
      // auau jet pt
      auau_hard_lead_pt[key]->Fill(*auau_cent, (*auau_jl).Pt());
      auau_hard_lead_rho[key]->Fill(*auau_cent, *auau_jlrho);
      auau_hard_lead_sig[key]->Fill(*auau_cent, *auau_jlsig);
      auau_hard_sub_pt[key]->Fill(*auau_cent, (*auau_js).Pt());
      auau_hard_sub_rho[key]->Fill(*auau_cent, *auau_jsrho);
      auau_hard_sub_sig[key]->Fill(*auau_cent, *auau_jssig);
      auau_match_lead_pt[key]->Fill(*auau_cent,(*auau_jlm).Pt());
      auau_match_lead_rho[key]->Fill(*auau_cent, *auau_jlmrho);
      auau_match_lead_sig[key]->Fill(*auau_cent, *auau_jlmsig);
      auau_match_sub_pt[key]->Fill(*auau_cent, (*auau_jsm).Pt());
      auau_match_sub_rho[key]->Fill(*auau_cent, *auau_jsmrho);
      auau_match_sub_sig[key]->Fill(*auau_cent, *auau_jsmsig);
      auau_hard_lead_const[key]->Fill(*auau_cent, *auau_jlconst);
      auau_hard_sub_const[key]->Fill(*auau_cent, *auau_jsconst);
      auau_match_lead_const[key]->Fill(*auau_cent, *auau_jlmconst);
      auau_match_sub_const[key]->Fill(*auau_cent, *auau_jsmconst);
      
      if (*auau_jloaconst != 0 && *auau_jsoaconst != 0) {
        auau_off_axis_lead_pt[key]->Fill(*auau_cent, (*auau_jloa).Pt());
        auau_off_axis_lead_rho[key]->Fill(*auau_cent, *auau_jloarho);
        auau_off_axis_lead_sig[key]->Fill(*auau_cent, *auau_jloasig);
        auau_off_axis_sub_pt[key]->Fill(*auau_cent, (*auau_jsoa).Pt());
        auau_off_axis_sub_rho[key]->Fill(*auau_cent, *auau_jsoarho);
        auau_off_axis_sub_sig[key]->Fill(*auau_cent, *auau_jsoasig);
        auau_off_axis_lead_const[key]->Fill(*auau_cent, *auau_jloaconst);
        auau_off_axis_sub_const[key]->Fill(*auau_cent, *auau_jsoaconst);
        auau_off_axis_aj[key]->Fill(*auau_cent,
                                    fabs((*auau_jloa).Pt() - (*auau_jsoa).Pt())/((*auau_jloa).Pt() + (*auau_jsoa).Pt()));
      }
      // auau Aj
      auau_hard_aj[key]->Fill(*auau_cent,
                           fabs((*auau_jl).Pt() - (*auau_js).Pt())/((*auau_jl).Pt() + (*auau_js).Pt()));
      auau_match_aj[key]->Fill(*auau_cent,
                            fabs((*auau_jlm).Pt() - (*auau_jsm).Pt())/((*auau_jlm).Pt() + (*auau_jsm).Pt()));
      
      // auau dphi
      double dphi = (*auau_jl).Phi() - (*auau_js).Phi();
      
      //rotate dphi to be in [0,2*pi]
      while (dphi < 0)
      dphi += 2 * TMath::Pi();
      while (dphi > 2.0 * TMath::Pi())
      dphi -= 2.0 * TMath::Pi();
      
      auau_dphi[key]->Fill(*auau_cent, dphi);
      
      // and grefmult/npart
      auau_npart[key]->Fill(*auau_cent, *auau_nprt);
    }
    
    
    // extract nPart in centrality bins
    auau_npart_cent[key] = SplitByCentralityNormalized(auau_npart[key]);
    
    auau_hard_lead_pt_cent[key] = SplitByCentralityNormalized(auau_hard_lead_pt[key]);
    auau_hard_sub_pt_cent[key] = SplitByCentralityNormalized(auau_hard_sub_pt[key]);
    auau_match_lead_pt_cent[key] = SplitByCentralityNormalized(auau_match_lead_pt[key]);
    auau_match_sub_pt_cent[key] = SplitByCentralityNormalized(auau_match_sub_pt[key]);
    auau_hard_aj_cent[key] = SplitByCentralityNormalized(auau_hard_aj[key]);
    auau_match_aj_cent[key] = SplitByCentralityNormalized(auau_match_aj[key]);
    auau_off_axis_lead_pt_cent[key] = SplitByCentralityNormalized(auau_off_axis_lead_pt[key]);
    auau_off_axis_sub_pt_cent[key] = SplitByCentralityNormalized(auau_off_axis_sub_pt[key]);
    auau_off_axis_aj_cent[key] = SplitByCentralityNormalized(auau_off_axis_aj[key]);
    
    auau_hard_lead_const_cent[key] = SplitByCentralityNormalized(auau_hard_lead_const[key]);
    auau_hard_sub_const_cent[key] = SplitByCentralityNormalized(auau_hard_sub_const[key]);
    auau_match_lead_const_cent[key] = SplitByCentralityNormalized(auau_match_lead_const[key]);
    auau_match_sub_const_cent[key] = SplitByCentralityNormalized(auau_match_sub_const[key]);
    auau_off_axis_lead_const_cent[key] = SplitByCentralityNormalized(auau_off_axis_lead_const[key]);
    auau_off_axis_sub_const_cent[key] = SplitByCentralityNormalized(auau_off_axis_sub_const[key]);
    
    auau_hard_lead_rho_cent[key] = SplitByCentralityNormalized(auau_hard_lead_rho[key]);
    auau_hard_sub_rho_cent[key] = SplitByCentralityNormalized(auau_hard_sub_rho[key]);
    auau_match_lead_rho_cent[key] = SplitByCentralityNormalized(auau_match_lead_rho[key]);
    auau_match_sub_rho_cent[key] = SplitByCentralityNormalized(auau_match_sub_rho[key]);
    
    auau_hard_lead_sig_cent[key] = SplitByCentralityNormalized(auau_hard_lead_sig[key]);
    auau_hard_sub_sig_cent[key] = SplitByCentralityNormalized(auau_hard_sub_sig[key]);
    auau_match_lead_sig_cent[key] = SplitByCentralityNormalized(auau_match_lead_sig[key]);
    auau_match_sub_sig_cent[key] = SplitByCentralityNormalized(auau_match_sub_sig[key]);
    
    auau_off_axis_lead_rho_cent[key] = SplitByCentralityNormalized(auau_off_axis_lead_rho[key]);
    auau_off_axis_sub_rho_cent[key] = SplitByCentralityNormalized(auau_off_axis_sub_rho[key]);
    auau_off_axis_lead_sig_cent[key] = SplitByCentralityNormalized(auau_off_axis_lead_sig[key]);
    auau_off_axis_sub_sig_cent[key] = SplitByCentralityNormalized(auau_off_axis_sub_sig[key]);
    
    for (int i = 0; i < auau_hard_aj_cent[key].size(); ++i) {
      auau_hard_aj_cent[key][i]->RebinX(2);
      auau_match_aj_cent[key][i]->RebinX(2);
      auau_off_axis_aj_cent[key][i]->RebinX(2);
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
    
    Overlay1D(auau_hard_lead_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_lead_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_hard_sub_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_sub_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_match_lead_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_lead_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_match_sub_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_sub_rho",
              "", "#rho", "fraction");
    
    Overlay1D(auau_hard_lead_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_lead_sig",
              "", "#sigma", "fraction");
    Overlay1D(auau_hard_sub_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_hard_sub_sig",
              "", "#sigma", "fraction");
    Overlay1D(auau_match_lead_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_lead_sig",
              "", "#sigma", "fraction");
    Overlay1D(auau_match_sub_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_match_sub_sig",
              "", "#sigma", "fraction");
    
    Overlay1D(auau_off_axis_lead_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_lead_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_off_axis_sub_rho_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_sub_sub_rho",
              "", "#rho", "fraction");
    Overlay1D(auau_off_axis_lead_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_lead_sig",
              "", "#sigma", "fraction");
    Overlay1D(auau_off_axis_sub_sig_cent[key], refcent_string, hopts, copts, out_loc, "auau_off_axis_sub_sig",
              "", "#sigma", "fraction");
    
  }
  
  return 0;
}
