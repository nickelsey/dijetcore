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
#include "TLatex.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include "dijetcore/lib/map.h"
#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(pythia, "", "input file for Pythia data");
DIJETCORE_DEFINE_string(outputDir, "results", "directory for output");
DIJETCORE_DEFINE_bool(setScanning, false, "fixes the initial radius, and scans through matched radii")
DIJETCORE_DEFINE_double(initRadius, 0.2, "initial radius when set to scan");
DIJETCORE_DEFINE_string(radii, "0.2,0.25,0.3,0.35,0.4", "radii to put in grid");
DIJETCORE_DEFINE_string(constPt, "1.0,1.5,2.0,2.5,3.0", "radii to put in grid");

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

int main(int argc, char* argv[]) {

  string usage = "Pythia differential di-jet imbalance print routine";

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
  if (!boost::filesystem::exists(FLAGS_pythia)) {
    std::cout << "input file " << FLAGS_pythia;
    std::cout << "doesn't exist: exiting" << std::endl;
    return 1;
  }

  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);

  // read in the file
  TFile pythia_file(FLAGS_pythia.c_str(), "READ");
  
  // and the radii & constituent pT we want to use
  std::vector<double> radii = dijetcore::ParseArgStringToVec<double>(FLAGS_radii);
  std::sort(radii.begin(), radii.end());
  std::vector<string> radii_string;
  for (auto& val : radii) {
    std::stringstream stream;
    stream << val;
    radii_string.push_back(stream.str());
  }
  std::vector<double> constpt = dijetcore::ParseArgStringToVec<double>(FLAGS_constPt);
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
  }
  else {
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
  std::unordered_map<string, TTree*> pythia_trees;

  GetTreesFromFile(pythia_file, pythia_trees);
  // match keys
  for (auto entry : pythia_trees) {
    keys.push_back(entry.first);
    parsed_keys[entry.first] = dijetcore::ParseStringToDijetKey(entry.first);
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
        if (FLAGS_setScanning &&
            key.lead_init_r == FLAGS_initRadius &&
            key.sub_init_r == FLAGS_initRadius &&
            key.lead_match_r == radii[rad] &&
            key.sub_match_r == radii[rad] &&
            key.lead_init_const_pt == constpt[pt] &&
            key.sub_init_const_pt == constpt[pt]) {
          grid_keys[rad].push_back(keys[key_index]);
          grid_key_params[rad].push_back(key);
        }
        else if (!FLAGS_setScanning &&
                 key.lead_init_r == radii[rad] &&
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

  std::unordered_map<string, TH1D*> hard_lead_eta;
  std::unordered_map<string, TH1D*> hard_lead_phi;
  std::unordered_map<string, TH1D*> hard_lead_pt;
  std::unordered_map<string, TH1D*> hard_lead_const;
  std::unordered_map<string, TH1D*> hard_lead_rho;
  std::unordered_map<string, TH1D*> hard_lead_sig;
  std::unordered_map<string, TH1D*> match_lead_eta;
  std::unordered_map<string, TH1D*> match_lead_phi;
  std::unordered_map<string, TH1D*> match_lead_pt;
  std::unordered_map<string, TH1D*> match_lead_const;
  std::unordered_map<string, TH1D*> match_lead_rho;
  std::unordered_map<string, TH1D*> match_lead_sig;

  std::unordered_map<string, TH1D*> hard_sub_eta;
  std::unordered_map<string, TH1D*> hard_sub_phi;
  std::unordered_map<string, TH1D*> hard_sub_pt;
  std::unordered_map<string, TH1D*> hard_sub_const;
  std::unordered_map<string, TH1D*> hard_sub_rho;
  std::unordered_map<string, TH1D*> hard_sub_sig;
  std::unordered_map<string, TH1D*> match_sub_eta;
  std::unordered_map<string, TH1D*> match_sub_phi;
  std::unordered_map<string, TH1D*> match_sub_pt;
  std::unordered_map<string, TH1D*> match_sub_const;
  std::unordered_map<string, TH1D*> match_sub_rho;
  std::unordered_map<string, TH1D*> match_sub_sig;

  std::unordered_map<string, TH1D*> npart;
  std::unordered_map<string, TH1D*> pthat;
  std::unordered_map<string, TH1D*> hard_dphi;
  std::unordered_map<string, TH1D*> match_dphi;
  std::unordered_map<string, TH1D*> lead_dr;
  std::unordered_map<string, TH1D*> sub_dr;
  std::unordered_map<string, TH1D*> lead_dpt;
  std::unordered_map<string, TH1D*> sub_dpt;
  std::unordered_map<string, TH1D*> lead_dpt_frac;
  std::unordered_map<string, TH1D*> sub_dpt_frac;
  
  std::unordered_map<string, TH1D*> hard_aj;
  std::unordered_map<string, TH1D*> match_aj;
  std::unordered_map<string, TH1D*> hard_aj_test;
  std::unordered_map<string, TH1D*> match_aj_test;
  
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
    
    //increment entry counter
    entry++;
    LOG(INFO) << "loading: " << key;
    // make our output directory
    string out_loc = FLAGS_outputDir + "/" + key;
    boost::filesystem::path dir(out_loc.c_str());
    boost::filesystem::create_directories(dir);
    
    // build prefix for names so histograms don't get confused
    // (root fuckin sucks)
    std::string hist_prefix = "key_" + std::to_string(entry) + "_";
    
    // load our reader
    TTree* current_tree = pythia_trees[key];
    TTreeReader reader(current_tree);
    
    // create readervalues for auau first
    dijetcore::unique_ptr<TTreeReaderValue<int>> nprt = dijetcore::make_unique<TTreeReaderValue<int>>(reader, "npart");
    dijetcore::unique_ptr<TTreeReaderValue<double>> ptht = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "pthat");
    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jl = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "jl");
    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> js = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "js");
    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jlm = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "jlm");
    dijetcore::unique_ptr<TTreeReaderValue<TLorentzVector>> jsm = dijetcore::make_unique<TTreeReaderValue<TLorentzVector>>(reader, "jsm");
    dijetcore::unique_ptr<TTreeReaderValue<int>> jlconst = dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jlconst");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jlrho = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jlrho");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jlsig = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jlsig");
    dijetcore::unique_ptr<TTreeReaderValue<int>> jlmconst = dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jlmconst");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jlmrho = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jlmrho");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jlmsig = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jlmsig");
    dijetcore::unique_ptr<TTreeReaderValue<int>> jsconst = dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jsconst");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jsrho = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jsrho");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jssig = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jssig");
    dijetcore::unique_ptr<TTreeReaderValue<int>> jsmconst = dijetcore::make_unique<TTreeReaderValue<int>>(reader, "jsmconst");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jsmrho = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jsmrho");
    dijetcore::unique_ptr<TTreeReaderValue<double>> jsmsig = dijetcore::make_unique<TTreeReaderValue<double>>(reader, "jsmsig");
    
    // create the histograms
    hard_lead_eta[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardleadeta").c_str(), "",
                                              50, -1, 1);
    hard_lead_phi[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardleadphi").c_str(), "",
                                              50, -TMath::Pi(), TMath::Pi());
    hard_lead_pt[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardleadpt").c_str(), "",
                                             50, 0, 50);
    hard_lead_const[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardleadconst").c_str(), "",
                                                50, 0, 100);
    hard_lead_rho[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardleadrho").c_str(), "",
                                              50, 0, 100);
    hard_lead_sig[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardleadsig").c_str(), "",
                                              50, 0, 20);
    
    match_lead_eta[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchleadeta").c_str(), "",
                                               50, -1, 1);
    match_lead_phi[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchleadphi").c_str(), "",
                                               50, -TMath::Pi(), TMath::Pi());
    match_lead_pt[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchleadpt").c_str(), "",
                                              50, 0, 50);
    match_lead_const[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchleadconst").c_str(), "",
                                                 50, 0, 100);
    match_lead_rho[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchleadrho").c_str(), "",
                                               50, 0, 100);
    match_lead_sig[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchleadsig").c_str(), "",
                                               50, 0, 20);
    
    hard_sub_eta[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardsubeta").c_str(), "",
                                             50, -1, 1);
    hard_sub_phi[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardsubphi").c_str(), "",
                                             50, -TMath::Pi(), TMath::Pi());
    hard_sub_pt[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardsubpt").c_str(), "",
                                            50, 0, 50);
    hard_sub_const[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardsubconst").c_str(), "",
                                               50, 0, 100);
    hard_sub_rho[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardsubrho").c_str(), "",
                                             50, 0, 100);
    hard_sub_sig[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardsubsig").c_str(), "",
                                             50, 0, 20);
    
    match_sub_eta[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchsubeta").c_str(), "",
                                              50, -1, 1);
    match_sub_phi[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchsubphi").c_str(), "",
                                              50, -TMath::Pi(), TMath::Pi());
    match_sub_pt[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchsubpt").c_str(), "",
                                             50, 0, 50);
    match_sub_const[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchsubconst").c_str(), "",
                                                50, 0, 100);
    match_sub_rho[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchsubrho").c_str(), "",
                                              50, 0, 100);
    match_sub_sig[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchsubsig").c_str(), "",
                                              50, 0, 20);
    
    npart[key] = new TH1D(dijetcore::MakeString(hist_prefix, "npart").c_str(), "",
                                      100, 0, 1200);
    pthat[key] = new TH1D(dijetcore::MakeString(hist_prefix, "pthat").c_str(), "",
                          100, 0, 100);
    hard_dphi[key] = new TH1D(dijetcore::MakeString(hist_prefix, "harddphi").c_str(), "",
                                          50, 0, TMath::Pi());
    match_dphi[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchdphi").c_str(), "",
                                           50, 0, TMath::Pi());
    lead_dr[key] = new TH1D(dijetcore::MakeString(hist_prefix, "leaddr").c_str(), "",
                                        50, 0, 0.5);
    sub_dr[key] = new TH1D(dijetcore::MakeString(hist_prefix, "subdr").c_str(), "",
                                       50, 0, 0.5);
    lead_dpt[key] = new TH1D(dijetcore::MakeString(hist_prefix, "leaddpt").c_str(), "",
                                         50, -20, 20);
    sub_dpt[key] = new TH1D(dijetcore::MakeString(hist_prefix, "subdptfrac").c_str(), "",
                                        50, -20, 20);
    lead_dpt_frac[key] = new TH1D(dijetcore::MakeString(hist_prefix, "leaddptfrac").c_str(), "",
                                              50, -2.0, 2.0);
    sub_dpt_frac[key] = new TH1D(dijetcore::MakeString(hist_prefix, "subdpt").c_str(), "",
                                             50, -2.0, 2.0);
    
    hard_aj[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardaj").c_str(), "",
                                        15, 0, 0.9);
    match_aj[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchaj").c_str(), "",
                                         15, 0, 0.9);
    hard_aj_test[key] = new TH1D(dijetcore::MakeString(hist_prefix, "hardajtest").c_str(), "",
                                             10000, 0, 0.9);
    match_aj_test[key] = new TH1D(dijetcore::MakeString(hist_prefix, "matchajtest").c_str(), "",
                                              10000, 0, 0.9);
    
    // event loop
    while (reader.Next()) {
      
      
      // check jet eta
      if (FLAGS_setScanning &&
          (fabs((*jl)->Eta()) > max_eta ||
           fabs((*js)->Eta()) > max_eta))
        continue;
      
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
      
      npart[key]->Fill(**nprt);
      pthat[key]->Fill(**ptht);
      hard_dphi[key]->Fill(fabs((*jl)->DeltaPhi(**js)));
      match_dphi[key]->Fill(fabs((*jlm)->DeltaPhi(**jsm)));
      lead_dr[key]->Fill((*jl)->DeltaR(**jlm));
      sub_dr[key]->Fill((*js)->DeltaR(**jsm));
      lead_dpt[key]->Fill((*jl)->Pt() - (*jlm)->Pt());
      sub_dpt[key]->Fill((*js)->Pt() - (*jsm)->Pt());
      lead_dpt_frac[key]->Fill(((*jl)->Pt() - (*jlm)->Pt()) / (*jl)->Pt());
      sub_dpt_frac[key]->Fill(((*js)->Pt() - (*jsm)->Pt()) / (*js)->Pt());
      
      hard_aj[key]->Fill(fabs((*jl)->Pt() - (*js)->Pt()) / ((*jl)->Pt() + (*js)->Pt()));
      match_aj[key]->Fill(fabs((*jlm)->Pt() - (*jsm)->Pt()) / ((*jlm)->Pt() + (*jsm)->Pt()));
      hard_aj_test[key]->Fill(fabs((*jl)->Pt() - (*js)->Pt()) / ((*jl)->Pt() + (*js)->Pt()));
      match_aj_test[key]->Fill(fabs((*jlm)->Pt() - (*jsm)->Pt()) / ((*jlm)->Pt() + (*jsm)->Pt()));
      
    }
    
    // now we will compare all the other observables
    PrettyPrint1D(npart[key], hopts, copts, "pythia",  out_loc, "npart", "", "N_{part}", "fraction");
    PrettyPrint1D(pthat[key], hopts, copts, "pythia",  out_loc, "pthat", "", "N_{part}", "fraction");
    PrettyPrint1D(hard_dphi[key], hopts, copts, "pythia",  out_loc, "hard_dphi", "", "d#phi", "fraction");
    PrettyPrint1D(match_dphi[key], hopts, copts, "pythia",  out_loc, "match_dphi", "", "d#phi", "fraction");
    PrettyPrint1D(lead_dr[key], hopts, copts, "pythia",  out_loc, "lead_dr", "", "d#phi", "fraction");
    PrettyPrint1D(sub_dr[key], hopts, copts, "pythia",  out_loc, "sub_dr", "", "d#phi", "fraction");
    PrettyPrint1D(lead_dpt[key], hopts, copts, "pythia",  out_loc, "lead_dpt", "", "dp_{T}", "fraction");
    PrettyPrint1D(sub_dpt[key], hopts, copts, "pythia",  out_loc, "sub_dpt", "", "dp_{T}", "fraction");
    PrettyPrint1D(lead_dpt_frac[key], hopts, copts, "pythia",  out_loc, "lead_dpt_frac", "", "dp_{T}", "fraction");
    PrettyPrint1D(sub_dpt_frac[key], hopts, copts, "pythia",  out_loc, "sub_dpt_frac", "", "dp_{T}", "fraction");
    PrettyPrint1D(hard_lead_eta[key], hopts, copts, "pythia",  out_loc, "hard_lead_eta", "", "#eta", "fraction");
    PrettyPrint1D(match_lead_eta[key], hopts, copts, "pythia",  out_loc, "match_lead_eta", "", "#eta", "fraction");
    PrettyPrint1D(hard_sub_eta[key], hopts, copts, "pythia",  out_loc, "hard_sub_eta", "", "#eta", "fraction");
    PrettyPrint1D(match_sub_eta[key], hopts, copts, "pythia",  out_loc, "match_sub_eta", "", "#eta", "fraction");
    PrettyPrint1D(hard_lead_phi[key], hopts, copts, "pythia",  out_loc, "hard_lead_phi", "", "#phi", "fraction");
    PrettyPrint1D(match_lead_phi[key], hopts, copts, "pythia",  out_loc, "match_lead_phi", "", "#phi", "fraction");
    PrettyPrint1D(hard_sub_phi[key], hopts, copts, "pythia",  out_loc, "hard_sub_phi", "", "#phi", "fraction");
    PrettyPrint1D(match_sub_phi[key], hopts, copts, "pythia",  out_loc, "match_sub_phi", "", "#phi", "fraction");
    PrettyPrint1D(hard_lead_pt[key], hopts, copts, "pythia",  out_loc, "hard_lead_pt", "", "p_{T}", "fraction");
    PrettyPrint1D(match_lead_pt[key], hopts, copts, "pythia",  out_loc, "match_lead_pt", "", "p_{T}", "fraction");
    PrettyPrint1D(hard_sub_pt[key], hopts, copts, "pythia",  out_loc, "hard_sub_pt", "", "p_{T}", "fraction");
    PrettyPrint1D(match_sub_pt[key], hopts, copts, "pythia",  out_loc, "match_sub_pt", "", "p_{T}", "fraction");
    PrettyPrint1D(hard_lead_const[key], hopts, copts, "pythia",  out_loc, "hard_lead_const", "", "N_{part}", "fraction");
    PrettyPrint1D(match_lead_const[key], hopts, copts, "pythia",  out_loc, "match_lead_const", "", "N_{part}", "fraction");
    PrettyPrint1D(hard_sub_const[key], hopts, copts, "pythia",  out_loc, "hard_sub_const", "", "N_{part}", "fraction");
    PrettyPrint1D(match_sub_const[key], hopts, copts, "pythia",  out_loc, "match_sub_const", "", "N_{part}", "fraction");
    PrettyPrint1D(hard_lead_rho[key], hopts, copts, "pythia",  out_loc, "hard_lead_rho", "", "rho", "fraction");
    PrettyPrint1D(match_lead_rho[key], hopts, copts, "pythia",  out_loc, "match_lead_rho", "", "rho", "fraction");
    PrettyPrint1D(hard_sub_rho[key], hopts, copts, "pythia",  out_loc, "hard_sub_rho", "", "rho", "fraction");
    PrettyPrint1D(match_sub_rho[key], hopts, copts, "pythia",  out_loc, "match_sub_rho", "", "rho", "fraction");
    PrettyPrint1D(hard_lead_sig[key], hopts, copts, "pythia",  out_loc, "hard_lead_sig", "", "sigma", "fraction");
    PrettyPrint1D(match_lead_sig[key], hopts, copts, "pythia",  out_loc, "match_lead_sig", "", "sigma", "fraction");
    PrettyPrint1D(hard_sub_sig[key], hopts, copts, "pythia",  out_loc, "hard_sub_sig", "", "sigma", "fraction");
    PrettyPrint1D(match_sub_sig[key], hopts, copts, "pythia",  out_loc, "match_sub_sig", "", "sigma", "fraction");
    
  } // key

  // make a directory for our grid outputs
  string out_loc_grid = FLAGS_outputDir + "/grid";
  LOG(INFO) << "output directory: " << out_loc_grid;
  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path grid_dir(out_loc_grid.c_str());
  boost::filesystem::create_directories(grid_dir);
  TH2D* mean_hard_pt = new TH2D("mean_hard_pt", ";R;p_{T}^{const}",
                                radii.size(), - 0.5, radii.size() - 0.5,
                                constpt.size(), -0.5, constpt.size() - 0.5);
  TH2D* mean_match_pt = new TH2D("mean_match_pt", ";R;p_{T}^{const}",
                                 radii.size(), - 0.5, radii.size() - 0.5,
                                 constpt.size(), -0.5, constpt.size() - 0.5);
  TH2D* mean_pt_hat = new TH2D("mean_pt_hat", ";R;p_{T}^{const}",
                               radii.size(), - 0.5, radii.size() - 0.5,
                               constpt.size(), -0.5, constpt.size() - 0.5);
  for (int j = 0; j < radii_string.size(); ++j) {
    mean_hard_pt->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
    mean_match_pt->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
    mean_pt_hat->GetXaxis()->SetBinLabel(j+1, radii_string[j].c_str());
  }
  for (int j = 0; j < constpt_string.size(); ++j) {
    mean_hard_pt->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
    mean_match_pt->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
    mean_pt_hat->GetYaxis()->SetBinLabel(j+1, constpt_string[j].c_str());
  }
  
  for (int rad = 0; rad < radii.size(); ++rad) {
    for (int pt = 0; pt < constpt.size(); ++pt) {
      
      int x_bin = mean_hard_pt->GetXaxis()->FindBin(rad);
      int y_bin = mean_hard_pt->GetYaxis()->FindBin(pt);
      string key = grid_keys[rad][pt];
      dijetcore::DijetKey key_params = grid_key_params[rad][pt];
      
      TH1D* hist_hard = hard_lead_pt[key];
      mean_hard_pt->SetBinContent(x_bin, y_bin, hist_hard->GetMean());
      TH1D* hist_match = match_lead_pt[key];
      mean_match_pt->SetBinContent(x_bin, y_bin, hist_match->GetMean());
      TH1D* hist_pt = pthat[key];
      mean_pt_hat->SetBinContent(x_bin, y_bin, hist_pt->GetMean());
    }
  }
  Print2DSimple(mean_hard_pt, hopts, copts, out_loc_grid, "mean_pt_hard", "", "R", "p_{T}^{const}", "TEXT COLZ");
  Print2DSimple(mean_match_pt, hopts, copts, out_loc_grid, "mean_pt_match", "", "R", "p_{T}^{const}", "TEXT COLZ");
  Print2DSimple(mean_pt_hat, hopts, copts, out_loc_grid, "mean_pt_hat", "", "R", "p_{T}^{const}", "TEXT COLZ");
 
  return 0;
}
