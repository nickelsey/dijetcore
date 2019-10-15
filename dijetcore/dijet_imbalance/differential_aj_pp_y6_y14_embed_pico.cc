#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/data/centrality/centrality_run14.h"
#include "dijetcore/util/data/efficiency/run14_eff.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/trigger_lookup.h"
#include "dijetcore/util/data/vector_conversion.h"
#include "dijetcore/worker/dijet_worker/dijet_worker.h"

#include "jetreader/reader/centrality.h"
#include "jetreader/reader/reader.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TTree.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoReader.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetVector.h"
#include "TStarJetVectorContainer.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(embedInput, "", "input file for embedding Au+Au");
DIJETCORE_DEFINE_int(nEmbed, 1, "number of times to reuse a pp event");
DIJETCORE_DEFINE_string(efficiencyFile, "",
                        "file with efficiency curves if applicable");
DIJETCORE_DEFINE_string(outputDir, "tmp", "directory for output");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(towList,
                        "resources/bad_tower_lists/y14_y6_bad_tower.txt",
                        "bad tower list");
DIJETCORE_DEFINE_string(runList, "resources/bad_run_lists/y14_bad_run.txt",
                        "bad run list");
DIJETCORE_DEFINE_string(triggers, "y6ppht", "trigger selection");
DIJETCORE_DEFINE_string(embedTriggers, "y14vpdmb30",
                        "trigger selection for embedding");
DIJETCORE_DEFINE_string(constEta, "1.0",
                        "constitutent eta cuts (comma separated)");
DIJETCORE_DEFINE_string(leadConstPt, "2.0",
                        "leading constituent pT cut (comma separated)");
DIJETCORE_DEFINE_string(subConstPt, "2.0",
                        "subleading constituent pT cut (comma separated)");
DIJETCORE_DEFINE_string(leadConstPtMatch, "0.2",
                        "matched leading constituent pT cut (comma separated)");
DIJETCORE_DEFINE_string(
    subConstPtMatch, "0.2",
    "matched subleading constituent pT cut (comma separated)");
DIJETCORE_DEFINE_string(leadR, "0.4", "leading jet R (comma separated");
DIJETCORE_DEFINE_string(leadMatchR, "0.4",
                        "leading matched jet R (comma separated");
DIJETCORE_DEFINE_string(subR, "0.4", "subleading jet R (comma separated");
DIJETCORE_DEFINE_string(subMatchR, "0.4",
                        "subleading matched jet R (comma separated");
DIJETCORE_DEFINE_string(leadJetPt, "20.0",
                        "leading jet pT cut (comma separated)");
DIJETCORE_DEFINE_string(subJetPt, "10.0",
                        "subleading jet pT cut (comma separated)");
DIJETCORE_DEFINE_int(towerUnc, 0, "tower scaling for systematic uncertainties");
DIJETCORE_DEFINE_int(
    trackingUnc, 0,
    "tracking efficiency ratio scaling for systematic uncertainties");
DIJETCORE_DEFINE_int(seed, 0,
                     "seed for the random number generator, so that results "
                     "can be reproducible");
DIJETCORE_DEFINE_bool(useEfficiency, true,
                      "flag to scale pp efficiency to be run 7 AuAu-like");
DIJETCORE_DEFINE_int(minCentrality, 3,
                     "minimum cutoff for embedding centrality - should not be "
                     "greater than 2 if efficiency is on")
    DIJETCORE_DEFINE_bool(forceConstituentPtEquality, true,
                          "Only use DijetDefinitions where pT const is equal "
                          "in leading/subleading jets");
DIJETCORE_DEFINE_bool(forceConstituentEtaEquality, true,
                      "Only use DijetDefinitions where eta const is equal in "
                      "leading/subleading jets");
DIJETCORE_DEFINE_bool(
    forceJetResolutionEquality, true,
    "Only use DijetDefinitions where leading/subleading R are equivalent");
DIJETCORE_DEFINE_bool(
    forceMatchJetResolutionEquality, false,
    "Only use DijetDefinitions where initial and matched R are equivalent");

bool GetEmbedEvent(jetreader::Reader &reader);

void ConvertPPToPseudojets(TStarJetPicoReader *reader,
                           jetreader::Reader &embed_reader,
                           dijetcore::Run14Eff *efficiency, double tower_scale,
                           std::vector<fastjet::PseudoJet> &particles,
                           std::mt19937 &gen,
                           std::uniform_real_distribution<> &dis,
                           TProfile2D* eff_ratio);

std::array<fastjet::PseudoJet, 4>
ClusterPP(std::vector<fastjet::PseudoJet> &input,
          const dijetcore::DijetDefinition &def,
          const fastjet::PseudoJet &lead_jet,
          const fastjet::PseudoJet &sublead_jet);

int main(int argc, char *argv[]) {

  string usage = "Run 6 embedded in Run 14 pp differential di-jet imbalance "
                 "analysis routine";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(FLAGS_input)) {
    LOG(ERROR) << "input file does not exist: " << FLAGS_input;
    return 1;
  }

  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);

  // first, build our input chain
  TChain *chain = dijetcore::NewChainFromInput(FLAGS_input);

  // initialize the reader
  TStarJetPicoReader *reader = new TStarJetPicoReader();
  dijetcore::InitReaderWithDefaults(reader, chain, FLAGS_towList);

  // if requested, setup the embedding reader
  if (FLAGS_embedInput.empty()) {
    LOG(ERROR) << "No embedding file specified - exiting";
    return 1;
  }

  // create the jetreader for the embedding event
  jetreader::Reader embed_reader(FLAGS_embedInput);
  embed_reader.centrality().loadCentralityDef(jetreader::CentDefId::Run14);
  if (!FLAGS_runList.empty())
    embed_reader.eventSelector()->addBadRuns(FLAGS_runList);
  embed_reader.eventSelector()->setVzRange(-30, 30);
  embed_reader.eventSelector()->setdVzMax(3.0);
  if (!FLAGS_towList.empty())
    embed_reader.towerSelector()->addBadTowers(FLAGS_towList);
  embed_reader.towerSelector()->setEtMax(30.0);
  embed_reader.towerSelector()->setEtMin(0.2);
  embed_reader.trackSelector()->setDcaMax(3.0);
  embed_reader.trackSelector()->setPtMin(0.2);
  embed_reader.trackSelector()->setPtMax(30.0);
  embed_reader.trackSelector()->setNHitsMin(15);
  embed_reader.trackSelector()->setNHitsFracMin(0.52);

  // get the trigger IDs that will be used for the pp data
  std::set<unsigned> triggers;
  if (!FLAGS_triggers.empty()) {
    triggers = dijetcore::GetTriggerIDs(FLAGS_triggers);
    LOG(INFO) << "taking triggers: " << FLAGS_triggers << " for primary";
    LOG(INFO) << "trigger ids: " << triggers;
  }

  // and load the trigger IDs for the embedding event
  std::set<unsigned> embed_triggers;
  if (!FLAGS_embedTriggers.empty()) {
    embed_triggers = dijetcore::GetTriggerIDs(FLAGS_embedTriggers);
    LOG(INFO) << "taking triggers: " << FLAGS_embedTriggers << " for embedding";
    LOG(INFO) << "trigger ids: " << embed_triggers;
    embed_reader.eventSelector()->addTriggerIds(FLAGS_embedTriggers);
  }
  // initialize reader
  embed_reader.init();

  // initialize efficiency curves
  dijetcore::Run14Eff *efficiency;
  if (FLAGS_efficiencyFile.empty())
    efficiency = new dijetcore::Run14Eff();
  else
    efficiency = new dijetcore::Run14Eff(FLAGS_efficiencyFile);

  switch (FLAGS_trackingUnc) {
  case 0:
    efficiency->setSystematicUncertainty(dijetcore::TrackingUnc::NONE);
    break;
  case 1:
    efficiency->setSystematicUncertainty(dijetcore::TrackingUnc::POSITIVE);
    break;
  case -1:
    efficiency->setSystematicUncertainty(dijetcore::TrackingUnc::NEGATIVE);
    break;
  default:
    LOG(ERROR) << "undefined tracking efficiency setting, exiting";
    return 1;
  }

  // define our tower uncertainty scaling as well
  const double tower_scale_percent = 0.02;
  if (abs(FLAGS_towerUnc) > 1) {
    LOG(ERROR) << "undefined tower efficiency setting, exiting";
    return 1;
  }
  double tower_scale = 1.0 + tower_scale_percent * FLAGS_towerUnc;

  // and we'll need a random number generator for randomly throwing away tracks
  std::mt19937 gen(FLAGS_seed);
  std::uniform_real_distribution<> dis(0.0, 1.0);

  // define a selector to accept tracks with pT > 0.2 GeV, inside our nominal
  // eta acceptance of + / -1.0
  fastjet::Selector track_pt_min_selector =
      fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);

  // parse jetfinding variables
  // --------------------------

  // first, hard code the algorithm to be anti-kt
  std::set<fastjet::JetAlgorithm> alg{fastjet::antikt_algorithm};

  // constituent range
  std::set<double> const_eta =
      dijetcore::ParseArgString<double>(FLAGS_constEta);

  // leading jet
  std::set<double> lead_const_hard_pt =
      dijetcore::ParseArgString<double>(FLAGS_leadConstPt);
  std::set<double> lead_const_match_pt =
      dijetcore::ParseArgString<double>(FLAGS_leadConstPtMatch);
  std::set<double> lead_R = dijetcore::ParseArgString<double>(FLAGS_leadR);
  std::set<double> lead_R_match =
      dijetcore::ParseArgString<double>(FLAGS_leadMatchR);
  std::set<double> lead_hard_pt =
      dijetcore::ParseArgString<double>(FLAGS_leadJetPt);

  // subleading jet
  std::set<double> sublead_const_hard_pt =
      dijetcore::ParseArgString<double>(FLAGS_subConstPt);
  std::set<double> sublead_const_match_pt =
      dijetcore::ParseArgString<double>(FLAGS_subConstPtMatch);
  std::set<double> sublead_R = dijetcore::ParseArgString<double>(FLAGS_subR);
  std::set<double> sublead_R_match =
      dijetcore::ParseArgString<double>(FLAGS_subMatchR);
  std::set<double> sublead_hard_pt =
      dijetcore::ParseArgString<double>(FLAGS_subJetPt);

  LOG(INFO) << "grid variables ";
  LOG(INFO) << "constituent eta: " << const_eta;
  LOG(INFO) << "lead constituent hard pt cut: " << lead_const_hard_pt;
  LOG(INFO) << "lead constituent match pt cut: " << lead_const_match_pt;
  LOG(INFO) << "lead hard R: " << lead_R;
  LOG(INFO) << "lead match R: " << lead_R_match;
  LOG(INFO) << "lead hard jet pt: " << lead_hard_pt;
  LOG(INFO) << "sublead constituent hard pt cut: " << sublead_const_hard_pt;
  LOG(INFO) << "sublead constituent match pt cut: " << sublead_const_match_pt;
  LOG(INFO) << "sublead hard R: " << sublead_R;
  LOG(INFO) << "sublead match R: " << sublead_R_match;
  LOG(INFO) << "sublead hard jet pt: " << sublead_hard_pt;

  // here we can initialize the worker
  LOG(INFO) << "initializing worker...";
  dijetcore::DijetWorker worker(
      alg, lead_hard_pt, lead_R, lead_R_match, sublead_hard_pt, sublead_R,
      sublead_R_match, lead_const_hard_pt, lead_const_match_pt,
      sublead_const_hard_pt, sublead_const_match_pt, const_eta);
  worker.forceConstituentPtEquality(FLAGS_forceConstituentPtEquality);
  worker.forceConstituentEtaEquality(FLAGS_forceConstituentEtaEquality);
  worker.forceJetResolutionEquality(FLAGS_forceJetResolutionEquality);
  worker.forceMatchJetResolutionEquality(FLAGS_forceMatchJetResolutionEquality);
  worker.initialize();

  LOG(INFO) << "worker initialized - number of dijet definitions: "
            << worker.size();

  std::set<std::string> keys = worker.keys();
  for (auto key : keys)
    LOG(INFO) << key;

  // create output file from the given directory, name & id
  string outfile_name = FLAGS_outputDir + "/" + FLAGS_name +
                        dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // create an output tree for each definition
  // -----------------------------------------

  std::unordered_map<std::string, dijetcore::unique_ptr<TTree>> trees;

  std::unordered_map<std::string, int> run_id_dict;
  std::unordered_map<std::string, int> event_id_dict;
  std::unordered_map<std::string, double> vz_dict;
  std::unordered_map<std::string, int> refmult_dict;
  std::unordered_map<std::string, int> grefmult_dict;
  std::unordered_map<std::string, double> refmultcorr_dict;
  std::unordered_map<std::string, double> grefmultcorr_dict;
  std::unordered_map<std::string, int> cent_dict;
  std::unordered_map<std::string, double> zdcrate_dict;
  std::unordered_map<std::string, double> reactionplane_dict;
  std::unordered_map<std::string, int> nglobal_dict;
  std::unordered_map<std::string, int> npart_dict;
  std::unordered_map<std::string, TLorentzVector> lead_hard_jet_dict;
  std::unordered_map<std::string, int> lead_hard_jet_nconst_dict;
  std::unordered_map<std::string, double> lead_hard_rho_dict;
  std::unordered_map<std::string, double> lead_hard_sigma_dict;
  std::unordered_map<std::string, double> lead_hard_area_dict;
  std::unordered_map<std::string, TH1D *> lead_hard_jet_const_pt_dict;
  std::unordered_map<std::string, TH1D *> lead_hard_jet_const_dr_dict;
  std::unordered_map<std::string, TLorentzVector> lead_match_jet_dict;
  std::unordered_map<std::string, int> lead_match_jet_nconst_dict;
  std::unordered_map<std::string, double> lead_match_rho_dict;
  std::unordered_map<std::string, double> lead_match_sigma_dict;
  std::unordered_map<std::string, double> lead_match_area_dict;
  std::unordered_map<std::string, TH1D *> lead_match_jet_const_pt_dict;
  std::unordered_map<std::string, TH1D *> lead_match_jet_const_dr_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_hard_jet_dict;
  std::unordered_map<std::string, int> sublead_hard_jet_nconst_dict;
  std::unordered_map<std::string, double> sublead_hard_rho_dict;
  std::unordered_map<std::string, double> sublead_hard_sigma_dict;
  std::unordered_map<std::string, double> sublead_hard_area_dict;
  std::unordered_map<std::string, TH1D *> sublead_hard_jet_const_pt_dict;
  std::unordered_map<std::string, TH1D *> sublead_hard_jet_const_dr_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_match_jet_dict;
  std::unordered_map<std::string, int> sublead_match_jet_nconst_dict;
  std::unordered_map<std::string, double> sublead_match_rho_dict;
  std::unordered_map<std::string, double> sublead_match_sigma_dict;
  std::unordered_map<std::string, double> sublead_match_area_dict;
  std::unordered_map<std::string, TH1D *> sublead_match_jet_const_pt_dict;
  std::unordered_map<std::string, TH1D *> sublead_match_jet_const_dr_dict;
  // for embedding
  std::unordered_map<std::string, int> embed_runid_dict;
  std::unordered_map<std::string, int> embed_eventid_dict;
  std::unordered_map<std::string, int> embed_refmult_dict;
  std::unordered_map<std::string, int> embed_grefmult_dict;
  std::unordered_map<std::string, double> embed_refmultcorr_dict;
  std::unordered_map<std::string, double> embed_grefmultcorr_dict;
  std::unordered_map<std::string, int> embed_cent_dict;
  std::unordered_map<std::string, int> embed_npart_dict;
  std::unordered_map<std::string, double> embed_vz_dict;
  std::unordered_map<std::string, double> embed_zdc_dict;
  std::unordered_map<std::string, double> embed_rp_dict;

  // pp w/o embedding
  std::unordered_map<std::string, bool> found_jet_pp_dict;
  std::unordered_map<std::string, TLorentzVector> lead_hard_jet_pp_dict;
  std::unordered_map<std::string, TLorentzVector> lead_match_jet_pp_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_hard_jet_pp_dict;
  std::unordered_map<std::string, TLorentzVector> sublead_match_jet_pp_dict;

  // fill the maps first, so that they don't decide to resize/move themselves
  // after branch creation...
  for (auto key : keys) {
    run_id_dict.insert({key, 0});
    event_id_dict.insert({key, 0});
    vz_dict.insert({key, 0});
    refmult_dict.insert({key, 0});
    grefmult_dict.insert({key, 0});
    refmultcorr_dict.insert({key, 0});
    grefmultcorr_dict.insert({key, 0});
    cent_dict.insert({key, 0});
    zdcrate_dict.insert({key, 0});
    reactionplane_dict.insert({key, 0});
    nglobal_dict.insert({key, 0});
    npart_dict.insert({key, 0});
    lead_hard_jet_dict.insert({key, TLorentzVector()});
    lead_hard_jet_const_pt_dict.insert(
        {key, new TH1D(dijetcore::MakeString(key, "jlconstpt").c_str(), "", 90,
                       0, 30)});
    lead_hard_jet_const_dr_dict.insert(
        {key, new TH1D(dijetcore::MakeString(key, "jlconstdr").c_str(), "", 100,
                       0.0, 1.0)});
    lead_hard_jet_nconst_dict.insert({key, 0});
    lead_hard_rho_dict.insert({key, 0});
    lead_hard_sigma_dict.insert({key, 0});
    lead_hard_area_dict.insert({key, 0});
    lead_match_jet_dict.insert({key, TLorentzVector()});
    lead_match_jet_const_pt_dict.insert(
        {key, new TH1D(dijetcore::MakeString(key, "jlmconstpt").c_str(), "", 90,
                       0, 30)});
    lead_match_jet_const_dr_dict.insert(
        {key, new TH1D(dijetcore::MakeString(key, "jlmconstdr").c_str(), "",
                       100, 0.0, 1.0)});
    lead_match_jet_nconst_dict.insert({key, 0});
    lead_match_rho_dict.insert({key, 0});
    lead_match_sigma_dict.insert({key, 0});
    lead_match_area_dict.insert({key, 0});
    sublead_hard_jet_dict.insert({key, TLorentzVector()});
    sublead_hard_jet_const_pt_dict.insert(
        {key, new TH1D(dijetcore::MakeString(key, "jsconstpt").c_str(), "", 90,
                       0, 30)});
    sublead_hard_jet_const_dr_dict.insert(
        {key, new TH1D(dijetcore::MakeString(key, "jsconstdr").c_str(), "", 100,
                       0.0, 1.0)});
    sublead_hard_jet_nconst_dict.insert({key, 0});
    sublead_hard_rho_dict.insert({key, 0});
    sublead_hard_sigma_dict.insert({key, 0});
    sublead_hard_area_dict.insert({key, 0});
    sublead_match_jet_dict.insert({key, TLorentzVector()});
    sublead_match_jet_const_pt_dict.insert(
        {key, new TH1D(dijetcore::MakeString(key, "jsmconstpt").c_str(), "", 90,
                       0, 30)});
    sublead_match_jet_const_dr_dict.insert(
        {key, new TH1D(dijetcore::MakeString(key, "jsmconstdr").c_str(), "",
                       100, 0.0, 1.0)});
    sublead_match_jet_nconst_dict.insert({key, 0});
    sublead_match_rho_dict.insert({key, 0});
    sublead_match_sigma_dict.insert({key, 0});
    sublead_match_area_dict.insert({key, 0});

    embed_runid_dict.insert({key, 0});
    embed_eventid_dict.insert({key, 0});
    embed_vz_dict.insert({key, 0});
    embed_zdc_dict.insert({key, 0});
    embed_rp_dict.insert({key, 0});
    embed_cent_dict.insert({key, 0});
    embed_refmult_dict.insert({key, 0});
    embed_grefmult_dict.insert({key, 0});
    embed_refmultcorr_dict.insert({key, 0});
    embed_grefmultcorr_dict.insert({key, 0});
    embed_npart_dict.insert({key, 0});

    found_jet_pp_dict.insert({key, false});
    lead_hard_jet_pp_dict.insert({key, TLorentzVector()});
    lead_match_jet_pp_dict.insert({key, TLorentzVector()});
    sublead_hard_jet_pp_dict.insert({key, TLorentzVector()});
    sublead_match_jet_pp_dict.insert({key, TLorentzVector()});
  }

  for (auto key : keys) {
    dijetcore::unique_ptr<TTree> tmp =
        dijetcore::make_unique<TTree>(key.c_str(), key.c_str());
    // create branches for the tree
    tmp->Branch("runid", &run_id_dict[key]);
    tmp->Branch("eventid", &event_id_dict[key]);
    tmp->Branch("vz", &vz_dict[key]);
    tmp->Branch("refmult", &refmult_dict[key]);
    tmp->Branch("grefmult", &grefmult_dict[key]);
    tmp->Branch("refmultcorr", &refmultcorr_dict[key]);
    tmp->Branch("grefmultcorr", &grefmultcorr_dict[key]);
    tmp->Branch("cent", &cent_dict[key]);
    tmp->Branch("zdcrate", &zdcrate_dict[key]);
    tmp->Branch("rp", &reactionplane_dict[key]);
    tmp->Branch("nglobal", &nglobal_dict[key]);
    tmp->Branch("npart", &npart_dict[key]);
    tmp->Branch("jl", &lead_hard_jet_dict[key]);
    tmp->Branch("js", &sublead_hard_jet_dict[key]);
    tmp->Branch("jlm", &lead_match_jet_dict[key]);
    tmp->Branch("jsm", &sublead_match_jet_dict[key]);
    tmp->Branch("jlconst", &lead_hard_jet_nconst_dict[key]);
    tmp->Branch("jlrho", &lead_hard_rho_dict[key]);
    tmp->Branch("jlsig", &lead_hard_sigma_dict[key]);
    tmp->Branch("jlarea", &lead_hard_area_dict[key]);
    tmp->Branch("jlmconst", &lead_match_jet_nconst_dict[key]);
    tmp->Branch("jlmrho", &lead_match_rho_dict[key]);
    tmp->Branch("jlmsig", &lead_match_sigma_dict[key]);
    tmp->Branch("jlmarea", &lead_match_area_dict[key]);
    tmp->Branch("jsconst", &sublead_hard_jet_nconst_dict[key]);
    tmp->Branch("jsrho", &sublead_hard_rho_dict[key]);
    tmp->Branch("jssig", &sublead_hard_sigma_dict[key]);
    tmp->Branch("jsarea", &sublead_hard_area_dict[key]);
    tmp->Branch("jsmconst", &sublead_match_jet_nconst_dict[key]);
    tmp->Branch("jsmrho", &sublead_match_rho_dict[key]);
    tmp->Branch("jsmsig", &sublead_match_sigma_dict[key]);
    tmp->Branch("jsmarea", &sublead_match_area_dict[key]);

    tmp->Branch("embed_eventid", &embed_eventid_dict[key]);
    tmp->Branch("embed_runid", &embed_runid_dict[key]);
    tmp->Branch("embed_refmult", &embed_refmult_dict[key]);
    tmp->Branch("embed_grefmult", &embed_grefmult_dict[key]);
    tmp->Branch("embed_refmultcorr", &embed_refmultcorr_dict[key]);
    tmp->Branch("embed_grefmultcorr", &embed_grefmultcorr_dict[key]);
    tmp->Branch("embed_cent", &embed_cent_dict[key]);
    tmp->Branch("embed_npart", &embed_npart_dict[key]);
    tmp->Branch("embed_rp", &embed_rp_dict[key]);
    tmp->Branch("embed_zdcrate", &embed_zdc_dict[key]);
    tmp->Branch("embed_vz", &embed_vz_dict[key]);

    tmp->Branch("foundpp", &found_jet_pp_dict[key]);
    tmp->Branch("ppjl", &lead_hard_jet_pp_dict[key]);
    tmp->Branch("ppjs", &sublead_hard_jet_pp_dict[key]);
    tmp->Branch("ppjlm", &lead_match_jet_pp_dict[key]);
    tmp->Branch("ppjsm", &sublead_match_jet_pp_dict[key]);

    trees[key] = std::move(tmp);
    trees[key]->SetDirectory(0);
  }

  // histograms
  // ----------

  std::unordered_map<string, dijetcore::unique_ptr<TH1D>> lead_jet_count_dict;
  std::unordered_map<string, dijetcore::unique_ptr<TH1D>>
      sublead_jet_count_dict;

  dijetcore::unique_ptr<TProfile2D> eff_ratio =
      dijetcore::make_unique<TProfile2D>(
          "pp_eff_ratio", "average pp efficiency ratio;p_{T}; #eta; ratio ",
          200, 0, 5, 10, -1.0, 1.0);
  eff_ratio->SetDirectory(0);

  for (auto key : keys) {
    // create a unique histogram name for each key
    string lead_name = key + "_lead_count";
    string sublead_name = key + "_sublead_count";
    dijetcore::unique_ptr<TH1D> lead_tmp = dijetcore::make_unique<TH1D>(
        lead_name.c_str(), "count lead jets", 800, 0.5, 800.5);
    dijetcore::unique_ptr<TH1D> sublead_tmp = dijetcore::make_unique<TH1D>(
        sublead_name.c_str(), "count sublead jets", 800, 0.5, 800.5);

    lead_jet_count_dict[key] = std::move(lead_tmp);
    sublead_jet_count_dict[key] = std::move(sublead_tmp);
    lead_jet_count_dict[key]->SetDirectory(0);
    sublead_jet_count_dict[key]->SetDirectory(0);
  }

  // start the analysis loop
  // -----------------------
  try {
    while (reader->NextEvent()) {

      // Print out reader status every 10 seconds
      reader->PrintStatus(10);

      // headers for convenience
      TStarJetPicoEventHeader *header = reader->GetEvent()->GetHeader();

      // check if event fired a trigger we will use
      if (triggers.size() != 0) {
        bool use_event = false;
        for (auto trigger : triggers)
          if (header->HasTriggerId(trigger))
            use_event = true;
        if (!use_event)
          continue;
      }

      // we're using this p+p event - start the embedding loop
      for (int emb = 0; emb < FLAGS_nEmbed; ++emb) {

        if (!GetEmbedEvent(embed_reader)) {
          if (!GetEmbedEvent(embed_reader)) {
            LOG(ERROR) << "No embedding events satisfy criteria: exiting";
            return 1;
          }
        }

        int refmult = header->GetReferenceMultiplicity();
        double refmultcorr = refmult;
        int centrality_bin = embed_reader.centrality16();
        int embed_centrality = centrality_bin;

        std::vector<fastjet::PseudoJet> particles;
        std::vector<fastjet::PseudoJet> embed_particles;
        std::vector<fastjet::PseudoJet> primary_particles;

        // if successful, load the embedding event
        for (auto &j : embed_reader.pseudojets()) {
          particles.push_back(j);
          embed_particles.push_back(j);
        }

        // and now convert the pp - if there are any efficiency curves to
        // apply, do it now
        ConvertPPToPseudojets(reader, embed_reader, efficiency, tower_scale,
                              primary_particles, gen, dis, eff_ratio.get());

        // add pp tracks to the full event
        particles.insert(particles.end(), primary_particles.begin(),
                         primary_particles.end());

        // select tracks above the minimum pt threshold
        particles = track_pt_min_selector(particles);

        // run the worker
        auto &worker_out = worker.run(particles);

        // process any found di-jet pairs
        for (auto &result : worker_out) {
          std::string key = result.first;
          dijetcore::ClusterOutput &out = *result.second.get();

          if (out.found_lead)
            lead_jet_count_dict[key]->Fill(header->GetReferenceMultiplicity());
          if (out.found_sublead)
            sublead_jet_count_dict[key]->Fill(
                header->GetReferenceMultiplicity());

          // now fill dijet results
          if (out.found_match) {
            // fill all branches for that key
            run_id_dict[key] = header->GetRunId();
            event_id_dict[key] = header->GetEventId();
            vz_dict[key] = header->GetPrimaryVertexZ();
            refmult_dict[key] = header->GetReferenceMultiplicity();
            grefmult_dict[key] = header->GetGReferenceMultiplicity();
            refmultcorr_dict[key] = refmultcorr;
            grefmultcorr_dict[key] =
                header->GetCorrectedGReferenceMultiplicity();
            cent_dict[key] = centrality_bin;
            zdcrate_dict[key] = header->GetZdcCoincidenceRate();
            reactionplane_dict[key] = header->GetReactionPlaneAngle();
            nglobal_dict[key] = header->GetNGlobalTracks();
            npart_dict[key] = primary_particles.size();

            embed_runid_dict[key] = embed_reader.picoDst()->event()->runId();
            embed_eventid_dict[key] = embed_reader.picoDst()->event()->eventId();
            embed_refmult_dict[key] = embed_reader.picoDst()->event()->refMult();
            embed_grefmult_dict[key] =
                embed_reader.picoDst()->event()->grefMult();
            embed_refmultcorr_dict[key] =
                embed_reader.centrality().refMultCorr();
            embed_grefmultcorr_dict[key] =
                embed_reader.picoDst()->event()->grefMult();
            embed_cent_dict[key] = embed_centrality;
            embed_vz_dict[key] = embed_reader.picoDst()->event()->primaryVertex().Z();
            embed_zdc_dict[key] = embed_reader.picoDst()->event()->ZDCx();
            embed_rp_dict[key] = 0.0;
            embed_npart_dict[key] = embed_particles.size();

            // set the four jets
            // set the four jets
            lead_hard_jet_dict[key] =
                TLorentzVector(out.lead_hard.px(), out.lead_hard.py(),
                               out.lead_hard.pz(), out.lead_hard.E());
            int l_const = 0;
            for (int i = 0; i < out.lead_hard.constituents().size(); ++i) {
              if (out.lead_hard.constituents()[i].pt() > 0.01) {
                l_const++;
                lead_hard_jet_const_pt_dict[key]->Fill(
                    out.lead_hard.constituents()[i].pt());
                lead_hard_jet_const_dr_dict[key]->Fill(
                    out.lead_hard.constituents()[i].delta_R(out.lead_hard));
              }
            }
            lead_hard_jet_nconst_dict[key] = l_const;
            lead_hard_rho_dict[key] = out.lead_hard_rho;
            lead_hard_sigma_dict[key] = out.lead_hard_sigma;
            lead_hard_area_dict[key] = out.lead_hard.area();
            lead_match_jet_dict[key] =
                TLorentzVector(out.lead_match.px(), out.lead_match.py(),
                               out.lead_match.pz(), out.lead_match.E());
            int lm_const = 0;
            for (int i = 0; i < out.lead_match.constituents().size(); ++i) {
              if (out.lead_match.constituents()[i].pt() > 0.01) {
                lm_const++;
                lead_match_jet_const_pt_dict[key]->Fill(
                    out.lead_match.constituents()[i].pt());
                lead_match_jet_const_dr_dict[key]->Fill(
                    out.lead_match.constituents()[i].delta_R(out.lead_match));
              }
            }
            lead_match_jet_nconst_dict[key] = lm_const;
            lead_match_rho_dict[key] = out.lead_match_rho;
            lead_match_sigma_dict[key] = out.lead_match_sigma;
            lead_match_area_dict[key] = out.lead_match.area();
            sublead_hard_jet_dict[key] =
                TLorentzVector(out.sublead_hard.px(), out.sublead_hard.py(),
                               out.sublead_hard.pz(), out.sublead_hard.E());
            int s_const = 0;
            for (int i = 0; i < out.sublead_hard.constituents().size(); ++i) {
              if (out.sublead_hard.constituents()[i].pt() > 0.01) {
                s_const++;
                sublead_hard_jet_const_pt_dict[key]->Fill(
                    out.sublead_hard.constituents()[i].pt());
                sublead_hard_jet_const_dr_dict[key]->Fill(
                    out.sublead_hard.constituents()[i].delta_R(
                        out.sublead_hard));
              }
            }
            sublead_hard_jet_nconst_dict[key] = s_const;
            sublead_hard_rho_dict[key] = out.sublead_hard_rho;
            sublead_hard_sigma_dict[key] = out.sublead_hard_sigma;
            sublead_hard_area_dict[key] = out.sublead_hard.area();
            sublead_match_jet_dict[key] =
                TLorentzVector(out.sublead_match.px(), out.sublead_match.py(),
                               out.sublead_match.pz(), out.sublead_match.E());
            int sm_const = 0;
            for (int i = 0; i < out.sublead_match.constituents().size(); ++i) {
              if (out.sublead_match.constituents()[i].pt() > 0.01) {
                sm_const++;
                sublead_match_jet_const_pt_dict[key]->Fill(
                    out.sublead_match.constituents()[i].pt());
                sublead_match_jet_const_dr_dict[key]->Fill(
                    out.sublead_match.constituents()[i].delta_R(
                        out.sublead_match));
              }
            }
            sublead_match_jet_nconst_dict[key] = sm_const;
            sublead_match_rho_dict[key] = out.sublead_match_rho;
            sublead_match_sigma_dict[key] = out.sublead_match_sigma;
            sublead_match_area_dict[key] = out.sublead_match.area();

            // now run clustering on only the pp
            auto pp_result = ClusterPP(primary_particles, *out.dijet_def,
                                       out.lead_hard, out.sublead_hard);
            if (pp_result[0].pt() > 0) {
              found_jet_pp_dict[key] = true;
              lead_hard_jet_pp_dict[key] =
                  TLorentzVector(pp_result[0].px(), pp_result[0].py(),
                                 pp_result[0].pz(), pp_result[0].E());
              lead_match_jet_pp_dict[key] =
                  TLorentzVector(pp_result[2].px(), pp_result[2].py(),
                                 pp_result[2].pz(), pp_result[2].E());
              sublead_hard_jet_pp_dict[key] =
                  TLorentzVector(pp_result[1].px(), pp_result[1].py(),
                                 pp_result[1].pz(), pp_result[1].E());
              sublead_match_jet_pp_dict[key] =
                  TLorentzVector(pp_result[3].px(), pp_result[3].py(),
                                 pp_result[3].pz(), pp_result[3].E());
            } else {
              found_jet_pp_dict[key] = false;
              lead_hard_jet_pp_dict[key] = TLorentzVector();
              lead_match_jet_pp_dict[key] = TLorentzVector();
              sublead_hard_jet_pp_dict[key] = TLorentzVector();
              sublead_match_jet_pp_dict[key] = TLorentzVector();
            }

            trees[key]->Fill();
          }
        }
      }
    }
  } catch (std::exception &e) {
    LOG(ERROR) << "Caught: " << e.what() << " during analysis loop.";
  }

  out.cd();
  for (auto &entry : trees) {
    entry.second->Write();
  }
  out.Close();

  return 0;
}

bool GetEmbedEvent(jetreader::Reader &reader) {
  while (reader.next()) {

    int centrality = reader.centrality16();

    if (centrality < 0)
      continue;
    if (centrality > FLAGS_minCentrality)
      continue;
    return true;
  }
  return false;
}

void ConvertPPToPseudojets(TStarJetPicoReader *reader,
                           jetreader::Reader &embed_reader,
                           dijetcore::Run14Eff *efficiency, double tower_scale,
                           std::vector<fastjet::PseudoJet> &particles,
                           std::mt19937 &gen,
                           std::uniform_real_distribution<> &dis,
                           TProfile2D* eff_ratio) {
  TStarJetVectorContainer<TStarJetVector> *container =
      reader->GetOutputContainer();
  for (int i = 0; i < container->GetEntries(); ++i) {
    TStarJetVector *sv = container->Get(i);

    if (fabs(sv->Eta()) > 1.0)
      continue;

    double scale = 1.0;

    if (sv->GetCharge() && FLAGS_useEfficiency) {
      // if the track is charged and we are using efficiencies,
      // then we get the efficiency ratio, and use that as the
      // probability to keep the track
      double ratio = efficiency->ratio(sv->Pt(), sv->Eta(), embed_reader.centrality16(),
                                       embed_reader.picoDst()->event()->ZDCx());
      if (!std::isfinite(ratio))
        continue;

      eff_ratio->Fill(sv->Pt(), sv->Eta(), ratio);

      double random_ = dis(gen);
      if (random_ > ratio)
        continue;
    } else {
      // the track is neutral - we keep it,
      // but we scale it by the tower_scale
      scale = tower_scale;
    }

    // now create the pseudojet, set the user index to the charge, and
    // scale
    fastjet::PseudoJet tmpPJ = fastjet::PseudoJet(*sv);
    tmpPJ *= scale;
    tmpPJ.set_user_index(sv->GetCharge());
    particles.push_back(tmpPJ);
  }
}

std::array<fastjet::PseudoJet, 4>
ClusterPP(std::vector<fastjet::PseudoJet> &input,
          const dijetcore::DijetDefinition &def,
          const fastjet::PseudoJet &lead_jet,
          const fastjet::PseudoJet &sublead_jet) {
  std::array<fastjet::PseudoJet, 4> ret;
  bool finish = true;

  // first, cluster with initial jet definitions (no area estimation necessary)
  fastjet::ClusterSequence seq_lead(
      def.lead->initialJetDef().constituentSelector()(input),
      def.lead->initialJetDef());
  fastjet::Selector circle_lead =
      fastjet::SelectorCircle(def.lead->initialJetDef().R());
  circle_lead.set_reference(lead_jet);
  std::vector<fastjet::PseudoJet> matched_lead =
      fastjet::sorted_by_pt(circle_lead(seq_lead.inclusive_jets()));
  if (matched_lead.size() != 0)
    ret[0] = matched_lead[0];
  else
    finish = false;

  // now do the same for subleading jet
  fastjet::ClusterSequence seq_sub(
      def.sub->initialJetDef().constituentSelector()(input),
      def.sub->initialJetDef());
  fastjet::Selector circle_sub =
      fastjet::SelectorCircle(def.sub->initialJetDef().R());
  fastjet::Selector recoil_selector =
      !fastjet::SelectorRectangle(2.1, TMath::Pi() - def.dPhi);
  circle_sub.set_reference(sublead_jet);
  fastjet::Selector sub_selector = circle_sub && recoil_selector;
  std::vector<fastjet::PseudoJet> matched_sub =
      fastjet::sorted_by_pt(circle_sub(seq_sub.inclusive_jets()));
  if (matched_sub.size() != 0)
    ret[1] = matched_sub[0];
  else
    finish = false;

  // if we have two dijets we will do matching as well - still, no background
  if (finish == false)
    return ret;

  fastjet::ClusterSequence seq_lead_match(
      def.lead->matchedJetDef().constituentSelector()(input),
      def.lead->matchedJetDef());
  fastjet::ClusterSequence seq_sub_match(
      def.lead->matchedJetDef().constituentSelector()(input),
      def.lead->matchedJetDef());
  fastjet::Selector circle_lead_match =
      fastjet::SelectorCircle(def.lead->matchedJetDef().R());
  circle_lead_match.set_reference(ret[0]);
  fastjet::Selector circle_sub_match =
      fastjet::SelectorCircle(def.sub->matchedJetDef().R());
  circle_sub_match.set_reference(ret[1]);

  std::vector<fastjet::PseudoJet> matched_lead_pp =
      fastjet::sorted_by_pt(circle_lead_match(seq_lead_match.inclusive_jets()));
  std::vector<fastjet::PseudoJet> matched_sub_pp =
      fastjet::sorted_by_pt(circle_sub_match(seq_sub_match.inclusive_jets()));
  if (matched_lead_pp.size() > 0 && matched_sub_pp.size() > 0) {
    ret[2] = matched_lead_pp[0];
    ret[3] = matched_sub_pp[0];
  }

  return ret;
}