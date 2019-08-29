#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/data/centrality/centrality_run7.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/trigger_lookup.h"
#include "dijetcore/util/data/vector_conversion.h"
#include "dijetcore/worker/dijet_worker/dijet_worker.h"
#include "dijetcore/worker/dijet_worker/off_axis_worker.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "dijetcore/lib/root_ext/TTreeReader.h"
#include "TTreeReaderValue.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoReader.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetVector.h"
#include "TStarJetVectorContainer.h"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "kolja_tree",       // root file containing Kolja's ResultTree
    "output_directory", // directory for output root file to written to
    "bad_tower_list",   // bad tower list for TStarJetPicoReader
    "bad_run_list",     // bad run list for TStarJetPicoReader
    "triggers",         // trigger id set (defined in
                        // dijetcore/util/data/trigger_lookup.h)
    "jet_radius",       // jet radius
    "hard_core_cut",    // jet hard core cut
    "lead_jet_pt",      // leading hard core jet pt cut
    "sublead_jet_pt"    // subleading hard core jet pt cut
};

fastjet::PseudoJet MatchJet(fastjet::PseudoJet &ref, double jet_radius,
                            std::vector<fastjet::PseudoJet> candidates) {
  fastjet::Selector circle = fastjet::SelectorCircle(jet_radius);
  circle.set_reference(ref);
  auto selected = fastjet::sorted_by_pt(circle(candidates));
  if (selected.size() > 0)
    return selected[0];
  return fastjet::PseudoJet();
}

string FormatJet(TLorentzVector &v) {
  return dijetcore::MakeString(std::fixed, std::setprecision(3), "pt: ", v.Pt(),
                               ", eta: ", v.Eta(), ", phi: ", v.Phi());
}

string FormatJet(fastjet::PseudoJet &v) {
  return dijetcore::MakeString(std::fixed, std::setprecision(3), "pt: ", v.pt(),
                               ", eta: ", v.eta(), ", phi: ", v.phi_std());
}

int main(int argc, char *argv[]) {
  string usage = "Simple y7 Aj analysis";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // parse configuration file
  dijetcore::json config;
  try {
    config = dijetcore::LoadConfig(FLAGS_config, required_params);
  } catch (std::exception &e) {
    LOG(ERROR) << "error loading config: " << e.what();
    return 1;
  }

  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(FLAGS_input)) {
    LOG(ERROR) << "input file does not exist: " << FLAGS_input;
    return 1;
  }

  // find output directory from configuration file
  string output_dir = config["output_directory"];

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(output_dir);
  boost::filesystem::create_directories(dir);

  // copy config file to output directory
  boost::filesystem::path input_file(FLAGS_config.c_str());
  boost::filesystem::path copy_path(dir);
  copy_path /= input_file.filename();
  boost::filesystem::copy_file(
      input_file, copy_path,
      boost::filesystem::copy_option::overwrite_if_exists);

  // read in kolja's tree
  std::string kolja_filename = config["kolja_tree"];
  TFile kolja_file(kolja_filename.c_str(), "READ");

  TTreeReader kolja_tree("ResultTree", &kolja_file);
  TTreeReaderValue<unsigned> k_run(kolja_tree, "runid");
  TTreeReaderValue<unsigned> k_event(kolja_tree, "eventid");
  TTreeReaderValue<double> k_refmult(kolja_tree, "refmult");
  TTreeReaderValue<TLorentzVector> k_jl(kolja_tree, "j1");
  TTreeReaderValue<TLorentzVector> k_js(kolja_tree, "j2");
  TTreeReaderValue<TLorentzVector> k_jlm(kolja_tree, "jm1");
  TTreeReaderValue<TLorentzVector> k_jsm(kolja_tree, "jm2");
  TTreeReaderValue<float> k_rho(kolja_tree, "rho");

  // create a map where the key is the event and run ID, and the stored value is
  // the tree index of the event, so that we can quickly access each event in
  // kolja's tree
  dijetcore::pair_map<unsigned, unsigned, unsigned> kolja_events;
  while (kolja_tree.Next()) {
    std::pair<unsigned, unsigned> key = {*k_run, *k_event};
    kolja_events[key] = kolja_tree.GetTree()->GetReadEntry();
  }
  // reset te tree
  kolja_tree.Restart();

  // create output file from the given directory, name & id
  string outfile_name =
      output_dir + "/" + FLAGS_name + dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // build our input chain
  TChain *chain = dijetcore::NewChainFromInput(FLAGS_input);

  // initialize the reader
  TStarJetPicoReader *reader = new TStarJetPicoReader();
  dijetcore::InitReaderWithDefaults(reader, chain, config["bad_tower_list"],
                                    config["bad_run_list"]);
  reader->GetEventCuts()->SetVertexZDiffCut(99999);

  // initialize Run 7 centrality
  dijetcore::CentralityRun7 centrality;

  // get the trigger IDs that will be used
  std::set<unsigned> triggers = dijetcore::GetTriggerIDs(config["triggers"]);
  LOG(INFO) << "taking triggers: " << config["triggers"] << " for analysis";
  LOG(INFO) << "trigger ids: " << triggers;

  // output tree
  TTree *output = new TTree("MyTree", "aj for auau");

  // saved output for TTree (event level)
  unsigned runid = 0;
  unsigned eventid = 0;
  unsigned gref = 0;
  unsigned ref = 0;
  unsigned cent = 0;
  double vz = 0.0;
  double rho = 0.0;
  TLorentzVector jl;
  TLorentzVector js;
  TLorentzVector jlm;
  TLorentzVector jsm;
  double jl_area;
  double js_area;
  double jlm_area;
  double jsm_area;
  unsigned jl_const;
  unsigned js_const;
  unsigned jlm_const;
  unsigned jsm_const;

  output->Branch("runid", &runid);
  output->Branch("eventid", &eventid);
  output->Branch("gref", &gref);
  output->Branch("ref", &ref);
  output->Branch("cent", &cent);
  output->Branch("rho", &rho);
  output->Branch("vz", &vz);
  output->Branch("jl", &jl);
  output->Branch("js", &js);
  output->Branch("jlm", &jlm);
  output->Branch("jsm", &jsm);
  output->Branch("jl_area", &jl_area);
  output->Branch("js_area", &js_area);
  output->Branch("jlm_area", &jlm_area);
  output->Branch("jsm_area", &jsm_area);
  output->Branch("jl_const", &jl_const);
  output->Branch("js_const", &js_const);
  output->Branch("jlm_const", &jlm_const);
  output->Branch("jsm_const", &jsm_const);

  // fastjet definitions
  double jet_radius = static_cast<double>(config["jet_radius"]);

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_radius);
  fastjet::JetDefinition bkg_jet_def(fastjet::kt_algorithm, jet_radius);
  fastjet::GhostedAreaSpec area_spec(1.0 + jet_radius, 1, 0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,
                                   area_spec);

  // define track and jet selectors
  fastjet::Selector track_pt_min_selector =
      fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsRapMax(1.0);
  fastjet::Selector hard_core_track_selector =
      fastjet::SelectorPtMin(config["hard_core_cut"]) &&
      fastjet::SelectorAbsRapMax(1.0);
  fastjet::Selector jet_selector =
      fastjet::SelectorPtMin(config["sublead_jet_pt"]) &&
      fastjet::SelectorAbsRapMax(1.0 - jet_radius);
  fastjet::Selector match_jet_selector =
      fastjet::SelectorPtMin(0.0) && fastjet::SelectorAbsRapMax(1.0);
  fastjet::Selector bkg_selector =
      fastjet::SelectorAbsRapMax(1.0 - jet_radius) *
      (!fastjet::SelectorNHardest(2));

  int dijets = 0;
  try {
    while (reader->NextEvent()) {
      // Print out reader status every 10 seconds
      reader->PrintStatus(10);

      // headers for convenience
      TStarJetPicoEventHeader *header = reader->GetEvent()->GetHeader();

      // check if event fired a trigger we will use
      if (triggers.size() != 0) {
        bool use_event = false;
        for (auto trigger : triggers) {
          if (header->HasTriggerId(trigger)) {
            use_event = true;
          }
        }
        if (!use_event)
          continue;
      }

      // get centrality
      unsigned centrality_bin = centrality.Centrality9(
          reader->GetEvent()->GetHeader()->GetGReferenceMultiplicity());

      // set event quantities for tree
      runid = header->GetRunId();
      eventid = header->GetEventId();
      ref = header->GetReferenceMultiplicity();
      gref = header->GetGReferenceMultiplicity();
      vz = header->GetPrimaryVertexZ();
      cent = centrality_bin;

      bool kolja_found = false;
      if (kolja_events.count({runid, eventid}) && gref >= 269) {
        LOG(INFO) << "";
        kolja_found = true;
        LOG(INFO) << "kolja had a dijet in this event: "
                  << std::pair<unsigned, unsigned>(runid, eventid);
        unsigned tree_idx = kolja_events[{runid, eventid}];
        kolja_tree.GetTree()->GetEntry(tree_idx);
        LOG(INFO) << "leading jet: " << FormatJet(*k_jl);
        LOG(INFO) << "subleading jet: " << FormatJet(*k_js);
      }

      // get the vector container
      TStarJetVectorContainer<TStarJetVector> *container =
          reader->GetOutputContainer();
      std::vector<fastjet::PseudoJet> primary_particles;
      dijetcore::ConvertTStarJetVector(container, primary_particles);

      // select tracks above the minimum pt threshold
      primary_particles = track_pt_min_selector(primary_particles);

      // select hard-core constituents
      std::vector<fastjet::PseudoJet> hard_particles =
          hard_core_track_selector(primary_particles);

      // do hard-core jetfinding
      fastjet::ClusterSequenceArea hc_cluster(hard_particles, jet_def,
                                              area_def);
      std::vector<fastjet::PseudoJet> hc_candidates =
          fastjet::sorted_by_pt(jet_selector(hc_cluster.inclusive_jets()));

      if (hc_candidates.size() < 2) {
        if (kolja_found) {
          if (hc_candidates.size() == 1) {
            LOG(INFO) << "MISSED: ----------------------------------";
            LOG(INFO) << "we did not find a jet pair to match, only one jet:";
            LOG(INFO) << "leading jet: " << FormatJet(hc_candidates[0]);
          } else {
            LOG(INFO) << "MISSED: ----------------------------------";
            LOG(INFO) << "we did not find a jet pair to match";
          }
        }
        continue;
      }

      if (hc_candidates[0].pt() < config["lead_jet_pt"]) {
        if (kolja_found) {
          LOG(INFO) << "MISSED: ----------------------------------";
          LOG(INFO) << "we found a pair but the leading jet's pT was below "
                       "threshold: ";
          LOG(INFO) << "leading jet: " << FormatJet(hc_candidates[0]);
        }
        continue;
      }

      double dphi = hc_candidates[0].delta_phi_to(hc_candidates[1]);

      if (fabs(dphi) < TMath::Pi() - 0.4) {
        if (kolja_found) {
          LOG(INFO) << "MISSED: ----------------------------------";
          LOG(INFO) << "found a pair but they weren't back to back: ";
          LOG(INFO) << "leading jet: " << FormatJet(hc_candidates[0]);
          LOG(INFO) << "subleading jet: " << FormatJet(hc_candidates[0]);
          LOG(INFO) << "dphi: "
                    << hc_candidates[0].delta_phi_to(hc_candidates[1]);
        }
        continue;
      }

      if (gref >= 269)
        dijets++;

      fastjet::PseudoJet lead = hc_candidates[0];
      fastjet::PseudoJet sublead = hc_candidates[1];

      if (kolja_found) {
        LOG(INFO) << "we found jets in this event as well: ";
        LOG(INFO) << "leading jet: " << FormatJet(lead);
        LOG(INFO) << "subleading jet: " << FormatJet(sublead);
      } else if (gref >= 269) {
        LOG(INFO) << "ERROR -------------------------";
        LOG(INFO) << "We found jets that kolja didn't";
        LOG(INFO) << "leading jet: " << FormatJet(lead);
        LOG(INFO) << "subleading jet: " << FormatJet(sublead);
      }

      // fill this part of the tree
      jl = TLorentzVector(lead.px(), lead.py(), lead.pz(), lead.E());
      js =
          TLorentzVector(sublead.px(), sublead.py(), sublead.pz(), sublead.E());
      jl_area = lead.area();
      js_area = sublead.area();

      int lead_constituents = 0, sublead_constituents = 0;
      for (auto &c : lead.constituents())
        if (c.pt() > 0.2)
          lead_constituents++;

      for (auto c : sublead.constituents())
        if (c.pt() > 0.2)
          sublead_constituents++;

      jl_const = lead_constituents;
      js_const = sublead_constituents;

      // perform clustering for matched jets
      fastjet::ClusterSequenceArea match_cluster(primary_particles, jet_def,
                                                 area_def);
      std::vector<fastjet::PseudoJet> match_candidates = fastjet::sorted_by_pt(
          match_jet_selector(match_cluster.inclusive_jets()));

      // do background subtraction
      fastjet::JetMedianBackgroundEstimator match_bkg_est(
          bkg_selector, bkg_jet_def, area_def);
      match_bkg_est.set_particles(primary_particles);
      // Subtract A*rho from the original pT
      fastjet::Subtractor match_bkg_sub(&match_bkg_est);
      match_candidates = fastjet::sorted_by_pt(match_bkg_sub(match_candidates));

      // do matching to hard-core in dR
      fastjet::PseudoJet lead_match =
          MatchJet(lead, jet_radius, match_candidates);
      fastjet::PseudoJet sub_match =
          MatchJet(sublead, jet_radius, match_candidates);
      if (lead_match.pt() < 0.001)
        continue;
      if (sub_match.pt() < 0.001)
        continue;
      jlm = TLorentzVector(lead_match.px(), lead_match.py(), lead_match.pz(),
                           lead_match.E());
      jsm = TLorentzVector(sub_match.px(), sub_match.py(), sub_match.pz(),
                           sub_match.E());

      jlm_area = lead_match.area();
      jsm_area = sub_match.area();
      rho = match_bkg_est.rho();

      int lead_match_constituents = 0, sub_match_constituents = 0;
      for (auto &c : lead_match.constituents())
        if (c.pt() > 0.2)
          lead_match_constituents++;

      for (auto &c : sub_match.constituents())
        if (c.pt() > 0.2)
          sub_match_constituents++;

      jlm_const = lead_match_constituents;
      jsm_const = sub_match_constituents;

      output->Fill();
    }
  } catch (std::exception &e) {
    LOG(ERROR) << "Caught: " << e.what() << " during analysis loop.";
  }

  output->Write();
  out.Close();

  std::cout << "number of dijets: " << dijets << std::endl;

  gflags::ShutDownCommandLineFlags();
  return 0;
}