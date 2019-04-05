#include <fstream>
#include <iostream>
#include <set>

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

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using dijetcore::dijetcore_map;
using dijetcore::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "output_directory",  // directory for output root file to written to
    "bad_tower_list",    // bad tower list for TStarJetPicoReader
    "bad_run_list",      // bad run list for TStarJetPicoReader
    "triggers",          // trigger id set (defined in
                         // dijetcore/util/data/trigger_lookup.h)
    "jet_radius",        // jet radius
    "hard_core_cut",     // jet hard core cut
    "lead_jet_pt",       // leading hard core jet pt cut
    "sublead_jet_pt"     // subleading hard core jet pt cut

};

fastjet::PseudoJet MatchJet(fastjet::PseudoJet& ref, double jet_radius, std::vector<fastjet::PseudoJet> candidates) {
  fastjet::Selector circle = fastjet::SelectorCircle(jet_radius);
  circle.set_reference(ref);
  auto selected = fastjet::sorted_by_pt(circle(candidates));
  if (selected.size() > 0)
    return selected[0];
  return fastjet::PseudoJet();
}

int main(int argc, char* argv[]) {
  string usage = "Differential di-jet imbalance secondary QA for analysis note";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // parse configuration file
  dijetcore::json config;
  try {
    config = dijetcore::LoadConfig(FLAGS_config, required_params);
  } catch (std::exception& e) {
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

  // create output file from the given directory, name & id
  string outfile_name =
      output_dir + "/" + FLAGS_name + dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // build our input chain
  TChain* chain = dijetcore::NewChainFromInput(FLAGS_input);

  // initialize the reader
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  dijetcore::InitReaderWithDefaults(reader, chain, config["bad_tower_list"],
                                    config["bad_run_list"]);

  // initialize Run 7 centrality
  dijetcore::CentralityRun7 centrality;

  // get the trigger IDs that will be used
  std::set<unsigned> triggers = dijetcore::GetTriggerIDs(config["triggers"]);
  LOG(INFO) << "taking triggers: " << config["triggers"] << " for analysis";
  LOG(INFO) << "trigger ids: " << triggers;

  // output tree
  TTree* output = new TTree("analysis_qa", "analysis_qa");

  // saved output for TTree (event level)
  unsigned runid, eventid, gref, ref, cent;
  double vx, vy, vz;

  output->Branch("runid", &runid);
  output->Branch("eventid", &eventid);
  output->Branch("gref", &gref);
  output->Branch("ref", &ref);
  output->Branch("cent", &cent);
  output->Branch("vz", &vz);
  output->Branch("vy", &vy);
  output->Branch("vx", &vx);

  // add booleans to count which events had which trigger
  dijetcore_map<unsigned, bool> trigger_flags;
  for (auto& trigger : triggers) trigger_flags[trigger] = false;

  for (auto& pair : trigger_flags)
    output->Branch(dijetcore::MakeString("t", pair.first).c_str(),
                   &pair.second);

  // jet info
  double leadpt, subpt, rho, sigma, leadarea, subarea, dphi;
  double leadptm, subptm, rhom, sigmam, leadaream, subaream, dphim;
  double leaddr, subdr;

  output->Branch("leadpt", &leadpt);
  output->Branch("subpt", &subpt);
  output->Branch("rho", &rho);
  output->Branch("sigma", &sigma);
  output->Branch("leadarea", &leadarea);
  output->Branch("subarea", &subarea);
  output->Branch("dphi", &dphi);
  output->Branch("leadptm", &leadptm);
  output->Branch("subptm", &subptm);
  output->Branch("rhom", &rhom);
  output->Branch("sigmam", &sigmam);
  output->Branch("leadaream", &leadaream);
  output->Branch("subaream", &subaream);
  output->Branch("dphim", &dphim);
  output->Branch("leaddr", &leaddr);
  output->Branch("subdr", &subdr);

  // fastjet definitions
  double jet_radius = static_cast<double>(config["jet_radius"]);

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_radius);
  fastjet::JetDefinition bkg_jet_def(fastjet::kt_algorithm, jet_radius);
  fastjet::GhostedAreaSpec area_spec(1.0 + jet_radius, 1, 0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,
                                   area_spec);

  // define track and jet selectors
  fastjet::Selector track_pt_min_selector =
      fastjet::SelectorPtMin(0.2) && fastjet::SelectorAbsEtaMax(1.0);
  fastjet::Selector hard_core_track_selector =
      fastjet::SelectorPtMin(config["hard_core_cut"]) &&
      fastjet::SelectorAbsEtaMax(1.0);
  fastjet::Selector jet_selector =
      fastjet::SelectorPtMin(config["sublead_jet_pt"]) &&
      fastjet::SelectorAbsEtaMax(1.0 - jet_radius);
  fastjet::Selector bkg_selector =
      fastjet::SelectorAbsEtaMax(1.0 - jet_radius) *
      (!fastjet::SelectorNHardest(2));

  try {
    while (reader->NextEvent()) {
      // Print out reader status every 10 seconds
      reader->PrintStatus(10);

      // headers for convenience
      TStarJetPicoEventHeader* header = reader->GetEvent()->GetHeader();

      // check if event fired a trigger we will use
      if (triggers.size() != 0) {
        bool use_event = false;
        for (auto trigger : triggers) {
          if (header->HasTriggerId(trigger)) {
            use_event = true;
            trigger_flags[trigger] = true;
          } else {
            trigger_flags[trigger] = false;
          }
        }
        if (!use_event) continue;
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
      vx = header->GetPrimaryVertexX();
      vy = header->GetPrimaryVertexY();
      cent = centrality_bin;

      // get the vector container
      TStarJetVectorContainer<TStarJetVector>* container =
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
          hc_cluster.inclusive_jets();

      // setup background subtraction
      fastjet::JetMedianBackgroundEstimator bkg_est(bkg_selector, bkg_jet_def,
                                                    area_def);
      bkg_est.set_particles(hard_particles);
      // Subtract A*rho from the original pT
      fastjet::Subtractor bkg_sub(&bkg_est);
      hc_candidates =
          fastjet::sorted_by_pt(jet_selector(bkg_sub(hc_candidates)));

      if (hc_candidates.size() < 2) {
        leadpt = 0.0;
        subpt = 0.0;
        rho = bkg_est.rho();
        sigma = bkg_est.sigma();
        leadarea = 0.0;
        subarea = 0.0;
        dphi = 0.0;

        if (hc_candidates.size() == 1) {
          leadpt = hc_candidates[0].pt();
          leadarea = hc_candidates[0].area();
        }
        output->Fill();
        continue;
      }

      leadpt = hc_candidates[0].pt();
      subpt = hc_candidates[1].pt();
      rho = bkg_est.rho();
      sigma = bkg_est.sigma();
      leadarea = hc_candidates[0].area();
      subarea = hc_candidates[1].area();
      dphi = hc_candidates[0].delta_phi_to(hc_candidates[1]);

      // perform clustering for matched jets
      fastjet::ClusterSequenceArea match_cluster(primary_particles, jet_def,
                                                 area_def);
      std::vector<fastjet::PseudoJet> match_candidates =
          match_cluster.inclusive_jets();

      // do background subtraction
      fastjet::JetMedianBackgroundEstimator match_bkg_est(
          bkg_selector, bkg_jet_def, area_def);
      match_bkg_est.set_particles(primary_particles);
      // Subtract A*rho from the original pT
      fastjet::Subtractor match_bkg_sub(&match_bkg_est);
      match_candidates = fastjet::sorted_by_pt(bkg_sub(match_candidates));

      // do matching to hard-core in dR
      fastjet::PseudoJet matched_lead = MatchJet(hc_candidates[0], jet_radius, match_candidates);
      fastjet::PseudoJet matched_sub = MatchJet(hc_candidates[1], jet_radius, match_candidates);

      double leadptm, subptm, rhom, sigmam, leadaream, subaream, dphim;
  double leaddr, subdr;
      if (matched_lead.pt() == 0.0) {
        leadptm = 0.0;
        leadaream = 0.0;
        leaddr = 0.0;
        dphim = 0.0;
      } 
      else {
        leadptm = matched_lead.pt();
        leadaream = matched_lead.area();
        leaddr = matched_lead.delta_R(hc_candidates[0]);
      }
      if (matched_sub.pt() == 0.0) {
        subptm = 0.0;
        subaream = 0.0;
        subdr = 0.0;
        dphim = 0.0;
      }
      else {
        subptm = matched_sub.pt();
        subaream = matched_sub.area();
        subdr = matched_sub.delta_R(hc_candidates[1]);
      }
      if (matched_lead.pt() && matched_sub.pt()) {
        dphim = matched_lead.delta_phi_to(matched_sub);
      }
      

      output->Fill();
    }
  } catch (std::exception& e) {
    LOG(ERROR) << "Caught: " << e.what() << " during analysis loop.";
  }

  output->Write();
  out.Close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}
