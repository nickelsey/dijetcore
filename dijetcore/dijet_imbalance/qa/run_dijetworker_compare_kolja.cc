#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

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
#include "TTreeReader.h"
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

double CompareVal(double val1, double val2, double tol) {
  if (fabs(val1 - val2) > tol)
    return false;
  return true;
}

std::vector<string> CompareJets(fastjet::PseudoJet &m, fastjet::PseudoJet &k) {
  std::vector<string> ret;
  if (!CompareVal(m.pt(), k.pt(), 0.01)) {
    ret.push_back(dijetcore::MakeString("my pt: ", m.pt(), " k pt: ", k.pt()));
  }
  if (!CompareVal(m.eta(), k.eta(), 0.01)) {
    ret.push_back(
        dijetcore::MakeString("my eta: ", m.eta(), " k eta: ", k.eta()));
  }
  if (!CompareVal(m.phi(), k.phi(), 0.01)) {
    ret.push_back(
        dijetcore::MakeString("my phi: ", m.phi(), " k phi: ", k.phi()));
  }
  return ret;
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

  // create worker
  dijetcore::DijetWorker worker(fastjet::antikt_algorithm, 20.0, 0.4, 0.4, 10.0,
                                0.4, 0.4, 2.0, 0.2, 2.0, 0.2, 1.0);
  worker.forceConstituentPtEquality(true);
  worker.forceConstituentEtaEquality(true);
  worker.forceJetResolutionEquality(true);
  worker.forceMatchJetResolutionEquality(true);

  fastjet::Selector lead_selector =
      fastjet::SelectorPtMin(config["lead_jet_pt"]) &&
      fastjet::SelectorAbsEtaMax(1.0 - jet_radius);
  fastjet::Selector sublead_selector =
      fastjet::SelectorPtMin(config["sublead_jet_pt"]) &&
      fastjet::SelectorAbsEtaMax(1.0 - jet_radius);

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

      if (gref < 269)
        continue;

      bool kolja_found = false;
      if (kolja_events.count({runid, eventid}) && gref >= 269) {
        kolja_found = true;
      }

      // get the vector container
      TStarJetVectorContainer<TStarJetVector> *container =
          reader->GetOutputContainer();
      std::vector<fastjet::PseudoJet> primary_particles;
      dijetcore::ConvertTStarJetVector(container, primary_particles);

      // run the worker
      auto &worker_out = worker.run(primary_particles);

      for (auto &result : worker_out) {
        std::string key = result.first;
        dijetcore::ClusterOutput &out = *result.second.get();
        if (out.foundDijet())
          dijets++;
        LOG(INFO) << "event: " << std::pair<unsigned, unsigned>(runid, eventid);
        if (kolja_found) {
          kolja_tree.GetTree()->GetEntry(kolja_events[{runid, eventid}]);
        }
        // case where kolja and I found di-jets
        if (kolja_found && out.found_lead && out.found_sublead) {
          fastjet::PseudoJet k_vl;
          k_vl.reset_momentum_PtYPhiM((*k_jl).Pt(), (*k_jl).Eta(),
                                      (*k_jl).Phi(), 0.0);
          fastjet::PseudoJet k_vs;
          k_vs.reset_momentum_PtYPhiM((*k_js).Pt(), (*k_js).Eta(),
                                      (*k_js).Phi(), 0.0);

          std::vector<string> lead_dif = CompareJets(out.lead_hard, k_vl);
          std::vector<string> sub_dif = CompareJets(out.sublead_hard, k_vs);

          if (lead_dif.size()) {
            LOG(INFO) << "difference in lead jets:";
            for (auto &e : lead_dif)
              LOG(INFO) << e;
            LOG(INFO) << "rho? : " << out.lead_hard_rho;
          }
          if (sub_dif.size()) {
            LOG(INFO) << "difference in sublead jets:";
            for (auto &e : sub_dif)
              LOG(INFO) << e;
            LOG(INFO) << "rho? : " << out.lead_hard_rho;
          }
        }

        // case where I found a di-jet pair and he didn't
        if (out.found_lead && out.found_sublead && !kolja_found) {
          LOG(INFO) << "We found di-jets that kolja did not";
          LOG(INFO) << "lead jet: " << FormatJet(out.lead_hard);
          LOG(INFO) << "sublead jet: " << FormatJet(out.sublead_hard);
          LOG(INFO) << "dphi: " << out.lead_hard.delta_phi_to(out.sublead_hard);
        }

        // case where He found a di-jet pair and I did not
        if (kolja_found && (!out.found_lead || !out.found_sublead)) {
          LOG(INFO) << "Kolja found jets that we didn't";
          LOG(INFO) << "lead jet: " << FormatJet(*k_jl);
          LOG(INFO) << "sublead jet: " << FormatJet(*k_js);

          // find most probable pair
          auto lead_cl = out.lead_hard_seq;
          auto sublead_cl = out.sublead_hard_seq;

          auto lead_jets =
              fastjet::sorted_by_pt(lead_selector(lead_cl->inclusive_jets()));
          auto sublead_jets = fastjet::sorted_by_pt(
              sublead_selector(sublead_cl->inclusive_jets()));

          if (lead_jets.size() == 0) {
            LOG(INFO) << "we don't have a lead jet";
            if (sublead_jets.size()) {
              LOG(INFO) << "highest pt jet: " << FormatJet(sublead_jets[0]);
            }
            continue;
          }
          auto lead_jet = lead_jets[0];
          double dphi = 0.0;
          fastjet::PseudoJet sublead_jet;
          for (auto & j : sublead_jets) {
            dphi = fabs(lead_jet.delta_phi_to(j));
            if (dphi > TMath::Pi() - 0.4) {
              sublead_jet = j;
              break;
            }
          }

          if (sublead_jet.pt() < config["sublead_jet_pt"]) {
            LOG(INFO) << "we don't have a sublead jet";
            continue;
          }

          LOG(INFO) << "my jets that are closest: ";
          LOG(INFO) << "lead: " << FormatJet(lead_jet);
          LOG(INFO) << "sublead: " << FormatJet(sublead_jet);
          LOG(INFO) << "dphi: " << dphi;
          LOG(INFO) << "rho? : " << out.lead_hard_rho;
        }
      }
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