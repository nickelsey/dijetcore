// analyzes the correlation between the change in
// subjet splits when track loss occurs in p+p collisions

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/worker/mc/detector_sim_worker/pythia_star_sim.h"

#include <algorithm>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/IteratedSoftDrop.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/SoftDrop.hh"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "output_directory", // directory for output root file to written to
    "n_events",         // number of events to generate
    "pthat_min",        // minimum pT hat for generator
    "pthat_max",        // maximum pT hat for generator
    "track_pt_min",     // minimum pT for tracks at detector level
    "track_pt_max",     // maximum pT for tracks at detector level
    "track_eta_max",    // maximum eta for tracks at detector level
    "jet_pt_min",       // minimum pT for reconstructed jets at detector level
    "jet_radius",       // radius used during clustering (det and gen)
    "jet_alg",          // jetfinding algorithm (kt, antikt, or ca)
    "det_sim_mode",     // detector simulation method (none, gauss, star)
    "match_radius",     // the maximum radial distance to match a detector to a
                        // gen jet
    // if you want to add optional pythia style strings for the generator,
    // add them as "opt1"... "optn"
    "beta",           // softdrop angular parameter beta
    "theta_cut",      // softdrop angular cut
    "z_cut",          // softdrop z cut
    "angular_scaling" // angular scaling for recursive softdrop
};

double formation_time(double z, double pt, double theta) {
  return 2.0 / (z * (1.0 - z) * pt * pow(theta, 2.0));
}

double rel_err(double truth, double measured) {
  return (truth - measured) / truth;
}

struct splitInfo {
  fastjet::PseudoJet leading_prong;
  fastjet::PseudoJet subleading_prong;
  double zg;
  double rg;
  double pt;

  splitInfo() {
    leading_prong = fastjet::PseudoJet();
    subleading_prong = fastjet::PseudoJet();
    zg = -1.0;
    rg = -1.0;
    pt = -1.0;
  }
};

splitInfo parse_split(fastjet::PseudoJet &j) {

  if (!j.has_structure() || j.pieces().size() != 2 ||
      !j.has_structure_of<fastjet::contrib::RecursiveSoftDrop>())
    return splitInfo();

  std::vector<fastjet::PseudoJet> pieces = j.pieces();

  const fastjet::contrib::RecursiveSoftDrop::StructureType &j_struct =
      j.structure_of<fastjet::contrib::RecursiveSoftDrop>();

  splitInfo ret;
  if (pieces[0].pt() < pieces[1].pt())
    std::reverse(pieces.begin(), pieces.begin() + 1);
  ret.leading_prong = pieces[0];
  ret.subleading_prong = pieces[1];
  ret.zg = j_struct.symmetry();
  ret.rg = j_struct.thetag();
  ret.pt = pieces[0].pt() + pieces[1].pt();

  return ret;
}

int main(int argc, char *argv[]) {

  string usage =
      "Analysis of variation in subjet quantities due to detector response";

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

  // create generator
  dijetcore::PythiaStarSim gen;
  gen.SetPtHatRange(config["pthat_min"].get<double>(),
                    config["pthat_max"].get<double>());
  gen.SetJetPtMin(config["jet_pt_min"].get<double>());
  gen.SetTrackPtMin(config["track_pt_min"].get<double>());
  gen.SetTrackPtMax(config["track_pt_max"].get<double>());
  gen.SetTrackEtaMax(config["track_eta_max"].get<double>());
  gen.SetJetRadius(config["jet_radius"].get<double>());

  double jet_eta_max = config["track_eta_max"].get<double>() -
                       config["jet_radius"].get<double>();
  if (jet_eta_max > 0.0) {
    gen.SetJetEtaMax(jet_eta_max);
  } else {
    LOG(ERROR) << "can not have a jet radius ("
               << config["jet_radius"].get<double>()
               << ") larger than track max eta ("
               << config["track_eta_max"].get<double>() << ")";
    return 1;
  }

  string opt_tag = "opt";
  int opt_idx = 1;
  while (config.count(dijetcore::MakeString(opt_tag, opt_idx)) > 0) {
    gen.AddSettingString(config[dijetcore::MakeString(opt_tag, opt_idx)]);
    opt_idx++;
  }

  string jet_alg = config["jet_alg"];

  if (jet_alg == "kt") {
    gen.SetJetAlgorithm(fastjet::kt_algorithm);
  } else if (jet_alg == "antikt") {
    gen.SetJetAlgorithm(fastjet::antikt_algorithm);
  } else if (jet_alg == "ca") {
    gen.SetJetAlgorithm(fastjet::cambridge_aachen_algorithm);
  } else {
    LOG(ERROR) << "unrecognized jetfinding algorithm: " << jet_alg
               << " select from kt, antikt or ca";
  }

  string sim_mode = config["det_sim_mode"];

  if (sim_mode == "none") {
    gen.SetEfficiencyMode(dijetcore::PythiaStarSim::EfficiencyMode::None);
  } else if (sim_mode == "gauss") {
    gen.SetEfficiencyMode(
        dijetcore::PythiaStarSim::EfficiencyMode::GaussianSmearing);
  } else if (sim_mode == "star") {
    gen.SetEfficiencyMode(dijetcore::PythiaStarSim::EfficiencyMode::STAR);
  } else {
    LOG(ERROR) << "unrecognized simulation mode: " << sim_mode
               << " select from none, gauss or star";
  }

  // initialize the generator
  gen.Initialize();

  // now create output histograms
  TH1D *matched_jet_count =
      new TH1D("matchedcount", ";N_{matched}", 5, -0.5, 4.5);

  // measure of the pT smearing and track loss
  TH2D *track_pt_smearing = new TH2D("ptsmear", ";p_{T};#Delta p_{T}/p_{T}", 50,
                                     0, 10, 50, -0.5, 0.5);
  TH2D *jet_mult = new TH2D("jetmult", ";det N_{part};gen N_{part}", 20, -0.5,
                            19, 20, -0.5, 19);

  TH2D *jet_pt =
      new TH2D("jetpt", ";det p_{T}; gen p_{T}", 50, 0, 50, 50, 0, 50);

  TH2D *jet_zg =
      new TH2D("jetzg", ";det z_{g}; gen z_{g}", 24, 0, 0.6, 24, 0, 0.6);
  TH2D *jet_rg =
      new TH2D("jetrg", ";det R_{g}; gen R_{g}", 25, 0, 1.0, 25, 0, 1.0);
  TH2D *jet_form_1 =
      new TH2D("jetform1", ";det t_{f} split 1; gen t_{f} split 1", 50, 0, 10.0,
               50, 0, 10.0);
  TH2D *jet_form_2 =
      new TH2D("jetform2", ";det t_{f} split 2; gen t_{f} split 2", 50, 0, 10.0,
               50, 0, 10.0);

  TH2D *jet_zg_rel_err_1_2 = new TH2D(
      "jetzgrelerr", ";#Delta z_{g}^{1}/z_{g}^{1}; #Delta z_{g}^{2}/z_{g}^{2}",
      25, -1.0, 1.0, 25, -1.0, 1.0);
  TH2D *jet_rg_rel_err_1_2 = new TH2D(
      "jetrgrelerr", ";#Delta R_{g}^{1}/R_{g}^{1}; #Delta R_{g}^{2}/R_{g}^{2}",
      25, -1.0, 1.0, 25, -1.0, 1.0);
  TH2D *jet_form_rel_err_1_2 =
      new TH2D("jetformrelerr",
               ";#Delta t_{f}^{1}/t_{f}^{1}; #Delta t_{f}^{2}/t_{f}^{2}", 25,
               -1.0, 1.0, 25, -1.0, 1.0);

  // get Iterated SoftDrop variables
  double beta = config["beta"].get<double>();
  double z_cut = config["z_cut"].get<double>();
  double theta_cut = config["theta_cut"].get<double>();
  double r0 = config["angular_scaling"].get<double>();

  // start event loop
  int n_events = config["n_events"].get<int>();
  for (int i = 0; i < n_events; ++i) {

    gen.Run();

    std::vector<fastjet::PseudoJet> gen_jets = gen.GenJets();
    std::vector<fastjet::PseudoJet> det_jets = gen.DetectorJets();

    // create all matches, everything else is a fake
    std::vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet>> jet_pairs;
    std::vector<fastjet::PseudoJet> fakes;
    for (auto &j : det_jets) {
      fastjet::Selector rad_selector =
          fastjet::SelectorCircle(config["match_radius"]);
      rad_selector.set_reference(j);

      std::vector<fastjet::PseudoJet> candidates =
          rad_selector(fastjet::sorted_by_pt(gen_jets));

      if (candidates.size() > 0) {
        fastjet::PseudoJet match = candidates[0];
        jet_pairs.push_back({j, match});
        for (int i = 0; i < gen_jets.size(); ++i) {
          if (gen_jets[i] == match) {
            gen_jets.erase(gen_jets.begin() + i);
          }
        }
      } else {
        fakes.push_back(j);
      }
    }
    // define misses as any generator level jets that had no match
    // (no cuts yet)
    std::vector<fastjet::PseudoJet> misses = gen_jets;

    // save number of matched jets per event
    if (det_jets.size() > 0)
      matched_jet_count->Fill(jet_pairs.size());

    for (auto &m : jet_pairs) {
      fastjet::PseudoJet det = m.first;
      fastjet::PseudoJet gen = m.second;

      // first, look at constituents
      jet_mult->Fill(det.constituents().size(), gen.constituents().size());
      for (auto &p_g : gen.constituents()) {
        if (p_g.user_index() != 0) {
          for (auto &d_g : det.constituents()) {
            if (p_g.delta_R(d_g) < 0.001) {
              track_pt_smearing->Fill(p_g.pt(),
                                      (p_g.pt() - d_g.pt()) / p_g.pt());
            }
          }
        }
      }

      jet_pt->Fill(det.pt(), gen.pt());

      // create RSD
      fastjet::contrib::RecursiveSoftDrop det_rsd(beta, z_cut, -1, r0);
      det_rsd.set_reclustering(true);
      det_rsd.set_verbose_structure(true);
      det_rsd.set_hardest_branch_only(true);
      fastjet::contrib::RecursiveSoftDrop gen_rsd(beta, z_cut, -1, r0);
      gen_rsd.set_reclustering(true);
      gen_rsd.set_verbose_structure(true);
      gen_rsd.set_hardest_branch_only(true);

      fastjet::PseudoJet det_result_jet = det_rsd(det);
      fastjet::PseudoJet gen_result_jet = gen_rsd(gen);
      
      splitInfo det_first_split = parse_split(det_result_jet);
      splitInfo det_second_split = parse_split(det_first_split.leading_prong);
      
      splitInfo gen_first_split = parse_split(gen_result_jet);
      splitInfo gen_second_split = parse_split(gen_first_split.leading_prong);
      
      if (det_first_split.zg > 0.0 && gen_first_split.zg > 0.0) {
        jet_zg->Fill(det_first_split.zg, gen_first_split.zg);
        jet_rg->Fill(det_first_split.rg, gen_first_split.rg);
        double form_det_1 = formation_time(
            det_first_split.zg, det_first_split.pt, det_first_split.rg);
        double form_gen_1 = formation_time(
            gen_first_split.zg, gen_first_split.pt, gen_first_split.rg);
        jet_form_1->Fill(form_det_1, form_gen_1);
      }
      
      if (det_second_split.zg > 0.0 && gen_second_split.zg > 0.0) {
        double form_det_1 = formation_time(
            det_first_split.zg, det_first_split.pt, det_first_split.rg);
        double form_gen_1 = formation_time(
            gen_first_split.zg, gen_first_split.pt, gen_first_split.rg);
        double form_det_2 = formation_time(
            det_second_split.zg, det_second_split.pt, det_second_split.rg);
        double form_gen_2 = formation_time(
            gen_second_split.zg, gen_second_split.pt, gen_second_split.rg);
        double rel_err_form_1 = rel_err(form_gen_1, form_det_1);
        double rel_err_form_2 = rel_err(form_gen_2, form_det_2);

        jet_form_2->Fill(form_det_2, form_gen_2);
        jet_form_rel_err_1_2->Fill(rel_err_form_1, rel_err_form_2);

        double rel_err_zg_1 = rel_err(gen_first_split.zg, det_first_split.zg);
        double rel_err_zg_2 = rel_err(gen_second_split.zg, det_second_split.zg);

        jet_zg_rel_err_1_2->Fill(rel_err_zg_1, rel_err_zg_2);

        double rel_err_rg_1 = rel_err(gen_first_split.rg, det_first_split.rg);
        double rel_err_rg_2 = rel_err(gen_first_split.rg, det_first_split.rg);

        jet_rg_rel_err_1_2->Fill(rel_err_rg_1, rel_err_rg_2);
      }
    }
  }

  out.Write();
  out.Close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}