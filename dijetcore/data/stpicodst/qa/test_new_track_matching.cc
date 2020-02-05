#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/data/centrality/centrality_run7.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/trigger_lookup.h"

#include "jetreader/reader/bemc_helper.h"
#include "jetreader/reader/event_selector.h"
#include "jetreader/reader/reader.h"
#include "jetreader/reader/tower_selector.h"
#include "jetreader/reader/track_selector.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TTree.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "output_directory", // directory for output root file to written to
    "bad_tower_list",   // bad tower list for TStarJetPicoReader
    "bad_run_list",     // bad run list for TStarJetPicoReader
    "triggers",         // trigger id set (defined in
                        // dijetcore/util/data/trigger_lookup.h)
};

double dphi(double p1, double p2) {
  double dp = p1 - p2;
  while (dp > TMath::Pi())
    dp -= 2.0 * TMath::Pi();
  while (dp < -TMath::Pi())
    dp += 2.0 * TMath::Pi();
  return dp;
}

int main(int argc, char *argv[]) {
  string usage = "y14 qa routine";

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

  // create output file from the given directory, name & id
  string outfile_name =
      output_dir + "/" + FLAGS_name + dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // initialize the reader
  jetreader::Reader reader(FLAGS_input);

  if (!std::string(config["bad_run_list"]).empty())
    reader.eventSelector()->addBadRuns(std::string(config["bad_run_list"]));
  if (!std::string(config["triggers"]).empty())
    reader.eventSelector()->addTriggerIds(config["triggers"]);
  reader.eventSelector()->setVzRange(-30, 30);
  reader.eventSelector()->setdVzMax(3.0);

  reader.trackSelector()->setDcaMax(3.0);
  reader.trackSelector()->setNHitsMin(15);
  reader.trackSelector()->setNHitsFracMin(0.52);
  reader.trackSelector()->setPtMax(30.0);
  reader.trackSelector()->setPtMin(0.2);
  reader.trackSelector()->rejectEventOnPtFailure();

  if (!std::string(config["bad_tower_list"]).empty())
    reader.towerSelector()->addBadTowers(std::string(config["bad_tower_list"]));
  reader.towerSelector()->setEtMax(30.0);
  reader.towerSelector()->setEtMin(0.2);
  reader.towerSelector()->rejectEventOnEtFailure();

  reader.init();

  // turn on centrality for low/mid luminosity
  reader.centrality().loadCentralityDef(jetreader::CentDefId::Run14LowMid);

  // create BEMC helper to match tracks and towers
  jetreader::BemcHelper bemc_helper;

  // histograms
  TH1D *hvz = new TH1D("vz", ";vz{z}", 70, -35, 35);
  TH1D *hrefmult = new TH1D("refmult", ";refmult", 800, 0, 800);
  TH1D *hgrefmult = new TH1D("grefmult", ";grefmult", 800, 0, 800);
  TH1D *hlumi = new TH1D("lumi", ";zdcX [kHz]", 100, 0, 100);
  TH1D *hcent = new TH1D("cent", ";Centrality (16)", 16, -0.5, 15.5);

  TH1D *track_phi = new TH1D("trackphi", ";#phi", 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_eta = new TH1D("tracketa", ";#eta", 70, -1.5, 1.5);
  TH2D *track_pt_dphi =
      new TH2D("trackptdphi", ";p_{T};#Delta#phi", 10, 0, 10, 30, -1.0, 1.0);
  TH2D *track_pt_deta =
      new TH2D("trackptdeta", ";p_{T};#Delta#eta", 10, 0, 10, 30, -1.0, 1.0);

  TH1D *track_phi_global =
      new TH1D("trackphiglobal", ";#phi", 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_eta_global = new TH1D("tracketaglobal", ";#eta", 70, -1.5, 1.5);
  TH2D *track_pt_dphi_global = new TH2D(
      "trackptdphiglobal", ";p_{T};#Delta#phi", 10, 0, 10, 30, -1.0, 1.0);
  TH2D *track_pt_deta_global = new TH2D(
      "trackptdetaglobal", ";p_{T};#Delta#eta", 10, 0, 10, 30, -1.0, 1.0);

  TH1D *primary_track_match_success =
      new TH1D("primMatch", "no match/match", 2, -2, 2);
  TH1D *global_track_match_success =
      new TH1D("globMatch", "no match/match", 2, -2, 2);
  TH1D *tow_idx = new TH1D("towidx", ";Tower Index", 4800, -0.5, 4799.5);

  TH3D* track_miss = new TH3D("trackmiss", ";p_{T};#eta;#phi", 100, 0, 10, 70, -1.5, 1.5, 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_phi_miss = new TH1D("trackphimiss", ";#phi", 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_phi_miss_inner = new TH1D("trackphimissinner", ";#phi", 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_phi_miss_inner_pt = new TH1D("trackphimissinnerpt", ";#phi", 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_eta_miss = new TH1D("tracketamiss", ";#eta", 70, -1.5, 1.5);

  TH3D* track_miss_global = new TH3D("trackmissglobal", ";p_{T};#eta;#phi", 100, 0, 10, 70, -1.5, 1.5, 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_phi_miss_global = new TH1D("trackphimissglobal", ";#phi", 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_phi_miss_inner_global = new TH1D("trackphimissinnerglobal", ";#phi", 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_phi_miss_inner_global_pt = new TH1D("trackphimissinnerglobalpt", ";#phi", 70, -TMath::Pi(), TMath::Pi());
  TH1D *track_eta_miss_global = new TH1D("tracketamissglobal", ";#eta", 70, -1.5, 1.5);

  TH1D* track_miss_pos_vz_eta = new TH1D("trackmissposvzeta", ";#eta;counts", 70, -1.5, 1.5);
  TH1D* track_miss_neg_vz_eta = new TH1D("trackmissnegvzeta", ";#eta;counts", 70, -1.5, 1.5);

  try {
    while (reader.next()) {

      // set event quantities for tree
      unsigned runid = reader.picoDst()->event()->runId();
      unsigned eventid = reader.picoDst()->event()->eventId();
      unsigned ref = reader.picoDst()->event()->refMult();
      unsigned gref = reader.picoDst()->event()->grefMult();
      double vz = reader.picoDst()->event()->primaryVertex().Z();
      unsigned cent = reader.centrality16();

      if (cent < 0 || cent >= 16)
        continue;

      hvz->Fill(vz);
      hrefmult->Fill(ref);
      hgrefmult->Fill(gref);
      hlumi->Fill(reader.picoDst()->event()->ZDCx() / 1000.0);
      hcent->Fill(cent);

      for (int i = 0; i < reader.picoDst()->numberOfTracks(); ++i) {
        StPicoTrack *track = reader.picoDst()->track(i);

        if (track->gPt() < 0.2)
          continue;
        if (track->gDCA(reader.picoDst()->event()->primaryVertex()).Mag() > 0.2)
          continue;
        if (track->nHits() < 15)
          continue;
        if ((double)track->nHits() / track->nHitsPoss() < 0.52)
          continue;
        if (track->pPt() > 0.2 && fabs(track->pMom().Eta()) > 1.0)
          continue;

        double gphi = track->gMom().Phi();
        double geta = track->gMom().Eta();

        track_eta_global->Fill(geta);
        track_phi_global->Fill(gphi);

        int bemc_index = track->bemcTowerIndex();
        StPicoBTowHit *tow_hit = nullptr;

        if (bemc_index > -1) {
          tow_idx->Fill(bemc_index);
          tow_hit = reader.picoDst()->btowHit(bemc_index);
          global_track_match_success->Fill(1);
          double tow_eta = bemc_helper.vertexCorrectedEta(bemc_index + 1, vz);
          double tow_phi = bemc_helper.towerPhi(bemc_index + 1);
          track_pt_deta_global->Fill(track->gPt(), tow_eta - geta);
          track_pt_dphi_global->Fill(track->gPt(), dphi(tow_phi, gphi));
        } else {
          track_miss_global->Fill(track->gPt(), geta, gphi);
          track_phi_miss_global->Fill(gphi);
          track_eta_miss_global->Fill(geta);
          if (abs(geta) < 0.5) {
            track_phi_miss_inner_global->Fill(gphi);
            if (track->gPt() > 2.0)
             track_phi_miss_inner_global_pt->Fill(gphi);
          }
          global_track_match_success->Fill(-1);
        }
        if (track->pPt() > 0.2) {
          double peta = track->pMom().Eta();
          double pphi = track->pMom().Phi();
          if (bemc_index > -1) {
            tow_hit = reader.picoDst()->btowHit(bemc_index);
            primary_track_match_success->Fill(1);
            double tow_eta = bemc_helper.vertexCorrectedEta(bemc_index + 1, vz);
            double tow_phi = bemc_helper.towerPhi(bemc_index + 1);
            track_eta->Fill(peta);
            track_phi->Fill(pphi);
            track_pt_deta->Fill(track->pPt(), tow_eta - peta);
            track_pt_dphi->Fill(track->pPt(), dphi(tow_phi, pphi));
          }
          else {
            track_miss->Fill(track->pPt(), peta, pphi);
            track_phi_miss->Fill(pphi);
            track_eta_miss->Fill(peta);
            if (abs(peta) < 0.5) {
              track_phi_miss_inner->Fill(pphi);
              if (track->pPt() > 2.0)
                track_phi_miss_inner_pt->Fill(gphi);
            }
            primary_track_match_success->Fill(-1);
            if (vz > 10.0) {
              track_miss_pos_vz_eta->Fill(peta);
            }
            else if (vz < 10.0) {
              track_miss_neg_vz_eta->Fill(peta);
            }
          }
        }
      }
    }
  } catch (std::exception &e) {
    LOG(ERROR) << "Caught: " << e.what() << " during analysis loop.";
  }

  out.Write();
  out.Close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}