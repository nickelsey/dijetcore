#include "dijetcore/worker/data/run14_qa_worker/run14_qa_worker.h"

#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/math.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"

#include "TDirectory.h"

namespace dijetcore {

  // defaults set in header
  Run14QAWorker::Run14QAWorker() {}

  bool Run14QAWorker::Init(const std::string& hist_prefix) {

    // binning for run id
    unsigned bin_runid = run_id_map_.size();
    double runid_low = 0.5;
    double runid_high = runid_low + bin_runid;

    // binning for refmult/grefmult
    unsigned bin_refmult = 100;
    double refmult_low = 0;
    double refmult_high = 800;

    // binning for nprim/nglobal
    unsigned bin_nprim = 100;
    unsigned bin_nglob = 100;
    double nprim_low = 0;
    double nprim_high = 1400;
    double nglob_low = 0;
    double nglob_high = 4500;

    // binning for zdc
    unsigned bin_zdc = 50;
    double zdc_low = 0;
    double zdc_high = 100;

    // binnning for bbc
    unsigned bin_bbc = 50;
    double bbc_low = 0;
    double bbc_high = 100;

    // binning for momentum
    unsigned bin_momentum = 200;
    double momentum_low = 0;
    double momentum_high = 100;

    // binning for energy
    unsigned bin_energy = 200;
    double energy_low = 0;
    double energy_high = 100;

    // binning for dca
    unsigned bin_dca = 50;
    double dca_low = 0;
    double dca_high = 3.0;

    // binning for nhit/nhitposs
    unsigned bin_nhit = 50;
    double nhit_low = 0;
    double nhit_high = 50;

    // binning for nhitfrac
    unsigned bin_nhitfrac = 50;
    double nhitfrac_low = 0;
    double nhitfrac_high = 1.0;

    // binning for ADC
    unsigned bin_adc = 50;
    double adc_low = 0;
    double adc_high = 500;

    // binning for eta
    unsigned bin_eta = 40;
    double eta_low = -1.0;
    double eta_high = 1.0;

    // binning for phi
    unsigned bin_phi = 40;
    double phi_low = -math::pi;
    double phi_high = math::pi;

    // vertex z binning
    unsigned bin_vz = 60;
    double vz_low = -30;
    double vz_high = 30;

    // vertex xy binning
    unsigned bin_vxy = 50;
    double vxy_low = -5.0;
    double vxy_high = 5.0;

    // vertex dvz binning
    unsigned bin_dvz = 30;
    double dvz_low = -3.0;
    double dvz_high = 3.0;

    // tower binning
    unsigned bin_tower = 4800;
    double tower_low = 0.5;
    double tower_high = 4800.5;

    // vertex binning
    unsigned bin_vertices = 30;
    unsigned vertices_low = 0;
    unsigned vertices_high = 30;

    // Histograms will calculate gaussian errors
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();

    // event QA
    zdc_vz_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcvz").c_str(), ";zdc [kHz];V_{z}", bin_zdc, zdc_low, zdc_high, bin_vz, vz_low, vz_high);
    vz_vx_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "vzvx").c_str(), ";V_{z};V_{x}", bin_vz, vz_low, vz_high, bin_vxy, vxy_low, vxy_high);
    vz_vy_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "vzvy").c_str(), ";V_{z};V_{y}", bin_vz, vz_low, vz_high, bin_vxy, vxy_low, vxy_high);
    zdc_refmult_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcref").c_str(), ";zdc [kHz];refmult", bin_zdc, zdc_low, zdc_high, bin_refmult, refmult_low, refmult_high);
    bbc_refmult_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "bbcref").c_str(), ";bbc [kHz];refmult", bin_bbc, bbc_low, bbc_high, bin_refmult, refmult_low, refmult_high);
    n_vertices_ = std::make_unique<TH1F>(dijetcore::MakeString(hist_prefix, "nvertex").c_str(), ";N_{vertices};fraction", bin_vertices, vertices_low, vertices_high);

    // initialize histograms
    if (bin_runid > 0) {
      run_id_refmult_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidref").c_str(), ";run ID;refmult", bin_runid, runid_low, runid_high, bin_refmult, refmult_low, refmult_high);
      run_id_grefmult_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidgref").c_str(), ";run ID;grefmult", bin_runid, runid_low, runid_high, bin_refmult, refmult_low, refmult_high);
      run_id_nprim_nglob_ = std::make_unique<TH3F>(dijetcore::MakeString(hist_prefix, "runidnprimnglob").c_str(), ";run ID;N_{prim};N_{glob}", bin_runid, runid_low, runid_high, bin_nprim, nprim_low, nprim_high, bin_nglob, nglob_low, nglob_high);
      run_id_zdc_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidzdc").c_str(), ";run ID; zdc [kHz]", bin_runid, runid_low, runid_high, bin_zdc, zdc_low, zdc_high);
      run_id_bbc_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidbbc").c_str(), ";run ID; bbc [kHz]", bin_runid, runid_low, runid_high, bin_bbc, bbc_low, bbc_high);
      run_id_vzvpdvz_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runiddvz").c_str(), ";run ID; dV_{z} [cm]", bin_runid, runid_low, runid_high, bin_dvz, dvz_low, dvz_high);
      run_id_vz_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidvz").c_str(), "run ID; V_{z}", bin_runid, runid_low, runid_high, bin_vz, vz_low, vz_high);
      run_id_vx_vy_ = std::make_unique<TH3F>(dijetcore::MakeString(hist_prefix, "runidvxvy").c_str(), ";run ID;V_{x};V_{y}", bin_runid, runid_low, runid_high, bin_vxy, vxy_low, vxy_high, bin_vxy, vxy_low, vxy_high);
    }

    // track QA
    if (do_track_qa_) {
      px_py_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "pxpy").c_str(), ";p_{x};p_{y}", bin_momentum, momentum_low, momentum_high, bin_momentum, momentum_low, momentum_high);
      pz_px_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "pzpx").c_str(), ";p_{z};p_{x}", bin_momentum, momentum_low, momentum_high, bin_momentum, momentum_low, momentum_high);
      pz_py_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "pzpy").c_str(), ";p_{z};p_{y}", bin_momentum, momentum_low, momentum_high, bin_momentum, momentum_low, momentum_high);
      zdc_px_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcpx").c_str(), ";zdc [kHz];p_{x}", bin_zdc, zdc_low, zdc_high, bin_momentum, momentum_low, momentum_high);
      zdc_py_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcpy").c_str(), ";zdc [kHz];p_{y}", bin_zdc, zdc_low, zdc_high, bin_momentum, momentum_low, momentum_high);
      zdc_pz_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcpz").c_str(), ";zdc [kHz];p_{z}", bin_zdc, zdc_low, zdc_high, bin_momentum, momentum_low, momentum_high);
      zdc_pt_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcpt").c_str(), ";zdc [kHz];p_{T}", bin_zdc, zdc_low, zdc_high, bin_momentum, momentum_low, momentum_high);
      zdc_dca_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcdca").c_str(), ";zdc [kHz];DCA [cm]", bin_zdc, zdc_low, zdc_high, bin_dca, dca_low, dca_high);
      zdc_nhit_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcnhit").c_str(), ";zdc [kHz];N_{hit}", bin_zdc, zdc_low, zdc_high, bin_nhit, nhit_low, nhit_high);
      zdc_nhitposs_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcnhitposs").c_str(), ";zdc [kHz];N_{hit} poss", bin_zdc, zdc_low, zdc_high, bin_nhit, nhit_low, nhit_high);
      zdc_nhitfrac_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcnhitfrac").c_str(), ";zdc [kHz];N_{hit} frac", bin_zdc, zdc_low, zdc_high, bin_nhitfrac, nhitfrac_low, nhitfrac_high);
      zdc_eta_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdceta").c_str(), ";zdc [kHz];eta", bin_zdc, zdc_low, zdc_high, bin_eta, eta_low, eta_high);
      zdc_phi_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcphi").c_str(), ";zdc [kHz];phi", bin_zdc, zdc_low, zdc_high, bin_phi, phi_low, phi_high);
      eta_phi_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "tracketaphi").c_str(), ";eta;phi", bin_eta, eta_low, eta_high, bin_phi, phi_low, phi_high);
      // track runID variables
      if (run_id_map_.size() > 0) {
        run_id_track_pt_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidpt").c_str(), ";run ID;p_{T}", bin_runid, runid_low, runid_high, bin_momentum, momentum_low, momentum_high);
        run_id_track_px_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidpx").c_str(), ";run ID;p_{x}", bin_runid, runid_low, runid_high, bin_momentum, momentum_low, momentum_high);
        run_id_track_py_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidpy").c_str(), ";run ID;p_{y}", bin_runid, runid_low, runid_high, bin_momentum, momentum_low, momentum_high);
        run_id_track_pz_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidpz").c_str(), ";run ID;p_{z}", bin_runid, runid_low, runid_high, bin_momentum, momentum_low, momentum_high);
        run_id_track_eta_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runideta").c_str(), ";run ID;eta", bin_runid, runid_low, runid_high, bin_eta, eta_low, eta_high);
        run_id_track_phi_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidphi").c_str(), ";run ID;phi", bin_runid, runid_low, runid_high, bin_phi, phi_low, phi_high);
        run_id_track_dca_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runiddca").c_str(), ";run ID;dca [cm]", bin_runid, runid_low, runid_high, bin_dca, dca_low, dca_high);
        run_id_track_nhit_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidnhit").c_str(), ";run ID;N_{hit}", bin_runid, runid_low, runid_high, bin_nhit, nhit_low, nhit_high);
        run_id_track_nhitposs_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidnhitposs").c_str(), ";run ID;N_{hit} poss", bin_runid, runid_low, runid_high, bin_nhit, nhit_low, nhit_high);
        run_id_track_nhitfrac_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "runidnhitfrac").c_str(), ";run ID;N_{hit} frac", bin_runid, runid_low, runid_high, bin_nhitfrac, nhitfrac_low, nhitfrac_high);
      }
    }

    // tower QA
    if (do_tower_qa_) {
      e_et_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "eet").c_str(), ";E;E_{T}", bin_energy, energy_low, energy_high, bin_energy, energy_low, energy_high);
      zdc_e_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdce").c_str(), ";zdc [kHz];E", bin_zdc, zdc_low, zdc_high, bin_energy, energy_low, energy_high);
      zdc_et_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcet").c_str(), ";zdc [kHz];E_{T}", bin_zdc, zdc_low, zdc_high, bin_energy, energy_low, energy_high);
      zdc_adc_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "zdcadc").c_str(), ";zdc [kHz];ADC", bin_zdc, zdc_low, zdc_high, bin_adc, adc_low, adc_high);
      tow_eta_phi_ = std::make_unique<TH2F>(dijetcore::MakeString(hist_prefix, "towetaphi").c_str(), ";eta;phi", bin_eta, eta_low, eta_high, bin_phi, phi_low, phi_high);
      // tower runID variables
      if (run_id_map_.size() > 0) {
        // create THnSparse
        int tow_dims = 3;
        int tow_bins[] = {bin_runid, bin_tower, bin_energy};
        double tow_dims_low[] = {runid_low, tower_low, energy_low};
        double tow_dims_high[] = {runid_high, tower_high, energy_high};
        int tow_bins_adc[] = {bin_runid, bin_tower, bin_adc};
        double tow_dims_low_adc[] = {runid_low, tower_low, adc_low};
        double tow_dims_high_adc[] = {runid_high, tower_high, adc_high};
        run_id_tower_e_ = std::make_unique<THnSparseS>(dijetcore::MakeString(hist_prefix, "runide").c_str(), ";run ID;tower ID;E", tow_dims, tow_bins, tow_dims_low, tow_dims_high);
        run_id_tower_et_ = std::make_unique<THnSparseS>(dijetcore::MakeString(hist_prefix, "runidet").c_str(), ";run ID;tower ID;E_{T}", tow_dims, tow_bins, tow_dims_low, tow_dims_high);
        run_id_tower_adc_ = std::make_unique<THnSparseS>(dijetcore::MakeString(hist_prefix, "runidadc").c_str(), ";run ID;tower ID;ADC", tow_dims, tow_bins_adc, tow_dims_low_adc, tow_dims_high_adc);

        run_id_tower_e_->Sumw2();
        run_id_tower_et_->Sumw2();
        run_id_tower_adc_->Sumw2();
      }
    }

    initialized_ = true;
    return true;
  }

  bool Run14QAWorker::WriteTo(TFile& file) {
    if (!file.IsOpen()) {
      LOG(ERROR) << "input file " << file.GetName() << " is not open, failed to write histograms to disk";
      return false;
    }
    
    // record current directory
    TDirectory* current_dir = TDirectory::CurrentDirectory();
    // change to TFile to write out
    file.cd();
    
    // write event histograms after checking for initialization
    if (zdc_vz_ == nullptr) {
      LOG(ERROR) << "histograms do not exist: was Init() called?";
      return false;
    }

    zdc_vz_->Write();
    vz_vx_->Write();
    vz_vy_->Write();
    zdc_refmult_->Write();
    bbc_refmult_->Write();
    n_vertices_->Write();

    // check if we have runid observables
    if (run_id_refmult_ != nullptr) {
      run_id_refmult_->Write();
      run_id_grefmult_->Write();
      run_id_nprim_nglob_->Write();
      run_id_zdc_->Write();
      run_id_bbc_->Write();
      run_id_vzvpdvz_->Write();
      run_id_vz_->Write();
      run_id_vx_vy_->Write();
    }

    // check if we have track observables 
    if (px_py_ != nullptr) {
      px_py_->Write();
      pz_px_->Write();
      pz_py_->Write();
      zdc_px_->Write();
      zdc_py_->Write();
      zdc_pz_->Write();
      zdc_pt_->Write();
      zdc_dca_->Write();
      zdc_nhit_->Write();
      zdc_nhitposs_->Write();
      zdc_nhitfrac_->Write();
      zdc_eta_->Write();
      zdc_phi_->Write();
      eta_phi_->Write();
      if (run_id_track_pt_ != nullptr) {
        run_id_track_pt_->Write();
        run_id_track_px_->Write();
        run_id_track_py_->Write();
        run_id_track_pz_->Write();
        run_id_track_eta_->Write();
        run_id_track_phi_->Write();
        run_id_track_dca_->Write();
        run_id_track_nhit_->Write();
        run_id_track_nhitposs_->Write();
        run_id_track_nhitfrac_->Write();
      }
    }

    // check if we have tower observables
    LOG(ERROR) << "writing tower?";
    if (e_et_ != nullptr) {
      LOG(ERROR) << "writing tower";
      e_et_->Write();
      zdc_e_->Write();
      zdc_et_->Write();
      zdc_adc_->Write();
      tow_eta_phi_->Write();
      if (run_id_tower_e_ != nullptr) {
        LOG(ERROR) << "writing tower runid?";
        run_id_tower_e_->Write();
        run_id_tower_et_->Write();
        run_id_tower_adc_->Write();
      }
    }

    // change back to current_dir
    current_dir->cd();

    return true;
  }

  void Run14QAWorker::DoRunQA(std::set<unsigned>& run_ids) {
    unsigned bin = 1;
    run_id_map_.clear();
    for (auto& runid : run_ids) {
      run_id_map_.insert({runid, bin});
      bin++;
    }
    if (run_id_map_.size() > 0)
      do_run_qa_ = true;
  }

  bool Run14QAWorker::Run(TStarJetPicoReader& reader) {

    // first do event level QA
    TStarJetPicoEventHeader* header = reader.GetEvent()->GetHeader();
    
    // make sure initialization has happened
    if (initialized_ == false)
      Init();
    
    zdc_vz_->Fill(header->GetZdcCoincidenceRate()/1000.0, header->GetPrimaryVertexZ());
    vz_vx_->Fill(header->GetPrimaryVertexZ(), header->GetPrimaryVertexX());
    vz_vy_->Fill(header->GetPrimaryVertexZ(), header->GetPrimaryVertexY());
    zdc_refmult_->Fill(header->GetZdcCoincidenceRate()/1000.0, header->GetReferenceMultiplicity());
    bbc_refmult_->Fill(header->GetBbcCoincidenceRate()/1000.0, header->GetReferenceMultiplicity());
    n_vertices_->Fill(header->GetNumberOfVertices());

    // check if we are doing run-by-run QA
    if (run_id_map_.size() > 0 && run_id_refmult_ != nullptr) {
      if (RunQA(reader) == false) {
        LOG(ERROR) << "failure in RunQA";
      }
    }
    
    // check if we are doing tower QA
    if (e_et_ != nullptr) {
      if (!TowerQA(reader)) {
        LOG(ERROR) << "failure in TowerQA";
      }
    }
    
    if (px_py_ != nullptr) {
      if (!TrackQA(reader)) {
        LOG(ERROR) << "failure in TrackQA";
      }
    }

    return true;
  }

  bool Run14QAWorker::RunQA(TStarJetPicoReader& event) {
    TStarJetPicoEventHeader* header = event.GetEvent()->GetHeader();

    unsigned runid_idx = run_id_map_[header->GetRunId()];

    run_id_refmult_->Fill(runid_idx, header->GetReferenceMultiplicity());
    run_id_grefmult_->Fill(runid_idx, header->GetGReferenceMultiplicity());
    run_id_nprim_nglob_->Fill(runid_idx, header->GetReferenceMultiplicity(), header->GetGReferenceMultiplicity());
    run_id_zdc_->Fill(runid_idx, header->GetZdcCoincidenceRate()/1000.0);
    run_id_bbc_->Fill(runid_idx, header->GetBbcCoincidenceRate()/1000.0);
    run_id_vzvpdvz_->Fill(runid_idx, header->GetPrimaryVertexZ() - header->GetVpdVz());
    run_id_vz_->Fill(runid_idx, header->GetPrimaryVertexZ());
    run_id_vx_vy_->Fill(runid_idx, header->GetPrimaryVertexX(), header->GetPrimaryVertexY());

    return true;
  }

  bool Run14QAWorker::TowerQA(TStarJetPicoReader& reader) {
    TStarJetPicoEventHeader* header = reader.GetEvent()->GetHeader();
    TList* towers = reader.GetListOfSelectedTowers();
    TIter nextTower(towers);

    while(TStarJetPicoTower* tower = (TStarJetPicoTower*) nextTower()) {
      e_et_->Fill(tower->GetEnergy(), tower->GetEnergy()/cosh(tower->GetEtaCorrected()));
      zdc_e_->Fill(header->GetZdcCoincidenceRate()/1000, tower->GetEnergy());
      zdc_et_->Fill(header->GetZdcCoincidenceRate()/1000, tower->GetEnergy()/cosh(tower->GetEtaCorrected()));
      zdc_adc_->Fill(header->GetZdcCoincidenceRate()/1000, tower->GetADC());
      tow_eta_phi_->Fill(tower->GetEta(), tower->GetPhi());
      if (run_id_tower_e_ != nullptr) {
        unsigned runid_idx = run_id_map_[header->GetRunId()];
        double e_[] = {runid_idx, tower->GetId(), tower->GetEnergy()};
        double et_[] = {runid_idx, tower->GetId(), tower->GetEnergy()/cosh(tower->GetEtaCorrected())};
        double adc_[] = {runid_idx, tower->GetId(),tower->GetADC()};
        run_id_tower_e_->Fill(e_);
        run_id_tower_et_->Fill(et_);
        run_id_tower_adc_->Fill(adc_);
      }
    }
    
    return true;
  }

  bool Run14QAWorker::TrackQA(TStarJetPicoReader& reader) {
    TStarJetPicoEventHeader* header = reader.GetEvent()->GetHeader();
    TList* tracks = reader.GetListOfSelectedTracks();
    TIter nextTrack(tracks);

    while(TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) nextTrack()) {
      px_py_->Fill(track->GetPx(), track->GetPy());
      pz_px_->Fill(track->GetPz(), track->GetPx());
      pz_py_->Fill(track->GetPz(), track->GetPy());
      zdc_px_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetPx());
      zdc_py_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetPy());
      zdc_pz_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetPz());
      zdc_pt_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetPt());
      zdc_dca_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetDCA());
      zdc_nhit_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetNOfFittedHits());
      zdc_nhitposs_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetNOfPossHits());
      zdc_nhitfrac_->Fill(header->GetZdcCoincidenceRate()/1000, ((double)track->GetNOfFittedHits())/track->GetNOfPossHits());
      zdc_eta_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetEta());
      zdc_phi_->Fill(header->GetZdcCoincidenceRate()/1000, track->GetPhi());
      eta_phi_->Fill(track->GetEta(), track->GetPhi());
      if (run_id_track_pt_ != nullptr) {
        unsigned runid_idx = run_id_map_[header->GetRunId()];
        run_id_track_pt_->Fill(runid_idx, track->GetPt());
        run_id_track_px_->Fill(runid_idx, track->GetPx());
        run_id_track_py_->Fill(runid_idx, track->GetPy());
        run_id_track_pz_->Fill(runid_idx, track->GetPz());
        run_id_track_eta_->Fill(runid_idx, track->GetEta());
        run_id_track_phi_->Fill(runid_idx, track->GetPhi());
        run_id_track_dca_->Fill(runid_idx, track->GetDCA());
        run_id_track_nhit_->Fill(runid_idx, track->GetNOfFittedHits());
        run_id_track_nhitposs_->Fill(runid_idx, track->GetNOfPossHits());
        run_id_track_nhitfrac_->Fill(runid_idx, ((double) track->GetNOfFittedHits())/track->GetNOfPossHits());
      }
    }
    return true;
  }

} // namespace dijetcore