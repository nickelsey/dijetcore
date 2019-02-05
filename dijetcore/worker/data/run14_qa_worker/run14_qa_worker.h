#ifndef DIJETCORE_WORKER_DATA_RUN14_QA_WORKER_RUN14_QA_WORKER_H
#define DIJETCORE_WORKER_DATA_RUN14_QA_WORKER_RUN14_QA_WORKER_H

#include "dijetcore/lib/memory.h"

#include <set>
#include <map>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"

#include "TStarJetPicoReader.h"

namespace dijetcore {

  class Run14QAWorker {
  public:
    // constructor does not initialize histograms - must call Init() before running
    Run14QAWorker();

    // Initializes all histograms. Must be called before Run() is called. Any options
    // such as turning off tower/track QA or turning on run QA must be called before
    //  Init(). hist_prefix is an optional name that can be tagged on to the front of 
    // each histogram's root identifier, to keep different jobs separate after adding
    // root files
    bool Init(const std::string& hist_prefix = "");

    // Writes all initialized histograms to specified TFile
    bool WriteTo(TFile& file);

    // by default, QA for tracks and towers are on. However, turning them off if they
    // are not needed can reduce runtime and memory footprint (which is large).
    void DoTowerQA(bool flag) {do_tower_qa_ = flag;}
    void DoTrackQA(bool flag) {do_track_qa_ = flag;}

    // To turn on run-by-run QA, the maker needs a list of run IDs. Supply that here
    void DoRunQA(std::set<unsigned>& run_ids);

    // analyze event
    bool Run(TStarJetPicoReader& reader);

  private:

    bool initialized_ = false;

    bool RunQA(TStarJetPicoReader& reader);
    bool TowerQA(TStarJetPicoReader& event);
    bool TrackQA(TStarJetPicoReader& event);

    bool do_run_qa_ = false;
    bool do_tower_qa_ = true;
    bool do_track_qa_ = true;

    std::map<unsigned, unsigned> run_id_map_ = {};

    // histograms

    // run-by-run
    unique_ptr<TH3F> run_id_ref_gref_;
    unique_ptr<TH3F> run_id_nprim_nglob_;
    unique_ptr<TH2F> run_id_zdc_;
    unique_ptr<TH2F> run_id_bbc_;
    unique_ptr<TH2F> run_id_vzvpdvz_;
    unique_ptr<TH2F> run_id_vz_;
    unique_ptr<TH3F> run_id_vx_vy_;
    unique_ptr<TH2F> run_id_track_pt_;
    unique_ptr<TH2F> run_id_track_px_;
    unique_ptr<TH2F> run_id_track_py_;
    unique_ptr<TH2F> run_id_track_pz_;
    unique_ptr<TH2F> run_id_track_eta_;
    unique_ptr<TH2F> run_id_track_phi_;
    unique_ptr<TH2F> run_id_track_dca_;
    unique_ptr<TH2F> run_id_track_nhit_;
    unique_ptr<TH2F> run_id_track_nhitposs_;
    unique_ptr<TH2F> run_id_track_nhitfrac_;
    unique_ptr<THnSparseS> run_id_tower_e_;
    unique_ptr<THnSparseS> run_id_tower_et_;
    unique_ptr<THnSparseS> run_id_tower_adc_;

    // event QA
    unique_ptr<TH2F> zdc_vz_;
    unique_ptr<TH2F> vz_vx_;
    unique_ptr<TH2F> vz_vy_;
    unique_ptr<TH2F> zdc_refmult_;
    unique_ptr<TH2F> bbc_refmult_;
    unique_ptr<TH1F> n_vertices_;

    // track QA
    unique_ptr<TH2F> px_py_;
    unique_ptr<TH2F> pz_px_;
    unique_ptr<TH2F> pz_py_;
    unique_ptr<TH2F> zdc_px_;
    unique_ptr<TH2F> zdc_py_;
    unique_ptr<TH2F> zdc_pz_;
    unique_ptr<TH2F> zdc_pt_;
    unique_ptr<TH2F> zdc_dca_;
    unique_ptr<TH2F> zdc_nhit_;
    unique_ptr<TH2F> zdc_nhitposs_;
    unique_ptr<TH2F> zdc_nhitfrac_;
    unique_ptr<TH2F> zdc_eta_;
    unique_ptr<TH2F> zdc_phi_;
    unique_ptr<TH2F> eta_phi_;

    // tower QA
    unique_ptr<TH2F> e_et_;
    unique_ptr<TH2F> zdc_e_;
    unique_ptr<TH2F> zdc_et_;
    unique_ptr<TH2F> zdc_adc_;
    unique_ptr<TH2F> tow_eta_phi_;

  };

} // namespace dijetcore

#endif // DIJETCORE_WORKER_DATA_RUN14_QA_WORKER_RUN14_QA_WORKER_H