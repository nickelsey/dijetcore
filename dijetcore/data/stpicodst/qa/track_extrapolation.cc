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

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "", "input file");
DIJETCORE_DEFINE_string(name, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(output_dir, "", "output directory name");

void DrawExtrapolation(std::string name, std::vector<TVector3> &i,
                       std::vector<TVector3> &o);

int main(int argc, char *argv[]) {
  string usage = "track extrapolation qa";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // check to make sure the input file paths are sane
  if (!boost::filesystem::exists(FLAGS_input)) {
    LOG(ERROR) << "input file does not exist: " << FLAGS_input;
    return 1;
  }

  TFile *input_file = new TFile(FLAGS_input.c_str(), "READ");

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_output_dir);
  boost::filesystem::create_directories(dir);

  // create output file from the given directory, name & id
  string outfile_name = FLAGS_output_dir + "/" + FLAGS_name +
                        dijetcore::MakeString(FLAGS_id) + ".root";
  TFile out(outfile_name.c_str(), "RECREATE");

  // create TTreeReaders for event and track trees
  TTreeReader events("events", input_file);
  TTreeReaderValue<unsigned> runid(events, "runid");
  TTreeReaderValue<unsigned> eventid(events, "eventid");
  TTreeReaderValue<TVector3> pvertex(events, "vertex");
  TTreeReaderValue<double> zdcx(events, "zdcx");
  TTreeReaderValue<double> b(events, "b");
  TTreeReaderValue<double> refmult(events, "refmult");
  TTreeReaderValue<double> grefmult(events, "grefmult");
  TTreeReaderValue<unsigned> ntrack(events, "ntrack");

  TTreeReader tracks("tracks", input_file);
  TTreeReaderValue<std::vector<TVector3>> tpi(tracks, "pointsi");
  TTreeReaderValue<std::vector<TVector3>> tpo(tracks, "pointso");
  TTreeReaderValue<std::vector<TVector3>> gpi(tracks, "gpointsi");
  TTreeReaderValue<std::vector<TVector3>> gpo(tracks, "gpointso");
  TTreeReaderValue<TLorentzVector> primary(tracks, "vec");
  TTreeReaderValue<TLorentzVector> global(tracks, "global");
  TTreeReaderValue<TVector3> tihelix(tracks, "ihelix");
  TTreeReaderValue<TVector3> tohelix(tracks, "ohelix");
  TTreeReaderValue<TVector3> gihelix(tracks, "gihelix");
  TTreeReaderValue<TVector3> gohelix(tracks, "gohelix");
  TTreeReaderValue<TVector3> timom(tracks, "imom");
  TTreeReaderValue<TVector3> tomom(tracks, "omom");
  TTreeReaderValue<TVector3> gimom(tracks, "gimom");
  TTreeReaderValue<TVector3> gomom(tracks, "gomom");
  TTreeReaderValue<double> dca(tracks, "dca");
  TTreeReaderValue<double> dcaz(tracks, "dcaz");
  TTreeReaderValue<double> chi2(tracks, "chi2");
  TTreeReaderValue<int> towid(tracks, "towid");
  TTreeReaderValue<int> towexitid(tracks, "towexitid");
  TTreeReaderValue<int> gtowid(tracks, "gtowid");
  TTreeReaderValue<int> gtowexitid(tracks, "gtowexitid");
  TTreeReaderValue<double> ticurve(tracks, "icurve");
  TTreeReaderValue<double> tocurve(tracks, "ocurve");
  TTreeReaderValue<double> gicurve(tracks, "gicurve");
  TTreeReaderValue<double> gocurve(tracks, "gocurve");

  unsigned current_track = 0;
  while (events.Next()) {

    tracks.SetEntriesRange(current_track, current_track + *ntrack);
    while (tracks.Next()) {

      std::string name = dijetcore::MakeString(
          FLAGS_output_dir, "/track_id_", current_track, "_pt_",
          (*primary).Pt(), "_eta_", (*primary).Eta(), "_phi_", (*primary).Phi(),
          ".pdf");
      DrawExtrapolation(name, *tpi, *tpo);
      current_track++;
    }
  }

  out.Write();
  out.Close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}

void DrawExtrapolation(std::string name, std::vector<TVector3> &i,
                       std::vector<TVector3> &o) {
  TCanvas c("c", "c", 0, 0, 600, 600);
  TGeoManager manager("manager", "manager");
  TGeoMaterial mat("Al", 26.98, 13, 2.7);
  TGeoMedium med("MED", 1, &mat);
  TGeoVolume *vol = manager.MakeTube("TUBE", &med, 200, 200, 400);
  manager.SetTopVolume(vol);
  vol->SetLineWidth(2);
  manager.CloseGeometry();
  manager.SetNsegments(80);
  vol->Draw();
  // TView *view = gPad->GetView();
  // view->ShowAxis();

  c.SaveAs(name.c_str());
}