// basic QA routine for embedding TStarJetPicoTrees

#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/trigger_lookup.h"

#include <set>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoReader.h"
#include "TStarJetPicoPrimaryTrack.h"

DIJETCORE_DEFINE_string(
    tree, "", "input root TFile (or txt file with TFile list) with JetTree(s)");
DIJETCORE_DEFINE_string(towList, "resources/bad_tower_lists/empty_list.txt",
                        "bad tower list, default is empty list");

int main(int argc, char* argv[]) {
  // setup command line flags
  dijetcore::SetUsageMessage("Generate a tree of all unique runids");
  dijetcore::ParseCommandLineFlags(&argc, argv);

  // initialize logging - search for logging related command line flags
  dijetcore::InitLogging(argv[0]);

  // build our input chain of ROOT tree(s) and init reader
  TChain* chain = dijetcore::NewChainFromInput(FLAGS_tree);
  TChain* mc_chain = dijetcore::NewChainFromInput(FLAGS_tree, "JetTreeMc");

  TStarJetPicoReader* reader = new TStarJetPicoReader();
  if (chain != nullptr) reader->SetInputChain(chain);
  reader->Init(-1);

  TStarJetPicoReader* mc_reader = new TStarJetPicoReader();
  if (chain != nullptr) mc_reader->SetInputChain(mc_chain);
  mc_reader->Init(-1);

  TH1D* dca = new TH1D("dca", "", 50, 0, 5);
  TH1D* mc_dca = new TH1D("mcdca", "", 50, 0, 5);

  for (int i = 0; i < reader->GetNOfEvents(); ++i) {
    reader->ReadEvent(i);
    mc_reader->ReadEvent(i);
    TClonesArray* r_array = reader->GetEvent()->GetPrimaryTracks();
    TClonesArray* mc_array = mc_reader->GetEvent()->GetPrimaryTracks();
    TIter r_iter(r_array);
    TIter mc_iter(mc_array);

    while (TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) r_iter.Next()) {
      dca->Fill(track->GetDCA());
    }
    while (TStarJetPicoPrimaryTrack* track = (TStarJetPicoPrimaryTrack*) mc_iter.Next()) {
      mc_dca->Fill(track->GetDCA());
    }
  }

  TCanvas c;
  dca->Scale(1.0/dca->Integral());
  mc_dca->Scale(1.0/mc_dca->Integral());
  
  dca->Draw();
  mc_dca->Draw();
  c.SaveAs("tmp.pdf");

  google::ShutdownGoogleLogging();
  google::ShutDownCommandLineFlags();
  return 0;
}