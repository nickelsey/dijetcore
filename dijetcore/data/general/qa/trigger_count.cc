// counts the number of triggers present in a jetpico file

#include "dijetcore/lib/filesystem.h"
#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/trigger_lookup.h"

#include <set>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoReader.h"

DIJETCORE_DEFINE_string(
    tree, "", "input root TFile (or txt file with TFile list) with JetTree(s)");
DIJETCORE_DEFINE_string(triggerset, "", "set of trigger ids to count in file");
DIJETCORE_DEFINE_string(towList, "resources/bad_tower_lists/empty_list.txt",
                        "bad tower list, default is empty list");

using dijetcore::dijetcore_map;

int main(int argc, char* argv[]) {
  // setup command line flags
  dijetcore::SetUsageMessage("Generate a tree of all unique runids");
  dijetcore::ParseCommandLineFlags(&argc, argv);

  // initialize logging - search for logging related command line flags
  dijetcore::InitLogging(argv[0]);

  // build our input chain of ROOT tree(s) and init reader
  TChain* chain = dijetcore::NewChainFromInput(FLAGS_tree);

  TStarJetPicoReader* reader = new TStarJetPicoReader();
  if (chain != nullptr) reader->SetInputChain(chain);
  reader->GetTowerCuts()->AddBadTowers(FLAGS_towList.c_str());
  reader->Init(-1);

  // get set of trigger IDs
  std::set<unsigned> triggerids = dijetcore::GetTriggerIDs(FLAGS_triggerset);

  // build counting structure
  dijetcore_map<unsigned, unsigned> trigger_counts;
  for (auto& id : triggerids) trigger_counts[id] = 0;

  while (reader->NextEvent()) {
    for (auto& trigger : triggerids) {
      if (reader->GetEvent()->GetHeader()->HasTriggerId(trigger))
        trigger_counts[trigger]++;
    }
  }

  // print results
  for (auto& entry : trigger_counts) {
    LOG(INFO) << "trigger: " << entry.first << " counts: " << entry.second;
  }

  google::ShutdownGoogleLogging();
  google::ShutDownCommandLineFlags();
  return 0;
}