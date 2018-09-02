// reader_util_test.cc


#include <set>
#include <string>
#include <vector>

#include "jet_analysis/util/reader_util.hh"

#include "boost/filesystem.hpp"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetPicoUtils.h"

int main() {
  
  // first test the csv parser
  
  // possible locations to find the csv file...
  std::string csv_file_source = "${CMAKE_SOURCE_DIR}/jet_analysis/test/csv_test.csv";
  std::string csv_file_bin = "${CMAKE_BINARY_DIR}/jet_analysis/test/csv_test.csv";
  std::string csv_file_install = "${CMAKE_INSTALL_DIR}";
  
  if (!csv_file_install.empty())
    csv_file_install += "/jet_analysis/test/csv_test.csv";
  
  std::set<int> csv_out;
  if (boost::filesystem::exists(csv_file_source)) {
    csv_out = ParseCSV<int>(csv_file_source);
  }
  else if (boost::filesystem::exists(csv_file_bin)) {
    csv_out = ParseCSV<int>(csv_file_source);
  }
  else if (boost::filesystem::exists(csv_file_install)) {
    csv_out = ParseCSV<int>(csv_file_source);
  }
  else {
    std::cerr << "no copy of example csv file found, checked: \n" <<
    csv_file_source << "\n" <<
    csv_file_bin << "\n";
    if (csv_file_install.empty())
      std::cout << "(no install location found)" << std::endl;
    else
      std::cout << csv_file_install << "\n" << std::endl;
    return 1;
  }

  if (csv_out != std::set<int>{1, 54, 65, 34, 24, 65})
    return 1;
  
  
  // test function to strip leading characters from a string
  std::string test_string = "  = hello";
  StripLeadingChars(test_string, " =");
  if (test_string[0] != 'h') return 1;
  
  // find the bad tower list & bad run list files, plus the
  // parameter file used for non-default initialization
  std::vector<std::string> run_list_paths = {"${CMAKE_SOURCE_DIR}/jet_analysis/test/example_bad_run_list.txt",
                                             "${CMAKE_BINARY_DIR}/jet_analysis/test/example_bad_run_list.txt"};
  std::vector<std::string> tower_list_paths = {"${CMAKE_SOURCE_DIR}/jet_analysis/test/example_bad_tower_list.txt",
                                               "${CMAKE_BINARY_DIR}/jet_analysis/test/example_bad_tower_list.txt"};
  std::vector<std::string> param_list_paths = {"${CMAKE_SOURCE_DIR}/jet_analysis/test/reader_util_test.txt",
                                               "${CMAKE_BINARY_DIR}/jet_analysis/test/reader_util_test.txt"};
  if (!std::string("${CMAKE_INSTALL_DIR}").empty()) {
    run_list_paths.push_back("${CMAKE_INSTALL_DIR}/jet_analysis/test/example_bad_run_list.txt");
    tower_list_paths.push_back("${CMAKE_INSTALL_DIR}/jet_analysis/test/example_bad_tower_list.txt");
    param_list_paths.push_back("${CMAKE_INSTALL_DIR}/jet_analysis/test/reader_util_test.txt");
  }
  
  std::string run_list;
  std::string tower_list;
  std::string param_file;
  
  for (std::string path : run_list_paths) {
    if (boost::filesystem::exists(path)) {
      run_list = path;
      break;
    }
  }
  for (std::string path : tower_list_paths) {
    if (boost::filesystem::exists(path)) {
      tower_list = path;
      break;
    }
  }
  for (std::string path : param_list_paths) {
    if (boost::filesystem::exists(path)) {
      param_file = path;
      break;
    }
  }
  
  // test reader default initialization
  TChain* chain1 = new TChain("JetTree");
  TStarJetPicoReader* reader1 = new TStarJetPicoReader();
  InitReaderWithDefaults(reader1, chain1, tower_list, run_list);
  
  // test that default settings have been applied
  // tower cuts don't have access to private members or getter functions
  
  // so we start with track cuts
  if (reader1->GetTrackCuts()->GetDCACut() != 1.0 ||
      reader1->GetTrackCuts()->GetMinNFitPointsCut() != 20 ||
      reader1->GetTrackCuts()->GetFitOverMaxPointsCut() != 0.52 ||
      reader1->GetTrackCuts()->GetMaxPtCut() != 1000)
    return 1;
  
  // and event cuts
  if (reader1->GetEventCuts()->GetTriggerSelection() != "All" ||
      reader1->GetEventCuts()->GetVertexZCut() != 30.0 ||
      reader1->GetEventCuts()->GetVertexZDiffCut() != 3.0 ||
      reader1->GetEventCuts()->GetMaxEventPtCut() != 30.0 ||
      reader1->GetEventCuts()->GetMaxEventEtCut() != 30)
    return 1;
  
  // now test the the reader initialization by parameter file
  TChain* chain2 = new TChain("JetTree");
  TStarJetPicoReader* reader2 = new TStarJetPicoReader();
  
  InitReader(reader2, chain2, param_file);
  
  // check track cuts
  if (reader2->GetTrackCuts()->GetDCACut() != 1.2 ||
      reader2->GetTrackCuts()->GetMinNFitPointsCut() != 21 ||
      reader2->GetTrackCuts()->GetFitOverMaxPointsCut() != 0.50 ||
      reader2->GetTrackCuts()->GetMaxPtCut() != 1000)
    return 1;
  
  // check tower cuts
  if (reader2->GetTowerCuts()->GetTowerStatusRejection() != 1 ||
      reader2->GetTowerCuts()->Gety8PythiaCut() != 0 ||
      reader2->GetTowerCuts()->GetMaxEtCut() != 1500)
    return 1;

  // and check event cuts
  if (reader2->GetEventCuts()->GetTriggerSelection() != "all" ||
      reader2->GetEventCuts()->GetVertexZCut() != 30.0 ||
      reader2->GetEventCuts()->GetVertexZDiffCut() != 3.0 ||
      reader2->GetEventCuts()->GetMaxEventPtCut() != 31.0 ||
      reader2->GetEventCuts()->GetMaxEventEtCut() != 32.0 ||
      reader2->GetEventCuts()->GetMinEventEtCut() != (Float_t)5.1 || // conversion to float needed for success for some reason
      reader2->GetEventCuts()->GetRefMultCutMin() != 10 ||
      reader2->GetEventCuts()->GetRefMultCutMax() != 800 ||
      reader2->GetEventCuts()->GetRefCentMin() != 1 ||
      reader2->GetEventCuts()->GetRefCentMax() != 7 ||
      reader2->GetEventCuts()->GetPVRankingCut() != 1.5 ||
      reader2->GetEventCuts()->GetBbceCutMin() != 100 ||
      reader2->GetEventCuts()->GetBbceCutMax() != 1000)
    return 1;

  return 0;
}

