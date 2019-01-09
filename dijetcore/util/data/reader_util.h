#ifndef DIJETCORE_UTIL_DATA_READER_UTIL_H
#define DIJETCORE_UTIL_DATA_READER_UTIL_H

#include <string>
#include <set>
#include <fstream>
#include <algorithm>

#include "TChain.h"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetPicoUtils.h"

#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/assert.h"
#include "dijetcore/lib/logging.h"

#include "boost/filesystem.hpp"

namespace dijetcore {
  
  bool GetNextValidEvent(TStarJetPicoReader* reader, std::set<unsigned> triggers) {
    // loop over events to find one with an accepted trigger
    while (reader->NextEvent()) {
      if (triggers.size() == 0)
        return true;
      for (auto trigger : triggers)
        if (reader->GetEvent()->GetHeader()->HasTriggerId(trigger))
          return true;
    }
    
    // if none found, reset to first event and try again
    reader->ReadEvent(0);
    do {
      if (triggers.size() == 0)
        return true;
      for (auto trigger : triggers)
        if (reader->GetEvent()->GetHeader()->HasTriggerId(trigger))
          return true;
    } while (reader->NextEvent());
    // if we still haven't found an acceptable event, then
    // there is either a problem with the reader or data,
    // or the trigger is not present in the current data, exit
    return false;
  }
  
  // builds a TChain ready to use with the
  // TStarJetPicoReader, if unknown input format,
  // returns nullptr
  // supports .root, .list, and .txt
  TChain* NewChainFromInput(std::string str, int seed = 0, int eventThreshold = 0) {
    TChain* chain = nullptr;
    
    if (str.empty())
      return chain;
    
    if (HasEnding(str, ".root")) {
      chain = new TChain("JetTree");
      chain->Add(str.c_str());
    }
    else if (HasEnding(str, ".list") || HasEnding(str, ".txt")) {
      
      std::vector<std::string> files;
      // first we will read in each line and make sure its a valid file
      std::ifstream stream(str);
      if (stream.is_open()) {
        std::string line;
        while (std::getline(stream, line)) {
          if (boost::filesystem::exists(line))
            files.push_back(line);
          else {
            LOG(ERROR) << "input file: " << line << "does not exist, skipping";
          }
        }
        
        // now shuffle the files and add to the chain
        chain = new TChain("JetTree");
        std::mt19937 g(seed);
        std::shuffle(files.begin(), files.end(), g);
        for (auto& file : files) {
          if (eventThreshold > 0 &&
              chain->GetEntries() > eventThreshold)
            break;
          LOG(INFO) << "Adding file: " << file << " to chain";
          chain->Add(file.c_str());
        }
      }
    }
    return chain;
  }
  
  void StripLeadingChars(std::string& line, std::string chars) {
    while (chars.find(line[0]) != std::string::npos)
      Consume(line, line[0]);
  }
  
  // parse a CSV file to a set of unique entries.
  // All comments must start on their own line, and be proceeded
  // by a pound sign (#)
  template <typename T>
  std::set<T> ParseCSV(std::string csv) {
    // return set
    std::set<T> ret;
    
    std::ifstream fs(csv);
    std::string line;
    // first, split by line
    while (std::getline(fs, line)) {
      if (line.size() == 0) // reject empty lines
        continue;
      if (line[0] == '#') // reject comments
        continue;
      // split the string by commas
      std::istringstream ss(line);
      while (ss) {
        std::string str_value;
        std::getline(ss, str_value, ',');
        if (CanCast<T>(str_value)) {
          ret.insert(CastTo<T>(str_value));
        }
      }
    }
    return ret;
  }
  
  // initialize the reader with decent defaults,
  // plus an user defined hot tower list,
  // and an optional bad run list
  void InitReaderWithDefaults(TStarJetPicoReader* reader,
                              TChain* chain,
                              std::string bad_tower_list = "",
                              std::string bad_run_list = "") {
    
    // add the chain
    if (chain != nullptr)
      reader->SetInputChain(chain);
    
    // set hadronic correction
    reader->SetApplyFractionHadronicCorrection( true );
    reader->SetFractionHadronicCorrection( 0.999 );
    reader->SetRejectTowerElectrons( kFALSE );
    
    // tower cuts
    reader->GetTowerCuts()->SetMaxEtCut(1000);            // essentially infinity - cut in eventcuts
    
    // track cuts
    reader->GetTrackCuts()->SetDCACut(1.0);                // distance of closest approach to primary vtx
    reader->GetTrackCuts()->SetMinNFitPointsCut(20);       // minimum fit points in track reco
    reader->GetTrackCuts()->SetFitOverMaxPointsCut(0.52);  // minimum ratio of fit points used over possible
    reader->GetTrackCuts()->SetMaxPtCut(1000);             // essentially infinity - cut in eventcuts
    
    // event cuts
    reader->GetEventCuts()->SetMaxEventPtCut(30);          // Set Maximum track Pt
    reader->GetEventCuts()->SetMaxEventEtCut(30);          // Set Maximum tower Et
    reader->GetEventCuts()->SetVertexZCut(30);             // vertex z range (z = beam axis)
    reader->GetEventCuts()->SetTriggerSelection("All");    // setting trigger selection - set to all, selected later
    reader->GetEventCuts()->SetVertexZDiffCut(3);          // cut on Vz - VPD Vz
    
    // if a bad tower list is specified, add to tower cuts
    if (!bad_tower_list.empty())
      reader->GetTowerCuts()->AddBadTowers(bad_tower_list.c_str());
    
    // if bad run list is specified, add to reader
    if (!bad_run_list.empty()) {
      std::set<int> bad_runs = ParseCSV<int>(bad_run_list);
      for (auto run : bad_runs) {
        reader->AddMaskedRun(run);
      }
    }
    
    reader->Init();
  }
  
  void InitReader(TStarJetPicoReader* reader, TChain* chain, std::string file,
                  std::string bad_tower_list = "", std::string bad_run_list = "") {
    
    // first, set the input chain
    reader->SetInputChain(chain);
    
    // setup the file stream
    std::ifstream fs(file);
    std::string line;
    
    // first, split by line
    while (std::getline(fs, line)) {
      
      // ignore comments and empty lines
      if (line[0] == '#') continue;
      if (line.empty())   continue;
      
      // lowercase the string...
      std::transform(line.begin(), line.end(), line.begin(), ::tolower);
      
      // keep an unmodified copy
      std::string tmp = line;
      
      // now split by which type of cut it is:
      // badrunlist
      // mip/hadronic correction
      // eventcuts
      // towercuts
      // trackcuts
      if (Consume(line,"badrunlist")) {
        StripLeadingChars(line, " =");
        std::set<int> runs = ParseCSV<int>(line);
        for (auto run : runs) {
          reader->AddMaskedRun(run);
        }
      }
      else if (Consume(line, "hadcorrfraction")) {
        StripLeadingChars(line, " =");
        reader->SetFractionHadronicCorrection(atof(line.c_str()));
      }
      else if (Consume(line, "hadroniccorrection")) {
        StripLeadingChars(line, " =");
        reader->SetApplyFractionHadronicCorrection(atoi(line.c_str()));
      }
      else if (Consume(line, "mip")) {
        StripLeadingChars(line, " =");
        reader->SetApplyMIPCorrection(atoi(line.c_str()));
      }
      else if (Consume(line, "event::")) {
        
        if (Consume(line, "trigger")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetTriggerSelection(line.c_str());
        }
        else if (Consume(line, "vzdif")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetVertexZDiffCut(atof(line.c_str()));
        }
        else if (Consume(line, "vz")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetVertexZCut(atof(line.c_str()));
        }
        else if (Consume(line, "vertexz")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetVertexZCut(atof(line.c_str()));
        }
        else if (Consume(line, "refmultmax")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetRefMultCut(reader->GetEventCuts()->GetRefMultCutMin(),
                                                atof(line.c_str()));
        }
        else if (Consume(line, "refmultmin")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetRefMultCut(atof(line.c_str()),
                                                reader->GetEventCuts()->GetRefMultCutMax());
        }
        else if (Consume(line, "bbcemin")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetBbceCut(atoi(line.c_str()),
                                             reader->GetEventCuts()->GetBbceCutMax());
        }
        else if (Consume(line, "bbcemax")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetBbceCut(reader->GetEventCuts()->GetBbceCutMin(),
                                             atoi(line.c_str()));
        }
        else if (Consume(line, "refcentmin")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetReferenceCentralityCut(atoi(line.c_str()),
                                                            reader->GetEventCuts()->GetRefCentMax());
        }
        else if (Consume(line, "refcentmax")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetReferenceCentralityCut(reader->GetEventCuts()->GetRefCentMin(),
                                                            atoi(line.c_str()));
        }
        else if (Consume(line, "pvrank")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetPVRankingCut(atof(line.c_str()));
        }
        else if (Consume(line, "ptmax")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetMaxEventPtCut(atof(line.c_str()));
        }
        else if (Consume(line, "etmin")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetMinEventEtCut(atof(line.c_str()));
        }
        else if (Consume(line, "etmax")) {
          StripLeadingChars(line, " =");
          reader->GetEventCuts()->SetMaxEventEtCut(atof(line.c_str()));
        }
        else {
          DIJETCORE_THROW("unrecognized configuration line: ", tmp);
        }
      }
      else if (Consume(line, "track::")) {
        
        if (Consume(line, "dca")) {
          StripLeadingChars(line, " =");
          reader->GetTrackCuts()->SetDCACut(atof(line.c_str()));
        }
        else if (Consume(line, "minfitpoint")) {
          StripLeadingChars(line, " =");
          reader->GetTrackCuts()->SetMinNFitPointsCut(atoi(line.c_str()));
        }
        else if (Consume(line, "minfitfrac")) {
          StripLeadingChars(line, " =");
          reader->GetTrackCuts()->SetFitOverMaxPointsCut(atof(line.c_str()));
        }
        else if (Consume(line, "ptmax")) {
          StripLeadingChars(line, " =");
          reader->GetTrackCuts()->SetMaxPtCut(atof(line.c_str()));
        }
        else if (Consume(line, "chi2")) {
          StripLeadingChars(line, " =");
          reader->GetTrackCuts()->SetMaxChi2Cut(atof(line.c_str()));
        }
        else if (Consume(line, "pct")) {
          StripLeadingChars(line, " =");
          reader->GetTrackCuts()->SetPCTCut(atoi(line.c_str()));
        }
        else if (Consume(line, "phirange")) {
          StripLeadingChars(line, " =");
          std::stringstream ss(line);
          std::vector<std::string> vec;
          std::string segment;
          while(std::getline(ss, segment, '-')) {
            vec.push_back(segment);
          }
          if (vec.size() != 2) {
            DIJETCORE_THROW("unrecognized format for track phi range");
          }
          reader->GetTrackCuts()->RestrictPhiRange(atof(vec[0].c_str()), atof(vec[1].c_str()));
        }
        else {
          DIJETCORE_THROW("unrecognized configuration line: ", tmp);
        }
      }
      else if (Consume(line, "tower::")) {
        
        if (Consume(line, "towerstatus")) {
          StripLeadingChars(line, " =");
          reader->GetTowerCuts()->UseTowerStatusRejection(atoi(line.c_str()));
        }
        else if (Consume(line, "badtowerlist")) {
          StripLeadingChars(line, " =");
          reader->GetTowerCuts()->AddBadTowers(line.c_str());
        }
        else if (Consume(line, "y8pythia")) {
          StripLeadingChars(line, " =");
          reader->GetTowerCuts()->Sety8PythiaCut(atoi(line.c_str()));
        }
        else if (Consume(line, "etmax")) {
          StripLeadingChars(line, " =");
          reader->GetTowerCuts()->SetMaxEtCut(atof(line.c_str()));
        }
        else if (Consume(line, "phirange")) {
          StripLeadingChars(line, " =");
          std::stringstream ss(line);
          std::vector<std::string> vec;
          std::string segment;
          while(std::getline(ss, segment, '-')) {
            vec.push_back(segment);
          }
          if (vec.size() != 2) {
            DIJETCORE_THROW("unrecognized format for tower phi range");
          }
          reader->GetTowerCuts()->RestrictPhiRange(atof(vec[0].c_str()), atof(vec[1].c_str()));
        }
        else {
          DIJETCORE_THROW("unrecognized configuration line: ", tmp);
        }
      }
      else {
        DIJETCORE_THROW("unrecognized configuration line: ", line);
      }
      
    }
    
    // if a bad tower list is specified, add to tower cuts
    if (!bad_tower_list.empty())
      reader->GetTowerCuts()->AddBadTowers(bad_tower_list.c_str());
    // if bad run list is specified, add to reader
    if (!bad_run_list.empty()) {
      std::set<int> bad_runs = ParseCSV<int>(bad_run_list);
      for (auto run : bad_runs) {
        reader->AddMaskedRun(run);
      }
    }
    reader->Init();
  }
  
} // namespace dijetcore

#endif // DIJETCORE_UTIL_DATA_READER_UTIL_H
