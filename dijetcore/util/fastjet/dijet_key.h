#ifndef DIJETCORE_UTIL_FASTJET_DIJET_KEY_H
#define DIJETCORE_UTIL_FASTJET_DIJET_KEY_H

#include <string>
#include <sstream>
#include <iomanip>

#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"

#include "dijetcore/util/fastjet/selector_compare.h"
#include "dijetcore/lib/types.h"

namespace dijetcore {
  
  string MakeKeyFromDijetDefinition(const DijetDefinition& def) {
    string ret;
    
    MatchDef* lead_match = def.lead;
    MatchDef* sub_match = def.sub;
    
    // get all the values we will need
    string leadInitR, leadInitAlg, leadInitPt, subInitR, subInitAlg, subInitPt;
    string leadInitConstPt, leadInitConstEta, leadMatchConstPt, leadMatchConstEta;
    string subInitConstPt, subInitConstEta, subMatchConstPt, subMatchConstEta;
    string leadMatchR, leadMatchAlg, leadMatchPt, subMatchR, subMatchAlg, subMatchPt;
    
    // extract the values...
    std::stringstream stream;
    
    // first the lead jet
    stream << lead_match->InitialJetDef().R(); stream >> leadInitR; stream.clear();
    stream << lead_match->InitialJetDef().jet_algorithm(); stream >> leadInitAlg; stream.clear();
    stream << lead_match->MatchedJetDef().R(); stream >> leadMatchR; stream.clear();
    stream << lead_match->MatchedJetDef().jet_algorithm(); stream >> leadMatchAlg; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->InitialJetDef().JetSelector(), "pt >="); stream >> leadInitPt; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->MatchedJetDef().JetSelector(), "pt >="); stream >> leadMatchPt; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->InitialJetDef().ConstituentSelector(), "|rap| <="); stream >> leadInitConstEta; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->MatchedJetDef().ConstituentSelector(), "|rap| <="); stream >> leadMatchConstEta; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->InitialJetDef().ConstituentSelector(), "pt >="); stream >> leadInitConstPt; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->MatchedJetDef().ConstituentSelector(), "pt >="); stream >> leadMatchConstPt; stream.clear();
    
    // now the sublead
    stream << sub_match->InitialJetDef().R(); stream >> subInitR; stream.clear();
    stream << sub_match->InitialJetDef().jet_algorithm(); stream >> subInitAlg; stream.clear();
    stream << sub_match->MatchedJetDef().R(); stream >> subMatchR; stream.clear();
    stream << sub_match->MatchedJetDef().jet_algorithm(); stream >> subMatchAlg; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->InitialJetDef().JetSelector(), "pt >="); stream >> subInitPt; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->MatchedJetDef().JetSelector(), "pt >="); stream >> subMatchPt; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->InitialJetDef().ConstituentSelector(), "|rap| <="); stream >> subInitConstEta; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->MatchedJetDef().ConstituentSelector(), "|rap| <="); stream >> subMatchConstEta; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->InitialJetDef().ConstituentSelector(), "pt >="); stream >> subInitConstPt; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->MatchedJetDef().ConstituentSelector(), "pt >="); stream >> subMatchConstPt; stream.clear();
    
    ret += "LEAD_INIT_R_" + leadInitR + "_alg_" + leadInitAlg;
    ret += "_pt_" + leadInitPt + "_const_eta_" + leadInitConstEta + "_const_pt_" + leadInitConstPt;
    ret += "_MATCH_R_" + leadMatchR + "_alg_" + leadMatchAlg;
    ret += "_pt_" + leadMatchPt + "_const_eta_" + leadMatchConstEta + "_const_pt_" + leadMatchConstPt;
    ret += "_SUB_INIT_R_" + subInitR + "_alg_" + subInitAlg;
    ret += "_pt_" + subInitPt + "_const_eta_" + subInitConstEta + "_const_pt_" + subInitConstPt;
    ret += "_MATCH_R_" + subMatchR + "_alg_" + subMatchAlg;
    ret += "_pt_" + subMatchPt + "_const_eta_" + subMatchConstEta + "_const_pt_" + subMatchConstPt;
    
    return ret;
  }
  
  string MakeSortedKeyFromDijetDefinition(const DijetDefinition& def) {
    string ret;
    
    MatchDef* lead_match = def.lead;
    MatchDef* sub_match = def.sub;
    
    // get all the values we will need
    string leadInitR, leadInitAlg, subInitR, subInitAlg;
    string leadInitConstPt, leadInitConstEta, leadMatchConstPt, leadMatchConstEta;
    string subInitConstPt, subInitConstEta, subMatchConstPt, subMatchConstEta;
    string leadMatchR, leadMatchAlg, subMatchR, subMatchAlg;
    
    // extract the values...
    std::stringstream stream;
    
    // first the lead jet
    stream << lead_match->InitialJetDef().R(); stream >> leadInitR; stream.clear();
    stream << lead_match->InitialJetDef().jet_algorithm(); stream >> leadInitAlg; stream.clear();
    stream << lead_match->MatchedJetDef().R(); stream >> leadMatchR; stream.clear();
    stream << lead_match->MatchedJetDef().jet_algorithm(); stream >> leadMatchAlg; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->InitialJetDef().ConstituentSelector(), "|rap| <="); stream >> leadInitConstEta; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->MatchedJetDef().ConstituentSelector(), "|rap| <="); stream >> leadMatchConstEta; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->InitialJetDef().ConstituentSelector(), "pt >="); stream >> leadInitConstPt; stream.clear();
    stream << ExtractDoubleFromSelector(lead_match->MatchedJetDef().ConstituentSelector(), "pt >="); stream >> leadMatchConstPt; stream.clear();
    
    // now the sublead
    stream << sub_match->InitialJetDef().R(); stream >> subInitR; stream.clear();
    stream << sub_match->InitialJetDef().jet_algorithm(); stream >> subInitAlg; stream.clear();
    stream << sub_match->MatchedJetDef().R(); stream >> subMatchR; stream.clear();
    stream << sub_match->MatchedJetDef().jet_algorithm(); stream >> subMatchAlg; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->InitialJetDef().ConstituentSelector(), "|rap| <="); stream >> subInitConstEta; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->MatchedJetDef().ConstituentSelector(), "|rap| <="); stream >> subMatchConstEta; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->InitialJetDef().ConstituentSelector(), "pt >="); stream >> subInitConstPt; stream.clear();
    stream << ExtractDoubleFromSelector(sub_match->MatchedJetDef().ConstituentSelector(), "pt >="); stream >> subMatchConstPt; stream.clear();
    
    ret += "LEAD_INIT_R_" + leadInitR + "_alg_" + leadInitAlg;
    ret += "_const_eta_" + leadInitConstEta + "_const_pt_" + leadInitConstPt;
    ret += "_MATCH_R_" + leadMatchR + "_alg_" + leadMatchAlg;
    ret += "_const_eta_" + leadMatchConstEta + "_const_pt_" + leadMatchConstPt;
    ret += "_SUB_INIT_R_" + subInitR + "_alg_" + subInitAlg;
    ret += "_const_eta_" + subInitConstEta + "_const_pt_" + subInitConstPt;
    ret += "_MATCH_R_" + subMatchR + "_alg_" + subMatchAlg;
    ret += "_const_eta_" + subMatchConstEta + "_const_pt_" + subMatchConstPt;
    
    return ret;
  }
  
} // namespace dijetcore

#endif // DIJETCORE_LIB_FASTJET_DIJET_KEY_H
