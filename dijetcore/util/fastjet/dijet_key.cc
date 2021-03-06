#include "dijetcore/util/fastjet/dijet_key.h"

namespace dijetcore {
DijetKey ParseStringToDijetKey(string key) {
  DijetKey ret;

  if (!Consume(key, "LEAD_INIT_R_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_init_r = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "alg_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_init_alg = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "pt_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_init_pt = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "const_eta_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_init_const_eta = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "const_pt_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_init_const_pt = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "MATCH_R_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_match_r = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "alg_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_match_alg = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "pt_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_match_pt = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "const_eta_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_match_const_eta = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "const_pt_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.lead_match_const_pt = CastTo<double>(SplitOnNextOccurence(key, "_"));

  if (!Consume(key, "SUB_INIT_R_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_init_r = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "alg_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_init_alg = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "pt_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_init_pt = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "const_eta_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_init_const_eta = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "const_pt_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_init_const_pt = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "MATCH_R_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_match_r = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "alg_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_match_alg = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "pt_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_match_pt = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "const_eta_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_match_const_eta = CastTo<double>(SplitOnNextOccurence(key, "_"));
  if (!Consume(key, "const_pt_"))
    DIJETCORE_THROW("Could not parse key string");
  ret.sub_match_const_pt = CastTo<double>(key);

  return ret;
}

string ParseDijetKeyToString(DijetKey &key) {
  string ret;

  ret +=
      MakeString("LEAD_INIT_R_", key.lead_init_r, "_alg_", key.lead_init_alg);
  ret +=
      MakeString("_pt_", key.lead_init_pt, "_const_eta_",
                 key.lead_init_const_eta, "_const_pt_", key.lead_init_const_pt);
  ret += MakeString("_MATCH_R_", key.lead_match_r, "_alg_", key.lead_match_alg);
  ret += MakeString("_pt_", key.lead_match_pt, "_const_eta_",
                    key.lead_match_const_eta, "_const_pt_",
                    key.lead_match_const_pt);
  ret += MakeString("_SUB_INIT_R_", key.sub_init_r, "_alg_", key.sub_init_alg);
  ret +=
      MakeString("_pt_", key.sub_init_pt, "_const_eta_", key.sub_init_const_eta,
                 "_const_pt_", key.sub_init_const_pt);
  ret += MakeString("_MATCH_R_", key.sub_match_r, "_alg_", key.sub_match_alg);
  ret +=
      MakeString("_pt_", key.sub_match_pt, "_const_eta_",
                 key.sub_match_const_eta, "_const_pt_", key.sub_match_const_pt);

  return ret;
}

string MakeKeyFromDijetDefinition(const DijetDefinition &def) {
  string ret;

  MatchDef *lead_match = def.lead;
  MatchDef *sub_match = def.sub;

  // get all the values we will need
  string leadInitR, leadInitAlg, leadInitPt, subInitR, subInitAlg, subInitPt;
  string leadInitConstPt, leadInitConstEta, leadMatchConstPt, leadMatchConstEta;
  string subInitConstPt, subInitConstEta, subMatchConstPt, subMatchConstEta;
  string leadMatchR, leadMatchAlg, leadMatchPt, subMatchR, subMatchAlg,
      subMatchPt;

  // extract the values...
  std::stringstream stream;

  // first the lead jet
  stream << lead_match->initialJetDef().R();
  stream >> leadInitR;
  stream.clear();
  stream << lead_match->initialJetDef().jet_algorithm();
  stream >> leadInitAlg;
  stream.clear();
  stream << lead_match->matchedJetDef().R();
  stream >> leadMatchR;
  stream.clear();
  stream << lead_match->matchedJetDef().jet_algorithm();
  stream >> leadMatchAlg;
  stream.clear();
  stream << ExtractDoubleFromSelector(lead_match->initialJetDef().jetSelector(),
                                      "pt >=");
  stream >> leadInitPt;
  stream.clear();
  stream << ExtractDoubleFromSelector(lead_match->matchedJetDef().jetSelector(),
                                      "pt >=");
  stream >> leadMatchPt;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      lead_match->initialJetDef().constituentSelector(), "|eta| <=");
  stream >> leadInitConstEta;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      lead_match->matchedJetDef().constituentSelector(), "|eta| <=");
  stream >> leadMatchConstEta;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      lead_match->initialJetDef().constituentSelector(), "pt >=");
  stream >> leadInitConstPt;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      lead_match->matchedJetDef().constituentSelector(), "pt >=");
  stream >> leadMatchConstPt;
  stream.clear();

  // now the sublead
  stream << sub_match->initialJetDef().R();
  stream >> subInitR;
  stream.clear();
  stream << sub_match->initialJetDef().jet_algorithm();
  stream >> subInitAlg;
  stream.clear();
  stream << sub_match->matchedJetDef().R();
  stream >> subMatchR;
  stream.clear();
  stream << sub_match->matchedJetDef().jet_algorithm();
  stream >> subMatchAlg;
  stream.clear();
  stream << ExtractDoubleFromSelector(sub_match->initialJetDef().jetSelector(),
                                      "pt >=");
  stream >> subInitPt;
  stream.clear();
  stream << ExtractDoubleFromSelector(sub_match->matchedJetDef().jetSelector(),
                                      "pt >=");
  stream >> subMatchPt;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      sub_match->initialJetDef().constituentSelector(), "|eta| <=");
  stream >> subInitConstEta;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      sub_match->matchedJetDef().constituentSelector(), "|eta| <=");
  stream >> subMatchConstEta;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      sub_match->initialJetDef().constituentSelector(), "pt >=");
  stream >> subInitConstPt;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      sub_match->matchedJetDef().constituentSelector(), "pt >=");
  stream >> subMatchConstPt;
  stream.clear();

  ret += "LEAD_INIT_R_" + leadInitR + "_alg_" + leadInitAlg;
  ret += "_pt_" + leadInitPt + "_const_eta_" + leadInitConstEta + "_const_pt_" +
         leadInitConstPt;
  ret += "_MATCH_R_" + leadMatchR + "_alg_" + leadMatchAlg;
  ret += "_pt_" + leadMatchPt + "_const_eta_" + leadMatchConstEta +
         "_const_pt_" + leadMatchConstPt;
  ret += "_SUB_INIT_R_" + subInitR + "_alg_" + subInitAlg;
  ret += "_pt_" + subInitPt + "_const_eta_" + subInitConstEta + "_const_pt_" +
         subInitConstPt;
  ret += "_MATCH_R_" + subMatchR + "_alg_" + subMatchAlg;
  ret += "_pt_" + subMatchPt + "_const_eta_" + subMatchConstEta + "_const_pt_" +
         subMatchConstPt;

  return ret;
}

string MakeSortedKeyFromDijetDefinition(const DijetDefinition &def) {
  string ret;

  MatchDef *lead_match = def.lead;
  MatchDef *sub_match = def.sub;

  // get all the values we will need
  string leadInitR, leadInitAlg, subInitR, subInitAlg;
  string leadInitConstPt, leadInitConstEta, leadMatchConstPt, leadMatchConstEta;
  string subInitConstPt, subInitConstEta, subMatchConstPt, subMatchConstEta;
  string leadMatchR, leadMatchAlg, subMatchR, subMatchAlg;

  // extract the values...
  std::stringstream stream;

  // first the lead jet
  stream << lead_match->initialJetDef().R();
  stream >> leadInitR;
  stream.clear();
  stream << lead_match->initialJetDef().jet_algorithm();
  stream >> leadInitAlg;
  stream.clear();
  stream << lead_match->matchedJetDef().R();
  stream >> leadMatchR;
  stream.clear();
  stream << lead_match->matchedJetDef().jet_algorithm();
  stream >> leadMatchAlg;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      lead_match->initialJetDef().constituentSelector(), "|eta| <=");
  stream >> leadInitConstEta;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      lead_match->matchedJetDef().constituentSelector(), "|eta| <=");
  stream >> leadMatchConstEta;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      lead_match->initialJetDef().constituentSelector(), "pt >=");
  stream >> leadInitConstPt;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      lead_match->matchedJetDef().constituentSelector(), "pt >=");
  stream >> leadMatchConstPt;
  stream.clear();

  // now the sublead
  stream << sub_match->initialJetDef().R();
  stream >> subInitR;
  stream.clear();
  stream << sub_match->initialJetDef().jet_algorithm();
  stream >> subInitAlg;
  stream.clear();
  stream << sub_match->matchedJetDef().R();
  stream >> subMatchR;
  stream.clear();
  stream << sub_match->matchedJetDef().jet_algorithm();
  stream >> subMatchAlg;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      sub_match->initialJetDef().constituentSelector(), "|eta| <=");
  stream >> subInitConstEta;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      sub_match->matchedJetDef().constituentSelector(), "|eta| <=");
  stream >> subMatchConstEta;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      sub_match->initialJetDef().constituentSelector(), "pt >=");
  stream >> subInitConstPt;
  stream.clear();
  stream << ExtractDoubleFromSelector(
      sub_match->matchedJetDef().constituentSelector(), "pt >=");
  stream >> subMatchConstPt;
  stream.clear();

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