#include "jet_analysis/dijet_worker/dijet_matrix/jet_def.hh"
#include "jet_analysis/dijet_worker/dijet_matrix/match_def.hh"
#include "jet_analysis/dijet_worker/dijet_matrix/dijet_definition.hh"
#include "jet_analysis/util/make_unique.h"

int main() {

  DijetDefinition default_def;
  
  if (default_def.lead != nullptr ||
      default_def.sub != nullptr)
    return 1;
  
  if (default_def.IsValid() != false ||
      default_def.DoMatching() != false)
    return 1;
  
  JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4);
  MatchDef* valid_matchdef1 = new MatchDef(valid_jetdef, valid_jetdef);
  MatchDef* valid_matchdef2 = new MatchDef(valid_jetdef, valid_jetdef);
  DijetDefinition valid_dijetdef(valid_matchdef1, valid_matchdef2);
  
  if (*valid_dijetdef.lead != *valid_matchdef1 ||
      *valid_dijetdef.sub != *valid_matchdef2)
    return 1;
  
  if (valid_dijetdef.IsValid() != true ||
      valid_dijetdef.DoMatching() != true)
    return 1;
  
  JetDef default_jetdef;
  MatchDef* valid_matchdef_no_matching = new MatchDef(valid_jetdef, default_jetdef);
  DijetDefinition valid_dijetdef_no_matching(valid_matchdef_no_matching, valid_matchdef_no_matching);
  if (valid_dijetdef_no_matching.IsValid() != true ||
      valid_dijetdef_no_matching.DoMatching() != false)
    return 1;
  
  DijetDefinition tmp(valid_dijetdef_no_matching);
  if (tmp != valid_dijetdef_no_matching)
    return 1;
  
  return 0;
}
