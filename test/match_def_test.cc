#include "jet_analysis/dijet_worker/dijet_matrix/match_def.hh"
#include "jet_analysis/dijet_worker/dijet_matrix/jet_def.hh"



int main() {
  
  MatchDef default_def;
  JetDef   default_jet_def;
  
  if (default_def.IsValid() ||
      default_def.CanMatch())
    return 1;
  
  if (default_def.InitialJetDef() != default_jet_def ||
      default_def.MatchedJetDef() != default_jet_def)
    return 1;
  
  fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
  fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
  fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);
  
  JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4,
                      area_def,
                      bkg_def,
                      area_def,
                      bkg_sel);
  
  MatchDef cant_match(valid_jetdef, JetDef());
  
  if (!cant_match.IsValid() ||
      cant_match.CanMatch() ||
      cant_match.dR() != 0.4)
    return 1;
  
  JetDef valid_jetdef_alt(fastjet::antikt_algorithm, 0.1,
                          area_def,
                          bkg_def,
                          area_def,
                          bkg_sel);
  
  MatchDef can_match(valid_jetdef, valid_jetdef_alt);
  
  if (!can_match.IsValid() ||
      !can_match.CanMatch() ||
      can_match.dR() != 0.1)
    return 1;
  
  if (can_match.InitialJetDef() != valid_jetdef ||
      can_match.MatchedJetDef() != valid_jetdef_alt)
    return 1;
  
  return 0;
}
