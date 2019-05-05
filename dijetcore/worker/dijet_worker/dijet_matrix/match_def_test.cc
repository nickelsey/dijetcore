#include "gtest/gtest.h"

#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"

TEST(MatchDef, DefaultSettings) {
    dijetcore::MatchDef default_matchdef;
    dijetcore::JetDef default_jetdef;

    EXPECT_FALSE(default_matchdef.isValid());
    EXPECT_FALSE(default_matchdef.canMatch());

    EXPECT_EQ(default_matchdef.initialJetDef(), default_jetdef);
    EXPECT_EQ(default_matchdef.matchedJetDef(), default_jetdef);
}

TEST(MatchDef, ValidMatchDef) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);
  
    dijetcore::JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4, area_def, bkg_def, area_def, bkg_sel);
  
    dijetcore::MatchDef cant_match(valid_jetdef, dijetcore::JetDef());
  
    EXPECT_TRUE(cant_match.isValid());
    EXPECT_FALSE(cant_match.canMatch());
    EXPECT_EQ(cant_match.initialJetDef(), valid_jetdef);
    EXPECT_EQ(cant_match.matchedJetDef(), dijetcore::JetDef());
}

TEST(MatchDef, MatchingMatchDef) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);
  
    dijetcore::JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4, area_def, bkg_def, area_def, bkg_sel);
    dijetcore::JetDef valid_jetdef_alt(fastjet::antikt_algorithm, 0.5, area_def, bkg_def, area_def, bkg_sel);
    dijetcore::MatchDef can_match(valid_jetdef, valid_jetdef_alt);
  
    EXPECT_TRUE(can_match.isValid());
    EXPECT_TRUE(can_match.canMatch());
    EXPECT_EQ(can_match.initialJetDef(), valid_jetdef);
    EXPECT_EQ(can_match.matchedJetDef(), valid_jetdef_alt);

}

TEST(MatchDef, EquivalentCluster) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);
  
    dijetcore::JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4, area_def, bkg_def, area_def, bkg_sel);
    dijetcore::JetDef valid_jetdef_alt(fastjet::antikt_algorithm, 0.5, area_def, bkg_def, area_def, bkg_sel);
    dijetcore::MatchDef can_match(valid_jetdef, valid_jetdef_alt);
    dijetcore::MatchDef can_match_alt(valid_jetdef_alt, valid_jetdef);

    EXPECT_TRUE(can_match.equivalentCluster(can_match));
    EXPECT_FALSE(can_match.equivalentCluster(can_match_alt));
}
