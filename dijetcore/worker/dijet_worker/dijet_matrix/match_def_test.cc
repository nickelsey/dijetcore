#include "gtest/gtest.h"

#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"

TEST(MatchDef, DefaultSettings) {
    dijetcore::MatchDef default_matchdef;
    dijetcore::JetDef default_jetdef;

    EXPECT_FALSE(default_matchdef.IsValid());
    EXPECT_FALSE(default_matchdef.CanMatch());

    EXPECT_EQ(default_matchdef.InitialJetDef(), default_jetdef);
    EXPECT_EQ(default_matchdef.MatchedJetDef(), default_jetdef);
}

TEST(MatchDef, ValidMatchDef) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);
  
    dijetcore::JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4, area_def, bkg_def, area_def, bkg_sel);
  
    dijetcore::MatchDef cant_match(valid_jetdef, dijetcore::JetDef());
  
    EXPECT_TRUE(cant_match.IsValid());
    EXPECT_FALSE(cant_match.CanMatch());
    EXPECT_EQ(cant_match.InitialJetDef(), valid_jetdef);
    EXPECT_EQ(cant_match.MatchedJetDef(), dijetcore::JetDef());
}

TEST(MatchDef, MatchingMatchDef) {
    fastjet::GhostedAreaSpec ghost_spec(0.6, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, ghost_spec);
    fastjet::JetDefinition bkg_def = fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
    fastjet::Selector bkg_sel = !fastjet::SelectorNHardest(2);
  
    dijetcore::JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4, area_def, bkg_def, area_def, bkg_sel);
    dijetcore::JetDef valid_jetdef_alt(fastjet::antikt_algorithm, 0.5, area_def, bkg_def, area_def, bkg_sel);
    dijetcore::MatchDef can_match(valid_jetdef, valid_jetdef_alt);
  
    EXPECT_TRUE(can_match.IsValid());
    EXPECT_TRUE(can_match.CanMatch());
    EXPECT_EQ(can_match.InitialJetDef(), valid_jetdef);
    EXPECT_EQ(can_match.MatchedJetDef(), valid_jetdef_alt);

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

    EXPECT_TRUE(can_match.EquivalentCluster(can_match));
    EXPECT_FALSE(can_match.EquivalentCluster(can_match_alt));
}
