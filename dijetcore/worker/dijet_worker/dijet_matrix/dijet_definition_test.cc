#include "gtest/gtest.h"

#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"

TEST(DijetDefinition, Default) {
    dijetcore::DijetDefinition default_def;

    EXPECT_EQ(default_def.lead, nullptr);
    EXPECT_EQ(default_def.sub, nullptr);
    EXPECT_FALSE(default_def.isValid());
    EXPECT_FALSE(default_def.doMatching());
}

TEST(DijetDefinition, ValidDijetDef) {
    dijetcore::JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4);
    dijetcore::MatchDef* valid_matchdef1 = new dijetcore::MatchDef(valid_jetdef, valid_jetdef);
    dijetcore::MatchDef* valid_matchdef2 = new dijetcore::MatchDef(valid_jetdef, valid_jetdef);
    dijetcore::DijetDefinition valid_dijetdef(valid_matchdef1, valid_matchdef2);

    EXPECT_EQ(*valid_dijetdef.lead, *valid_matchdef1);
    EXPECT_EQ(*valid_dijetdef.sub, *valid_matchdef2);
    EXPECT_EQ(valid_dijetdef.lead, valid_matchdef1);
    EXPECT_EQ(valid_dijetdef.sub, valid_matchdef2);

    EXPECT_TRUE(valid_dijetdef.isValid());
    EXPECT_TRUE(valid_dijetdef.doMatching());
}

TEST(DijetDefinition, ValidDijetDefNoMatching) {
    dijetcore::JetDef valid_jetdef(fastjet::antikt_algorithm, 0.4);
    dijetcore::JetDef invalid_jetdef;
    dijetcore::MatchDef* valid_matchdef1 = new dijetcore::MatchDef(valid_jetdef, valid_jetdef);
    dijetcore::MatchDef* valid_matchdef2 = new dijetcore::MatchDef(valid_jetdef, invalid_jetdef);
    dijetcore::DijetDefinition valid_dijetdef(valid_matchdef1, valid_matchdef2);

    EXPECT_EQ(*valid_dijetdef.lead, *valid_matchdef1);
    EXPECT_EQ(*valid_dijetdef.sub, *valid_matchdef2);
    EXPECT_EQ(valid_dijetdef.lead, valid_matchdef1);
    EXPECT_EQ(valid_dijetdef.sub, valid_matchdef2);

    EXPECT_TRUE(valid_dijetdef.isValid());
    EXPECT_FALSE(valid_dijetdef.doMatching());

}
