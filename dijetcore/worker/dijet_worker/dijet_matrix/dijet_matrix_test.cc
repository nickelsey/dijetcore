#include "gtest/gtest.h"

#include <string>

#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_matrix.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"
#include "dijetcore/util/fastjet/selector_compare.h"
bool CheckDijetDefinition(dijetcore::DijetDefinition def, fastjet::JetAlgorithm lead_alg, fastjet::JetAlgorithm sub_alg,
                          double lead_R, double sub_R, double lead_pt, double sub_pt, double lead_const_init_pt,
                          double lead_const_match_pt, double sub_const_init_pt, double sub_const_match_pt,
                          double const_eta, fastjet::RecombinationScheme scheme, fastjet::Strategy strategy,
                          fastjet::AreaType area_type, int ghost_repeat, double ghost_area, double grid_scatter,
                          double pt_scatter, double mean_ghost_pt, fastjet::JetDefinition bkg_def,
                          fastjet::AreaDefinition bkg_area_lead, fastjet::AreaDefinition bkg_area_sub) {
    dijetcore::MatchDef* lead = def.lead;
    dijetcore::MatchDef* sub = def.sub;
  
  // check the leading jet
  if (lead->InitialJetDef().R() != lead_R ||
      lead->MatchedJetDef().R() != lead_R ||
      lead->InitialJetDef().jet_algorithm() != lead_alg ||
      lead->MatchedJetDef().jet_algorithm() != lead_alg ||
      dijetcore::ExtractDoubleFromSelector(lead->InitialJetDef().ConstituentSelector(), "pt >=") != lead_const_init_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->MatchedJetDef().ConstituentSelector(), "pt >=") != lead_const_match_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->InitialJetDef().JetSelector(), "pt >=") != lead_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->MatchedJetDef().JetSelector(), "pt >=") != 0 ||
      lead->InitialJetDef().recombination_scheme() != scheme ||
      lead->MatchedJetDef().recombination_scheme() != scheme ||
      lead->InitialJetDef().strategy() != strategy ||
      lead->MatchedJetDef().strategy() != strategy ||
      lead->InitialJetDef().AreaDefinition().area_type() != area_type ||
      lead->MatchedJetDef().AreaDefinition().area_type() != area_type ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().ghost_maxrap() != (const_eta + lead_R) ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().ghost_maxrap() != (const_eta + lead_R) ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().repeat() != ghost_repeat ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().repeat() != ghost_repeat ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().ghost_area() != ghost_area ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().ghost_area() != ghost_area ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      lead->InitialJetDef().AreaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      lead->MatchedJetDef().AreaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      lead->InitialJetDef().BackgroundJetDef().R() != bkg_def.R() ||
      lead->MatchedJetDef().BackgroundJetDef().R() != bkg_def.R() ||
      lead->InitialJetDef().BackgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      lead->MatchedJetDef().BackgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      lead->InitialJetDef().BackgroundJetDef().strategy() != bkg_def.strategy() ||
      lead->MatchedJetDef().BackgroundJetDef().strategy() != bkg_def.strategy() ||
      lead->InitialJetDef().BackgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      lead->MatchedJetDef().BackgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      lead->InitialJetDef().BackgroundAreaDef().area_type() != bkg_area_lead.area_type() ||
      lead->MatchedJetDef().BackgroundAreaDef().area_type() != bkg_area_lead.area_type() ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_lead.ghost_spec().ghost_maxrap() ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_lead.ghost_spec().ghost_maxrap() ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().repeat() != bkg_area_lead.ghost_spec().repeat() ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().repeat() != bkg_area_lead.ghost_spec().repeat() ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().ghost_area() != bkg_area_lead.ghost_spec().ghost_area() ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().ghost_area() != bkg_area_lead.ghost_spec().ghost_area() ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_lead.ghost_spec().grid_scatter()  ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_lead.ghost_spec().grid_scatter()  ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_lead.ghost_spec().pt_scatter()  ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_lead.ghost_spec().pt_scatter()  ||
      lead->InitialJetDef().BackgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_lead.ghost_spec().mean_ghost_pt()  ||
      lead->MatchedJetDef().BackgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_lead.ghost_spec().mean_ghost_pt() )
    return false;
  
  // and check subleading jet
  if (sub->InitialJetDef().R() != sub_R ||
      sub->MatchedJetDef().R() != sub_R ||
      sub->InitialJetDef().jet_algorithm() != sub_alg ||
      sub->MatchedJetDef().jet_algorithm() != sub_alg ||
      dijetcore::ExtractDoubleFromSelector(sub->InitialJetDef().ConstituentSelector(), "pt >=") != sub_const_init_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->MatchedJetDef().ConstituentSelector(), "pt >=") != sub_const_match_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->InitialJetDef().JetSelector(), "pt >=") != sub_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->MatchedJetDef().JetSelector(), "pt >=") != 0 ||
      sub->InitialJetDef().recombination_scheme() != scheme ||
      sub->MatchedJetDef().recombination_scheme() != scheme ||
      sub->InitialJetDef().strategy() != strategy ||
      sub->MatchedJetDef().strategy() != strategy ||
      sub->InitialJetDef().AreaDefinition().area_type() != area_type ||
      sub->MatchedJetDef().AreaDefinition().area_type() != area_type ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().ghost_maxrap() != (const_eta + sub_R) ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().ghost_maxrap() != (const_eta + sub_R) ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().repeat() != ghost_repeat ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().repeat() != ghost_repeat ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().ghost_area() != ghost_area ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().ghost_area() != ghost_area ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      sub->InitialJetDef().AreaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      sub->MatchedJetDef().AreaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      sub->InitialJetDef().BackgroundJetDef().R() != bkg_def.R() ||
      sub->MatchedJetDef().BackgroundJetDef().R() != bkg_def.R() ||
      sub->InitialJetDef().BackgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      sub->MatchedJetDef().BackgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      sub->InitialJetDef().BackgroundJetDef().strategy() != bkg_def.strategy() ||
      sub->MatchedJetDef().BackgroundJetDef().strategy() != bkg_def.strategy() ||
      sub->InitialJetDef().BackgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      sub->MatchedJetDef().BackgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      sub->InitialJetDef().BackgroundAreaDef().area_type() != bkg_area_sub.area_type() ||
      sub->MatchedJetDef().BackgroundAreaDef().area_type() != bkg_area_sub.area_type() ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_sub.ghost_spec().ghost_maxrap() ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_sub.ghost_spec().ghost_maxrap() ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().repeat() != bkg_area_sub.ghost_spec().repeat() ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().repeat() != bkg_area_sub.ghost_spec().repeat() ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().ghost_area() != bkg_area_sub.ghost_spec().ghost_area() ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().ghost_area() != bkg_area_sub.ghost_spec().ghost_area() ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_sub.ghost_spec().grid_scatter()  ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_sub.ghost_spec().grid_scatter()  ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_sub.ghost_spec().pt_scatter()  ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_sub.ghost_spec().pt_scatter()  ||
      sub->InitialJetDef().BackgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_sub.ghost_spec().mean_ghost_pt()  ||
      sub->MatchedJetDef().BackgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_sub.ghost_spec().mean_ghost_pt() )
    return false;
  
  return true;
}

TEST(DijetMatrix, Default) {
    dijetcore::DijetMatrix default_matrix;
  
    default_matrix.Initialize();
    std::string expected_default_string = "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2";
  
    std::set<std::string> keys = default_matrix.Keys();
  
    EXPECT_NE(keys.find(expected_default_string), keys.end());
}


TEST(DijetMatrix, Initialization) {
    dijetcore::DijetMatrix matrix;
    matrix.ForceConstituentPtEquality(false);
    matrix.ForceConstituentEtaEquality(false);
    matrix.AddJetAlgorithm({fastjet::antikt_algorithm, fastjet::kt_algorithm});
    matrix.AddLeadJetR({0.4, 0.5});
    matrix.AddSubJetR({0.4, 0.5});
    matrix.AddLeadJetPt({20, 18});
    matrix.AddSubJetPt({10, 9});
    matrix.Initialize();
    
    EXPECT_EQ(matrix.Size(), 32);
    
    matrix.AddSubJetPt(40);
    matrix.Initialize();
    
    EXPECT_EQ(matrix.Size(), 32);

    matrix.AddLeadJetPt({2, 60});
    matrix.Initialize();
    
    EXPECT_EQ(matrix.Size(), 48 + 8);

    matrix.Clear();
    matrix.Initialize();
    
    EXPECT_EQ(matrix.Size(), 1);
}

TEST(DijetMatrix, ForceConstistuentEta) {
    dijetcore::DijetMatrix matrix;
    matrix.AddConstituentEta({1.0, 1.5});
    matrix.ForceConstituentEtaEquality(true);
    matrix.Initialize();
    
    EXPECT_EQ(matrix.Size(), 2);

    matrix.ForceConstituentEtaEquality(false);
    matrix.Initialize();
    
    EXPECT_EQ(matrix.Size(), 4);
}

TEST(DijetMatrix, ForceConstituentPt) {
    dijetcore::DijetMatrix matrix;

    matrix.ForceConstituentPtEquality(true);

    matrix.AddConstituentLeadInitialPt(2.0);
    matrix.AddConstituentSubInitialPt(2.1);

    matrix.Initialize();
    EXPECT_EQ(matrix.Size(), 0);

    matrix.ForceConstituentPtEquality(false);

    matrix.Initialize();
    EXPECT_EQ(matrix.Size(), 1);

    matrix.Clear();

    matrix.AddConstituentLeadInitialPt({2.0, 2.2});
    matrix.AddConstituentSubInitialPt({2.1, 2.2});
    matrix.Initialize();

    EXPECT_EQ(matrix.Size(), 4);

    matrix.ForceConstituentPtEquality(true);
    matrix.Initialize();

    EXPECT_EQ(matrix.Size(), 1);

    matrix.AddConstituentLeadMatchPt({0.2, 0.3});
    matrix.AddConstituentSubMatchPt(0.2);
    matrix.Initialize();

    EXPECT_EQ(matrix.Size(), 1);

    matrix.ForceConstituentPtEquality(false);
    matrix.Initialize();
    
    EXPECT_EQ(matrix.Size(), 8);

}

TEST(DijetMatrix, SingleCaseFullTest) {
    dijetcore::DijetMatrix matrix;
    matrix.ForceConstituentPtEquality(true);
    matrix.ForceConstituentEtaEquality(true);
    matrix.AddLeadJetPt(18);
    matrix.AddLeadJetR(0.4);
    matrix.AddLeadJetR(0.5);
    matrix.AddSubJetR(0.6);
    matrix.AddSubJetPt(9);
    matrix.AddConstituentEta(1.0);
    matrix.AddConstituentLeadInitialPt(2.5);
    matrix.AddConstituentLeadMatchPt(0.5);
    matrix.AddConstituentSubInitialPt(2.3);
    matrix.AddConstituentSubMatchPt(0.4);
    matrix.Initialize();

    EXPECT_EQ(matrix.Size(), 0);

    matrix.ForceConstituentPtEquality(false);
    matrix.ForceConstituentEtaEquality(false);
    matrix.Initialize();

    EXPECT_EQ(matrix.Size(), 2);

    for (auto key : matrix.Keys()) {
        double lead_R;
        if ((key.find("INIT_R_0.4") != std::string::npos) &&
            (key.find("MATCH_R_0.4") != std::string::npos)) {
            lead_R = 0.4;
        }
        else if ((key.find("INIT_R_0.5") != std::string::npos) &&
             (key.find("MATCH_R_0.5") != std::string::npos)) {
            lead_R = 0.5;
        }
        else {
            EXPECT_STREQ("bad key found in DijetMatrix", "");
        }
    
        // build the background area def
  
        fastjet::GhostedAreaSpec ghost_spec_lead(1.0 + 0.4, 1, 0.01, 1, 0.1, 1e-100);
        fastjet::AreaDefinition bkg_area_def_lead(fastjet::active_area_explicit_ghosts, ghost_spec_lead);
    
        fastjet::GhostedAreaSpec ghost_spec_sub(1.0 + 0.4, 1, 0.01, 1, 0.1, 1e-100);
        fastjet::AreaDefinition bkg_area_def_sub(fastjet::active_area_explicit_ghosts, ghost_spec_sub);
    
        EXPECT_TRUE(CheckDijetDefinition(*matrix.DijetDefinitions()[key], fastjet::antikt_algorithm, fastjet::antikt_algorithm,
                                      lead_R, 0.6, 18, 9, 2.5, 0.5, 2.3, 0.4, 1.0, fastjet::E_scheme, fastjet::Best,
                                      fastjet::active_area_explicit_ghosts, 1, 0.01, 1, 0.1, 1e-100,
                                      fastjet::JetDefinition(fastjet::kt_algorithm, 0.4), bkg_area_def_lead,
                                      bkg_area_def_sub));

    }
}

TEST(DijetMatrix, TestClusterSequenceSets) {
    dijetcore::DijetMatrix matrix;

    matrix.AddLeadJetPt({20.0, 18.0});
    matrix.AddSubJetPt({10.0, 9.0});
    matrix.AddLeadJetR(0.4);
    matrix.AddSubJetR(0.4);
    matrix.AddConstituentEta(0.4);
    matrix.AddConstituentLeadInitialPt(2.0);
    matrix.AddConstituentSubInitialPt(2.0);
    matrix.AddConstituentLeadMatchPt(0.2);
    matrix.AddConstituentSubMatchPt(0.2);
    matrix.Initialize();
    
    ASSERT_EQ(matrix.SortedDefinitions().size(), 1);
    EXPECT_EQ(matrix.SortedDefinitions()[0].size(), 4);

    matrix.AddLeadJetR(0.5);
    matrix.AddSubJetR(0.5);
    matrix.Initialize();

    ASSERT_EQ(matrix.SortedDefinitions().size(), 4);

    for (int i = 0; i < matrix.SortedDefinitions().size(); ++i) {
        EXPECT_EQ(matrix.SortedDefinitions()[i].size(), 4);
        
        auto& defs = matrix.SortedDefinitions()[i];

        std::set<double> lead_init_r;
        std::set<double> lead_match_r;
        std::set<double> sub_init_r;
        std::set<double> sub_match_r;

        for(auto& e : defs) {
            lead_init_r.insert(e.second->lead->InitialJetDef().R());
            lead_match_r.insert(e.second->lead->MatchedJetDef().R());
            sub_init_r.insert(e.second->sub->InitialJetDef().R());
            sub_match_r.insert(e.second->sub->MatchedJetDef().R());
        }

        EXPECT_EQ(lead_init_r.size(), 1);
        EXPECT_EQ(lead_match_r.size(), 1);
        EXPECT_EQ(sub_init_r.size(), 1);
        EXPECT_EQ(sub_match_r.size(), 1);
    }

    matrix.AddConstituentLeadInitialPt({1.0, 3.0});
    matrix.AddConstituentLeadMatchPt({0.3, 0.4});
    matrix.AddConstituentSubInitialPt({1.0, 3.0});
    matrix.AddConstituentSubMatchPt({0.3, 0.4});
    matrix.Initialize();

    for (int i = 0; i < matrix.SortedDefinitions().size(); ++i) {
        EXPECT_EQ(matrix.SortedDefinitions()[i].size(), 4);
   
        auto& defs = matrix.SortedDefinitions()[i];

        std::set<double> lead_init_r;
        std::set<double> lead_match_r;
        std::set<double> sub_init_r;
        std::set<double> sub_match_r;
        std::set<double> lead_init_const_pt;
        std::set<double> lead_match_const_pt;
        std::set<double> sub_init_const_pt;
        std::set<double> sub_match_const_pt;

        for(auto& e : defs) {
            lead_init_r.insert(e.second->lead->InitialJetDef().R());
            lead_match_r.insert(e.second->lead->MatchedJetDef().R());
            sub_init_r.insert(e.second->sub->InitialJetDef().R());
            sub_match_r.insert(e.second->sub->MatchedJetDef().R());
            lead_init_const_pt.insert(dijetcore::ExtractDoubleFromSelector(e.second->lead->InitialJetDef().ConstituentSelector(), "pt >="));
            lead_match_const_pt.insert(dijetcore::ExtractDoubleFromSelector(e.second->lead->MatchedJetDef().ConstituentSelector(), "pt >="));
            sub_init_const_pt.insert(dijetcore::ExtractDoubleFromSelector(e.second->sub->InitialJetDef().ConstituentSelector(), "pt >="));
            sub_match_const_pt.insert(dijetcore::ExtractDoubleFromSelector(e.second->sub->MatchedJetDef().ConstituentSelector(), "pt >="));
        }
        
        EXPECT_EQ(lead_init_r.size(), 1);
        EXPECT_EQ(lead_match_r.size(), 1);
        EXPECT_EQ(sub_init_r.size(), 1);
        EXPECT_EQ(sub_match_r.size(), 1);
        EXPECT_EQ(lead_init_const_pt.size(), 1);
        EXPECT_EQ(lead_match_const_pt.size(), 1);
        EXPECT_EQ(sub_init_const_pt.size(), 1);
        EXPECT_EQ(sub_match_const_pt.size(), 1);
    }
}


