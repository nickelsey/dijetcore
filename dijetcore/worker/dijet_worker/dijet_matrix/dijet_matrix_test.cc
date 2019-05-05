#include "gtest/gtest.h"

#include <string>

#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_matrix.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"
#include "dijetcore/util/fastjet/selector_compare.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/lib/logging.h"

bool CheckDijetDefinition(dijetcore::DijetDefinition def, fastjet::JetAlgorithm lead_alg, fastjet::JetAlgorithm sub_alg,
                          double lead_R, double lead_match_R, double sub_R, double sub_match_R, double lead_pt, double sub_pt,
                          double lead_const_init_pt, double lead_const_match_pt, double sub_const_init_pt, double sub_const_match_pt,
                          double const_eta, fastjet::RecombinationScheme scheme, fastjet::Strategy strategy,
                          fastjet::AreaType area_type, int ghost_repeat, double ghost_area, double grid_scatter, double pt_scatter,
                          double mean_ghost_pt, fastjet::JetDefinition bkg_def, fastjet::AreaDefinition bkg_area_lead,
                          fastjet::AreaDefinition bkg_area_sub) {
    dijetcore::MatchDef* lead = def.lead;
    dijetcore::MatchDef* sub = def.sub;
  
  // check the leading jet
  if (lead->initialJetDef().R() != lead_R ||
      lead->matchedJetDef().R() != lead_match_R ||
      lead->initialJetDef().jet_algorithm() != lead_alg ||
      lead->matchedJetDef().jet_algorithm() != lead_alg ||
      dijetcore::ExtractDoubleFromSelector(lead->initialJetDef().constituentSelector(), "pt >=") != lead_const_init_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->matchedJetDef().constituentSelector(), "pt >=") != lead_const_match_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->initialJetDef().jetSelector(), "pt >=") != lead_pt ||
      dijetcore::ExtractDoubleFromSelector(lead->matchedJetDef().jetSelector(), "pt >=") != 0 ||
      lead->initialJetDef().recombination_scheme() != scheme ||
      lead->matchedJetDef().recombination_scheme() != scheme ||
      lead->initialJetDef().strategy() != strategy ||
      lead->matchedJetDef().strategy() != strategy ||
      lead->initialJetDef().areaDefinition().area_type() != area_type ||
      lead->matchedJetDef().areaDefinition().area_type() != area_type ||
      lead->initialJetDef().areaDefinition().ghost_spec().ghost_maxrap() != (const_eta + lead_R) ||
      lead->matchedJetDef().areaDefinition().ghost_spec().ghost_maxrap() != (const_eta + lead_R) ||
      lead->initialJetDef().areaDefinition().ghost_spec().repeat() != ghost_repeat ||
      lead->matchedJetDef().areaDefinition().ghost_spec().repeat() != ghost_repeat ||
      lead->initialJetDef().areaDefinition().ghost_spec().ghost_area() != ghost_area ||
      lead->matchedJetDef().areaDefinition().ghost_spec().ghost_area() != ghost_area ||
      lead->initialJetDef().areaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      lead->matchedJetDef().areaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      lead->initialJetDef().areaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      lead->matchedJetDef().areaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      lead->initialJetDef().areaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      lead->matchedJetDef().areaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      lead->initialJetDef().backgroundJetDef().R() != bkg_def.R() ||
      lead->matchedJetDef().backgroundJetDef().R() != bkg_def.R() ||
      lead->initialJetDef().backgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      lead->matchedJetDef().backgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      lead->initialJetDef().backgroundJetDef().strategy() != bkg_def.strategy() ||
      lead->matchedJetDef().backgroundJetDef().strategy() != bkg_def.strategy() ||
      lead->initialJetDef().backgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      lead->matchedJetDef().backgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      lead->initialJetDef().backgroundAreaDef().area_type() != bkg_area_lead.area_type() ||
      lead->matchedJetDef().backgroundAreaDef().area_type() != bkg_area_lead.area_type() ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_lead.ghost_spec().ghost_maxrap() ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_lead.ghost_spec().ghost_maxrap() ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().repeat() != bkg_area_lead.ghost_spec().repeat() ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().repeat() != bkg_area_lead.ghost_spec().repeat() ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().ghost_area() != bkg_area_lead.ghost_spec().ghost_area() ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().ghost_area() != bkg_area_lead.ghost_spec().ghost_area() ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_lead.ghost_spec().grid_scatter()  ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_lead.ghost_spec().grid_scatter()  ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_lead.ghost_spec().pt_scatter()  ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_lead.ghost_spec().pt_scatter()  ||
      lead->initialJetDef().backgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_lead.ghost_spec().mean_ghost_pt()  ||
      lead->matchedJetDef().backgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_lead.ghost_spec().mean_ghost_pt() )
    return false;
  
  // and check subleading jet
  if (sub->initialJetDef().R() != sub_R ||
      sub->matchedJetDef().R() != sub_match_R ||
      sub->initialJetDef().jet_algorithm() != sub_alg ||
      sub->matchedJetDef().jet_algorithm() != sub_alg ||
      dijetcore::ExtractDoubleFromSelector(sub->initialJetDef().constituentSelector(), "pt >=") != sub_const_init_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->matchedJetDef().constituentSelector(), "pt >=") != sub_const_match_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->initialJetDef().jetSelector(), "pt >=") != sub_pt ||
      dijetcore::ExtractDoubleFromSelector(sub->matchedJetDef().jetSelector(), "pt >=") != 0 ||
      sub->initialJetDef().recombination_scheme() != scheme ||
      sub->matchedJetDef().recombination_scheme() != scheme ||
      sub->initialJetDef().strategy() != strategy ||
      sub->matchedJetDef().strategy() != strategy ||
      sub->initialJetDef().areaDefinition().area_type() != area_type ||
      sub->matchedJetDef().areaDefinition().area_type() != area_type ||
      sub->initialJetDef().areaDefinition().ghost_spec().ghost_maxrap() != (const_eta + sub_R) ||
      sub->matchedJetDef().areaDefinition().ghost_spec().ghost_maxrap() != (const_eta + sub_R) ||
      sub->initialJetDef().areaDefinition().ghost_spec().repeat() != ghost_repeat ||
      sub->matchedJetDef().areaDefinition().ghost_spec().repeat() != ghost_repeat ||
      sub->initialJetDef().areaDefinition().ghost_spec().ghost_area() != ghost_area ||
      sub->matchedJetDef().areaDefinition().ghost_spec().ghost_area() != ghost_area ||
      sub->initialJetDef().areaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      sub->matchedJetDef().areaDefinition().ghost_spec().grid_scatter() != grid_scatter ||
      sub->initialJetDef().areaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      sub->matchedJetDef().areaDefinition().ghost_spec().pt_scatter() != pt_scatter ||
      sub->initialJetDef().areaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      sub->matchedJetDef().areaDefinition().ghost_spec().mean_ghost_pt() != mean_ghost_pt ||
      sub->initialJetDef().backgroundJetDef().R() != bkg_def.R() ||
      sub->matchedJetDef().backgroundJetDef().R() != bkg_def.R() ||
      sub->initialJetDef().backgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      sub->matchedJetDef().backgroundJetDef().jet_algorithm() != bkg_def.jet_algorithm() ||
      sub->initialJetDef().backgroundJetDef().strategy() != bkg_def.strategy() ||
      sub->matchedJetDef().backgroundJetDef().strategy() != bkg_def.strategy() ||
      sub->initialJetDef().backgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      sub->matchedJetDef().backgroundJetDef().recombination_scheme() != bkg_def.recombination_scheme() ||
      sub->initialJetDef().backgroundAreaDef().area_type() != bkg_area_sub.area_type() ||
      sub->matchedJetDef().backgroundAreaDef().area_type() != bkg_area_sub.area_type() ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_sub.ghost_spec().ghost_maxrap() ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().ghost_maxrap() != bkg_area_sub.ghost_spec().ghost_maxrap() ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().repeat() != bkg_area_sub.ghost_spec().repeat() ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().repeat() != bkg_area_sub.ghost_spec().repeat() ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().ghost_area() != bkg_area_sub.ghost_spec().ghost_area() ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().ghost_area() != bkg_area_sub.ghost_spec().ghost_area() ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_sub.ghost_spec().grid_scatter()  ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().grid_scatter() != bkg_area_sub.ghost_spec().grid_scatter()  ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_sub.ghost_spec().pt_scatter()  ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().pt_scatter() != bkg_area_sub.ghost_spec().pt_scatter()  ||
      sub->initialJetDef().backgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_sub.ghost_spec().mean_ghost_pt()  ||
      sub->matchedJetDef().backgroundAreaDef().ghost_spec().mean_ghost_pt() != bkg_area_sub.ghost_spec().mean_ghost_pt() )
    return false;
  
  return true;
}

TEST(DijetMatrix, Default) {
    dijetcore::DijetMatrix default_matrix;
  
    default_matrix.initialize();
    std::string expected_default_string = "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2";
  
    std::set<std::string> keys = default_matrix.keys();
  
    EXPECT_NE(keys.find(expected_default_string), keys.end());
}


TEST(DijetMatrix, Initialization) {
    dijetcore::DijetMatrix matrix;
    matrix.forceConstituentPtEquality(false);
    matrix.forceConstituentEtaEquality(false);
    matrix.forceJetResolutionEquality(false);
    matrix.addJetAlgorithm({fastjet::antikt_algorithm, fastjet::kt_algorithm});
    matrix.addLeadJetR({0.4, 0.5});
    matrix.addSubJetR({0.4, 0.5});
    matrix.addLeadJetPt({20, 18});
    matrix.addSubJetPt({10, 9});
    matrix.initialize();
    
    EXPECT_EQ(matrix.size(), 32);
    
    matrix.addSubJetPt(40);
    matrix.initialize();
    
    EXPECT_EQ(matrix.size(), 32);

    matrix.addLeadJetPt({2, 60});
    matrix.initialize();
    
    EXPECT_EQ(matrix.size(), 48 + 8);

    matrix.clear();
    matrix.initialize();
    
    EXPECT_EQ(matrix.size(), 1);
}

TEST(DijetMatrix, ForceConstistuentEta) {
    dijetcore::DijetMatrix matrix;
    matrix.addConstituentEta({1.0, 1.5});
    matrix.forceConstituentEtaEquality(true);
    matrix.forceJetResolutionEquality(false);
    matrix.initialize();
    
    EXPECT_EQ(matrix.size(), 2);

    matrix.forceConstituentEtaEquality(false);
    matrix.initialize();
    
    EXPECT_EQ(matrix.size(), 4);
}

TEST(DijetMatrix, ForceConstituentPt) {
    dijetcore::DijetMatrix matrix;

    matrix.forceConstituentPtEquality(true);
    matrix.forceJetResolutionEquality(false);

    matrix.addConstituentLeadInitialPt(2.0);
    matrix.addConstituentSubInitialPt(2.1);

    matrix.initialize();
    EXPECT_EQ(matrix.size(), 0);

    matrix.forceConstituentPtEquality(false);

    matrix.initialize();
    EXPECT_EQ(matrix.size(), 1);

    matrix.clear();

    matrix.addConstituentLeadInitialPt({2.0, 2.2});
    matrix.addConstituentSubInitialPt({2.1, 2.2});
    matrix.initialize();

    EXPECT_EQ(matrix.size(), 4);

    matrix.forceConstituentPtEquality(true);
    matrix.initialize();

    EXPECT_EQ(matrix.size(), 1);

    matrix.addConstituentLeadMatchPt({0.2, 0.3});
    matrix.addConstituentSubMatchPt(0.2);
    matrix.initialize();

    EXPECT_EQ(matrix.size(), 1);

    matrix.forceConstituentPtEquality(false);
    matrix.initialize();
    
    EXPECT_EQ(matrix.size(), 8);

}

TEST(DijetMatrix, SingleCaseFullTest) {
    dijetcore::DijetMatrix matrix;
    matrix.forceConstituentPtEquality(true);
    matrix.forceConstituentEtaEquality(true);
    matrix.forceJetResolutionEquality(false);
    matrix.forceMatchJetResolutionEquality(true);
    matrix.addLeadJetPt(18);
    matrix.addLeadJetR(0.4);
    matrix.addLeadJetR(0.5);
    matrix.addSubJetR(0.4);
    matrix.addSubJetR(0.6);
    matrix.addSubJetPt(9);
    matrix.addConstituentEta(1.0);
    matrix.addConstituentLeadInitialPt(2.5);
    matrix.addConstituentLeadMatchPt(0.5);
    matrix.addConstituentSubInitialPt(2.3);
    matrix.addConstituentSubMatchPt(0.4);
    matrix.initialize();

    EXPECT_EQ(matrix.size(), 0);

    matrix.forceConstituentPtEquality(false);
    matrix.forceConstituentEtaEquality(false);
    matrix.forceMatchJetResolutionEquality(false);
    matrix.initialize();

    EXPECT_EQ(matrix.size(), 4);
  
  std::vector<std::pair<double, double>> lead_jet_pairs{{0.4, 0.4}, {0.5, 0.4}};
  std::vector<std::pair<double, double>> sub_jet_pairs{{0.4, 0.4}, {0.6, 0.4}};
  
  for (auto key : matrix.keys()) {
    dijetcore::DijetKey parsed_key = dijetcore::ParseStringToDijetKey(key);
    
    if (std::find(lead_jet_pairs.begin(), lead_jet_pairs.end(), std::pair<double, double>{parsed_key.lead_init_r, parsed_key.lead_match_r})
        == lead_jet_pairs.end()) {
      EXPECT_STREQ("Found bad DijetDefinition in DijetMatrix", "");
      continue;
    }
    double lead_R = parsed_key.lead_init_r;
    double lead_match_R = parsed_key.lead_match_r;
    
    if (std::find(sub_jet_pairs.begin(), sub_jet_pairs.end(), std::pair<double, double>{parsed_key.sub_init_r, parsed_key.sub_match_r})
        == sub_jet_pairs.end()) {
      EXPECT_STREQ("Found bad DijetDefinition in DijetMatrix", "");
      continue;
    }
    double sub_R = parsed_key.sub_init_r;
    double sub_match_R = parsed_key.sub_match_r;
    // build the background area def
    
    fastjet::GhostedAreaSpec ghost_spec_lead(1.0 + 0.4, 1, 0.01, 1, 0.1, 1e-100);
    fastjet::AreaDefinition bkg_area_def_lead(fastjet::active_area_explicit_ghosts, ghost_spec_lead);
    
    fastjet::GhostedAreaSpec ghost_spec_sub(1.0 + 0.4, 1, 0.01, 1, 0.1, 1e-100);
    fastjet::AreaDefinition bkg_area_def_sub(fastjet::active_area_explicit_ghosts, ghost_spec_sub);
    
    EXPECT_TRUE(CheckDijetDefinition(*matrix.dijetDefinitions()[key], fastjet::antikt_algorithm, fastjet::antikt_algorithm,
                                     lead_R, lead_match_R, sub_R, sub_match_R, 18, 9, 2.5, 0.5, 2.3, 0.4, 1.0, fastjet::E_scheme,
                                     fastjet::Best, fastjet::active_area_explicit_ghosts, 1, 0.01, 1, 0.1, 1e-100,
                                     fastjet::JetDefinition(fastjet::kt_algorithm, 0.4), bkg_area_def_lead,
                                     bkg_area_def_sub));
    
  }
}

TEST(DijetMatrix, TestClusterSequenceSets) {
    dijetcore::DijetMatrix matrix;
    matrix.forceJetResolutionEquality(false);
    matrix.addLeadJetPt({20.0, 18.0});
    matrix.addSubJetPt({10.0, 9.0});
    matrix.addLeadJetR(0.4);
    matrix.addSubJetR(0.4);
    matrix.addConstituentEta(0.4);
    matrix.addConstituentLeadInitialPt(2.0);
    matrix.addConstituentSubInitialPt(2.0);
    matrix.addConstituentLeadMatchPt(0.2);
    matrix.addConstituentSubMatchPt(0.2);
    matrix.initialize();
    
    ASSERT_EQ(matrix.sortedDefinitions().size(), 1);
    EXPECT_EQ(matrix.sortedDefinitions()[0].size(), 4);

    matrix.addLeadJetR(0.5);
    matrix.addSubJetR(0.5);
    matrix.initialize();

    ASSERT_EQ(matrix.sortedDefinitions().size(), 4);

    for (int i = 0; i < matrix.sortedDefinitions().size(); ++i) {
        EXPECT_EQ(matrix.sortedDefinitions()[i].size(), 4);
        
        auto& defs = matrix.sortedDefinitions()[i];

        std::set<double> lead_init_r;
        std::set<double> lead_match_r;
        std::set<double> sub_init_r;
        std::set<double> sub_match_r;

        for(auto& e : defs) {
            lead_init_r.insert(e.second->lead->initialJetDef().R());
            lead_match_r.insert(e.second->lead->matchedJetDef().R());
            sub_init_r.insert(e.second->sub->initialJetDef().R());
            sub_match_r.insert(e.second->sub->matchedJetDef().R());
        }

        EXPECT_EQ(lead_init_r.size(), 1);
        EXPECT_EQ(lead_match_r.size(), 1);
        EXPECT_EQ(sub_init_r.size(), 1);
        EXPECT_EQ(sub_match_r.size(), 1);
    }

    matrix.addConstituentLeadInitialPt({1.0, 3.0});
    matrix.addConstituentLeadMatchPt({0.3, 0.4});
    matrix.addConstituentSubInitialPt({1.0, 3.0});
    matrix.addConstituentSubMatchPt({0.3, 0.4});
    matrix.initialize();

    for (int i = 0; i < matrix.sortedDefinitions().size(); ++i) {
        EXPECT_EQ(matrix.sortedDefinitions()[i].size(), 4);
   
        auto& defs = matrix.sortedDefinitions()[i];

        std::set<double> lead_init_r;
        std::set<double> lead_match_r;
        std::set<double> sub_init_r;
        std::set<double> sub_match_r;
        std::set<double> lead_init_const_pt;
        std::set<double> lead_match_const_pt;
        std::set<double> sub_init_const_pt;
        std::set<double> sub_match_const_pt;

        for(auto& e : defs) {
            lead_init_r.insert(e.second->lead->initialJetDef().R());
            lead_match_r.insert(e.second->lead->matchedJetDef().R());
            sub_init_r.insert(e.second->sub->initialJetDef().R());
            sub_match_r.insert(e.second->sub->matchedJetDef().R());
            lead_init_const_pt.insert(dijetcore::ExtractDoubleFromSelector(e.second->lead->initialJetDef().constituentSelector(), "pt >="));
            lead_match_const_pt.insert(dijetcore::ExtractDoubleFromSelector(e.second->lead->matchedJetDef().constituentSelector(), "pt >="));
            sub_init_const_pt.insert(dijetcore::ExtractDoubleFromSelector(e.second->sub->initialJetDef().constituentSelector(), "pt >="));
            sub_match_const_pt.insert(dijetcore::ExtractDoubleFromSelector(e.second->sub->matchedJetDef().constituentSelector(), "pt >="));
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


