#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_matrix.h"

#include "dijetcore/util/fastjet/selector_compare.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/lib/logging.h"

namespace dijetcore {
  
  // default constructor
  DijetMatrix::DijetMatrix() :
  force_constituent_pt_equality_(true),
  force_constituent_eta_equality_(true),
  force_jet_resolution_equality_(true),
  force_match_jet_resolution_equality_(false),
  strategy_(fastjet::Best),
  scheme_(fastjet::E_scheme),
  area_type_(fastjet::active_area_explicit_ghosts),
  ghost_repeat_(fastjet::gas::def_repeat),
  ghost_area_(fastjet::gas::def_ghost_area),
  grid_scatter_(fastjet::gas::def_grid_scatter),
  pt_scatter_(fastjet::gas::def_pt_scatter),
  mean_ghost_pt_(fastjet::gas::def_mean_ghost_pt),
  bkg_definition_(fastjet::JetDefinition(fastjet::kt_algorithm, 0.4)) { }
  
  // single entry construction
  DijetMatrix::DijetMatrix(fastjet::JetAlgorithm jet_alg_in,
                           double lead_pt_in,
                           double lead_R_in,
                           double lead_R_match_in,
                           double sub_pt_in,
                           double sub_R_in,
                           double sub_R_match_in,
                           double const_lead_pt_init_in,
                           double const_lead_pt_match_in,
                           double const_sub_pt_init_in,
                           double const_sub_pt_match_in,
                           double eta_in) :
  const_eta_(std::set<double>{eta_in}),
  const_lead_pt_init_(std::set<double>{const_lead_pt_init_in}),
  const_lead_pt_match_(std::set<double>{const_lead_pt_match_in}),
  const_sub_pt_init_(std::set<double>{const_sub_pt_init_in}),
  const_sub_pt_match_(std::set<double>{const_sub_pt_match_in}),
  jet_algorithm_(std::set<fastjet::JetAlgorithm>{jet_alg_in}),
  lead_pt_(std::set<double>{lead_pt_in}),
  lead_R_(std::set<double>{lead_R_in}),
  lead_R_match_(std::set<double>{lead_R_match_in}),
  sub_pt_(std::set<double>{sub_pt_in}),
  sub_R_(std::set<double>{sub_R_in}),
  sub_R_match_(std::set<double>{sub_R_match_in}),
  force_constituent_pt_equality_(true),
  force_constituent_eta_equality_(true),
  force_jet_resolution_equality_(true),
  force_match_jet_resolution_equality_(false),
  strategy_(fastjet::Best),
  scheme_(fastjet::E_scheme),
  area_type_(fastjet::active_area_explicit_ghosts),
  ghost_repeat_(fastjet::gas::def_repeat),
  ghost_area_(fastjet::gas::def_ghost_area),
  grid_scatter_(fastjet::gas::def_grid_scatter),
  pt_scatter_(fastjet::gas::def_pt_scatter),
  mean_ghost_pt_(fastjet::gas::def_mean_ghost_pt),
  bkg_definition_(fastjet::JetDefinition(fastjet::kt_algorithm, 0.4))
  { }
  
  
  // set construction
  DijetMatrix::DijetMatrix(std::set<fastjet::JetAlgorithm> jet_alg_in,
                           std::set<double> lead_pt_in,
                           std::set<double> lead_R_in,
                           std::set<double> lead_R_match_in,
                           std::set<double> sub_pt_in,
                           std::set<double> sub_R_in,
                           std::set<double> sub_R_match_in,
                           std::set<double> const_lead_pt_init_in,
                           std::set<double> const_lead_pt_match_in,
                           std::set<double> const_sub_pt_init_in,
                           std::set<double> const_sub_pt_match_in,
                           std::set<double> eta_in) :
  const_eta_(eta_in),
  const_lead_pt_init_(const_lead_pt_init_in),
  const_lead_pt_match_(const_lead_pt_match_in),
  const_sub_pt_init_(const_sub_pt_init_in),
  const_sub_pt_match_(const_sub_pt_match_in),
  jet_algorithm_(jet_alg_in),
  lead_pt_(lead_pt_in),
  lead_R_(lead_R_in),
  lead_R_match_(lead_R_match_in),
  sub_pt_(sub_pt_in),
  sub_R_(sub_R_in),
  sub_R_match_(sub_R_match_in),
  force_constituent_pt_equality_(true),
  force_constituent_eta_equality_(true),
  force_jet_resolution_equality_(true),
  force_match_jet_resolution_equality_(false),
  strategy_(fastjet::Best),
  scheme_(fastjet::E_scheme),
  area_type_(fastjet::active_area_explicit_ghosts),
  ghost_repeat_(fastjet::gas::def_repeat),
  ghost_area_(fastjet::gas::def_ghost_area),
  grid_scatter_(fastjet::gas::def_grid_scatter),
  pt_scatter_(fastjet::gas::def_pt_scatter),
  mean_ghost_pt_(fastjet::gas::def_mean_ghost_pt),
  bkg_definition_(fastjet::JetDefinition(fastjet::kt_algorithm, 0.4))
  { }
  
  DijetMatrix::DijetMatrix(const DijetMatrix& rhs) :
  const_eta_(rhs.const_eta_),
  const_lead_pt_init_(rhs.const_lead_pt_init_),
  const_lead_pt_match_(rhs.const_lead_pt_match_),
  const_sub_pt_init_(rhs.const_sub_pt_init_),
  const_sub_pt_match_(rhs.const_sub_pt_init_),
  jet_algorithm_(rhs.jet_algorithm_),
  lead_pt_(rhs.lead_pt_),
  lead_R_(rhs.lead_R_),
  lead_R_match_(rhs.lead_R_match_),
  sub_pt_(rhs.sub_pt_),
  sub_R_(rhs.sub_R_),
  sub_R_match_(rhs.sub_R_match_),
  force_constituent_pt_equality_(rhs.force_constituent_pt_equality_),
  force_constituent_eta_equality_(rhs.force_constituent_eta_equality_),
  force_jet_resolution_equality_(rhs.force_jet_resolution_equality_),
  force_match_jet_resolution_equality_(rhs.force_match_jet_resolution_equality_),
  strategy_(rhs.strategy_),
  scheme_(rhs.scheme_),
  area_type_(rhs.area_type_),
  ghost_repeat_(rhs.ghost_repeat_),
  ghost_area_(rhs.ghost_area_),
  grid_scatter_(rhs.grid_scatter_),
  pt_scatter_(rhs.pt_scatter_),
  mean_ghost_pt_(rhs.mean_ghost_pt_),
  bkg_definition_(rhs.bkg_definition_)
  { }
  
  void DijetMatrix::Clear() {
    const_eta_.clear();
    const_lead_pt_init_.clear();
    const_sub_pt_init_.clear();
    const_lead_pt_match_.clear();
    const_sub_pt_match_.clear();
    jet_algorithm_.clear();
    lead_pt_.clear();
    lead_R_.clear();
    lead_R_match_.clear();
    sub_pt_.clear();
    sub_R_.clear();
    sub_R_match_.clear();
    
    ClearDijetDefs();
  }
  
  void DijetMatrix::ClearDijetDefs() {
    dijet_defs_.clear();
    ordered_defs_.clear();
    lead_matchdefs_.clear();
    sub_matchdefs_.clear();
    keys_.clear();
    ordered_keys_.clear();
  }
  
  void DijetMatrix::ForceConstituentPtEquality(bool flag) {
    force_constituent_pt_equality_ = flag;
    CheckToUpdate();
  }
  
  void DijetMatrix::ForceConstituentEtaEquality(bool flag) {
    force_constituent_eta_equality_ = flag;
    CheckToUpdate();
  }
  
  void DijetMatrix::ForceJetResolutionEquality(bool flag) {
    force_jet_resolution_equality_ = flag;
    CheckToUpdate();
  }
  
  void DijetMatrix::ForceMatchJetResolutionEquality(bool flag) {
    force_match_jet_resolution_equality_ = flag;
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentEta(double eta) {
    const_eta_.insert(eta);
    CheckToUpdate();
  }
  void DijetMatrix::AddConstituentEta(std::set<double> eta) {
    for (auto& val : eta) {
      const_eta_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentLeadInitialPt(double pt) {
    const_lead_pt_init_.insert(pt);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentLeadInitialPt(std::set<double> pt) {
    for (auto& val : pt) {
      const_lead_pt_init_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentLeadMatchPt(double pt) {
    const_lead_pt_match_.insert(pt);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentLeadMatchPt(std::set<double> pt) {
    for (auto& val : pt) {
      const_lead_pt_match_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentSubInitialPt(double pt) {
    const_sub_pt_init_.insert(pt);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentSubInitialPt(std::set<double> pt) {
    for (auto& val : pt) {
      const_sub_pt_init_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentSubMatchPt(double pt) {
    const_sub_pt_match_.insert(pt);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddConstituentSubMatchPt(std::set<double> pt) {
    for (auto& val : pt) {
      const_sub_pt_match_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddJetAlgorithm(fastjet::JetAlgorithm alg) {
    jet_algorithm_.insert(alg);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddJetAlgorithm(std::set<fastjet::JetAlgorithm> alg) {
    for (auto& val : alg) {
      jet_algorithm_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddLeadJetPt(double pt) {
    lead_pt_.insert(pt);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddLeadJetPt(std::set<double> pt) {
    for (auto& val : pt) {
      lead_pt_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddLeadJetR(double R) {
    lead_R_.insert(R);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddLeadJetR(std::set<double> R) {
    for (auto& val : R) {
      lead_R_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddSubJetPt(double pt) {
    sub_pt_.insert(pt);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddSubJetPt(std::set<double> pt) {
    for (auto& val : pt) {
      sub_pt_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::AddSubJetR(double R) {
    sub_R_.insert(R);
    CheckToUpdate();
  }
  
  void DijetMatrix::AddSubJetR(std::set<double> R) {
    for (auto& val : R) {
      sub_R_.insert(val);
    }
    CheckToUpdate();
  }
  
  void DijetMatrix::Initialize() {
    // if already initialzed, remove old definitions
    if (dijet_defs_.size() != 0) {
      ClearDijetDefs();
    }
    
    // check to make sure there is at least one valid dijet
    // definition, either through custom definitions,
    // or parameters. Create default parameters otherwise
    InitializeEmptyFields();
    
    // first step - get our leading & subleading jet definitions
    std::vector<fastjet::JetDefinition> lead_jet_defs = FillLeadJetDefinitions();
    std::vector<fastjet::JetDefinition> lead_match_jet_defs = FillLeadMatchJetDefinitions();
    std::vector<fastjet::JetDefinition> sub_jet_defs = FillSubJetDefinitions();
    std::vector<fastjet::JetDefinition> sub_match_jet_defs = FillSubMatchJetDefinitions();
    
    // start creating the leading jet MatchDefs
    for (auto& lead_fj_def : lead_jet_defs) {
      for(auto& lead_match_fj_def : lead_match_jet_defs) {
        
        // we always use the same algorithm - if they are different, then discard
        if (lead_fj_def.jet_algorithm() != lead_match_fj_def.jet_algorithm())
          continue;
        
        // if we are forcing matched jets to have similar R, check here
        if (force_match_jet_resolution_equality_ &&
            lead_fj_def.R() != lead_match_fj_def.R())
          continue;
        
        for (auto& lead_jet_pt : lead_pt_) {
          for (auto& lead_const_pt_init : const_lead_pt_init_) {
            for (auto& lead_const_pt_match : const_lead_pt_match_) {
              for (auto& lead_const_eta : const_eta_) {
                
                // get a few needed parameters
                double R = lead_fj_def.R();
                double bkg_R = bkg_definition_.R();
                double jet_eta_max = lead_const_eta - R;
                double bkg_jet_eta_max = lead_const_eta - bkg_R;
                
                // build constituent & jet selectors
                fastjet::Selector init_const_selector = fastjet::SelectorAbsRapMax(lead_const_eta)
                && fastjet::SelectorPtMin(lead_const_pt_init);
                fastjet::Selector match_const_selector = fastjet::SelectorAbsRapMax(lead_const_eta)
                && fastjet::SelectorPtMin(lead_const_pt_match);
                fastjet::Selector init_jet_selector = fastjet::SelectorAbsRapMax(jet_eta_max)
                && fastjet::SelectorPtMin(lead_jet_pt);
                fastjet::Selector match_jet_selector = fastjet::SelectorIdentity();
                
                fastjet::GhostedAreaSpec ghost_def(jet_eta_max + 2 * R, ghost_repeat_, ghost_area_,
                                                   grid_scatter_, pt_scatter_, mean_ghost_pt_);
                fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,
                                                 ghost_def);
                
                fastjet::Selector bkg_selector = fastjet::SelectorAbsRapMax(lead_const_eta - bkg_R)
                * (!fastjet::SelectorNHardest(2));
                fastjet::GhostedAreaSpec bkg_ghost_def(bkg_jet_eta_max + 2 * bkg_R, ghost_repeat_, ghost_area_,
                                                       grid_scatter_, pt_scatter_, mean_ghost_pt_);
                fastjet::AreaDefinition bkg_area_def(fastjet::active_area_explicit_ghosts,
                                                     bkg_ghost_def);
                // create the initial JetDef & the
                // matched JetDef def
                JetDef init_def(lead_fj_def, area_def, bkg_definition_, bkg_area_def, bkg_selector);
                init_def.SetConstituentSelector(init_const_selector);
                init_def.SetJetSelector(init_jet_selector);
                JetDef match_def(lead_match_fj_def, area_def, bkg_definition_, bkg_area_def, bkg_selector);
                match_def.SetConstituentSelector(match_const_selector);
                match_def.SetJetSelector(match_jet_selector);
                
                // create the MatchDef
                lead_matchdefs_.insert(make_unique<MatchDef>(init_def, match_def));
              }
            }
          }
        }
      }
    }
    
    // do the same for subleading
    for (auto& sub_fj_def : sub_jet_defs) {
      for(auto& sub_match_fj_def : sub_match_jet_defs) {
        
        // we always use the same algorithm - if they are different, then discard
        if (sub_fj_def.jet_algorithm() != sub_match_fj_def.jet_algorithm())
        continue;
        
        // if we are forcing matched jets to have similar R, check here
        if (force_match_jet_resolution_equality_ &&
            sub_fj_def.R() != sub_match_fj_def.R())
        continue;
        
        for (auto& sub_jet_pt : sub_pt_) {
          for (auto& sub_const_pt_init : const_sub_pt_init_) {
            for (auto& sub_const_pt_match : const_sub_pt_match_) {
              for (auto& sub_const_eta : const_eta_) {
                
                // get a few needed parameters
                double R = sub_fj_def.R();
                double bkg_R = bkg_definition_.R();
                double jet_eta_max = sub_const_eta - R;
                double bkg_jet_eta_max = sub_const_eta - bkg_R;
                
                // build constituent & jet selectors
                fastjet::Selector init_const_selector = fastjet::SelectorAbsRapMax(sub_const_eta)
                && fastjet::SelectorPtMin(sub_const_pt_init);
                fastjet::Selector match_const_selector = fastjet::SelectorAbsRapMax(sub_const_eta)
                && fastjet::SelectorPtMin(sub_const_pt_match);
                fastjet::Selector init_jet_selector = fastjet::SelectorAbsRapMax(jet_eta_max)
                && fastjet::SelectorPtMin(sub_jet_pt);
                fastjet::Selector match_jet_selector = fastjet::SelectorIdentity();
                fastjet::Selector bkg_selector = fastjet::SelectorAbsRapMax(sub_const_eta - bkg_R)
                * (!fastjet::SelectorNHardest(2));
                
                fastjet::GhostedAreaSpec ghost_def(jet_eta_max + 2 * R, ghost_repeat_, ghost_area_,
                                                   grid_scatter_, pt_scatter_, mean_ghost_pt_);
                fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,
                                                 ghost_def);
                
                fastjet::GhostedAreaSpec bkg_ghost_def(bkg_jet_eta_max + 2 * bkg_R, ghost_repeat_, ghost_area_,
                                                       grid_scatter_, pt_scatter_, mean_ghost_pt_);
                fastjet::AreaDefinition bkg_area_def(fastjet::active_area_explicit_ghosts,
                                                     bkg_ghost_def);
                
                // create the initial JetDef & the
                // matched JetDef def
                JetDef init_def(sub_fj_def, area_def, bkg_definition_, bkg_area_def, bkg_selector);
                init_def.SetConstituentSelector(init_const_selector);
                init_def.SetJetSelector(init_jet_selector);
                JetDef match_def(sub_match_fj_def, area_def, bkg_definition_, bkg_area_def, bkg_selector);
                match_def.SetConstituentSelector(match_const_selector);
                match_def.SetJetSelector(match_jet_selector);
                
                // create the MatchDef
                sub_matchdefs_.insert(make_unique<MatchDef>(init_def, match_def));
                
              }
            }
          }
        }
      }
    }
    
    // create the dijet definitions
    for (auto& lead : lead_matchdefs_) {
      for (auto& sub : sub_matchdefs_) {
        
        // make sure we're looking at similar jets
        if (lead->InitialJetDef().jet_algorithm() != sub->InitialJetDef().jet_algorithm())
          continue;
        
        // make sure that the pt lead > pt sub
        if (SelectorPtMinLessThan(lead->InitialJetDef().JetSelector(), sub->InitialJetDef().JetSelector()))
          continue;
        
        // if force_constituent_pt_equality is on, force
        // pt to be equal. Same for eta
        if (force_constituent_pt_equality_) {
          double lead_pt_init = ExtractDoubleFromSelector(lead->InitialJetDef().ConstituentSelector(), "pt >=");
          double sub_pt_init = ExtractDoubleFromSelector(sub->InitialJetDef().ConstituentSelector(), "pt >=");
          if (lead_pt_init != sub_pt_init)
            continue;
          double lead_pt_match = ExtractDoubleFromSelector(lead->MatchedJetDef().ConstituentSelector(), "pt >=");
          double sub_pt_match = ExtractDoubleFromSelector(sub->MatchedJetDef().ConstituentSelector(), "pt >=");
          if (lead_pt_match != sub_pt_match)
            continue;
        }
        if (force_constituent_eta_equality_) {
          double lead_eta = ExtractDoubleFromSelector(lead->InitialJetDef().ConstituentSelector(), "|rap| <=");
          double sub_eta = ExtractDoubleFromSelector(sub->InitialJetDef().ConstituentSelector(), "|rap| <=");
          if (lead_eta != sub_eta)
            continue;
        }
        
        // if the force_jet_resolution_equality_ is on, force
        // R to be equal
        if (force_jet_resolution_equality_) {
          if (lead->InitialJetDef().R() != sub->InitialJetDef().R() ||
              lead->MatchedJetDef().R() != sub->MatchedJetDef().R())
            continue;
        }
        
        auto tmp = make_unique<DijetDefinition>(lead.get(), sub.get(),
                                                     std::min(lead->InitialJetDef().R(),
                                                              sub->InitialJetDef().R()));
        
        string key = MakeKeyFromDijetDefinition(*tmp);
        dijet_defs_[key] = std::move(tmp);
        keys_.insert(key);
      }
    }
    
    // finally, make sets of equivalent dijet definitions
    // loop over each dijet_definition
    for (auto& new_dijet : dijet_defs_) {
      const string& new_key = new_dijet.first;
      DijetDefinition* new_def = new_dijet.second.get();
      
      bool match_success = false;
      
      // look through each sorted set of dijet definitions
      // if it matched
      for (int i = 0; i < ordered_defs_.size(); ++i) {
        const string& set_key = ordered_defs_[i].begin()->first;
        DijetDefinition* set_def = ordered_defs_[i].begin()->second;
      
        if (set_def->EquivalentCluster(*new_def)) {
      
          ordered_defs_[i][new_key] = new_def;
          match_success = true;
      
          double lead_pt_min = ExtractDoubleFromSelector(new_def->lead->InitialJetDef().JetSelector(), "pt >=");
          double sub_pt_min = ExtractDoubleFromSelector(new_def->sub->InitialJetDef().JetSelector(), "pt >=");
      
          if (lead_pt_min < ordered_min_pt_[i].first)
            ordered_min_pt_[i].first = lead_pt_min;
          if (sub_pt_min < ordered_min_pt_[i].second)
            ordered_min_pt_[i].second = sub_pt_min;
  
          break;
        }
      }
      
      if (!match_success) {
        std::unordered_map<string, DijetDefinition*> tmp;
        tmp.insert({new_key, new_def});
        ordered_defs_.push_back(std::move(tmp));
        ordered_keys_.push_back(MakeSortedKeyFromDijetDefinition(*new_def));
        double lead_pt_min = ExtractDoubleFromSelector(new_def->lead->InitialJetDef().JetSelector(), "pt >=");
        double sub_pt_min = ExtractDoubleFromSelector(new_def->sub->InitialJetDef().JetSelector(), "pt >=");
        ordered_min_pt_.push_back({lead_pt_min, sub_pt_min});
      }
    }
    
  }
  
  void DijetMatrix::CheckToUpdate() {
    if (dijet_defs_.size() != 0) {
      ClearDijetDefs();
    }
  }
  
  void DijetMatrix::InitializeEmptyFields() {
    if (const_eta_.empty())
      const_eta_.insert(1.0);
    if (const_lead_pt_init_.empty())
      const_lead_pt_init_.insert(2.0);
    if (const_lead_pt_match_.empty())
      const_lead_pt_match_.insert(0.2);
    if (const_sub_pt_init_.empty())
      const_sub_pt_init_.insert(2.0);
    if (const_sub_pt_match_.empty())
      const_sub_pt_match_.insert(0.2);
    if (jet_algorithm_.empty())
      jet_algorithm_.insert(fastjet::antikt_algorithm);
    if (lead_pt_.empty())
      lead_pt_.insert(20.0);
    if (lead_R_.empty())
      lead_R_.insert(0.4);
    if (lead_R_match_.empty())
      lead_R_match_.insert(0.4);
    if (sub_pt_.empty())
      sub_pt_.insert(10.0);
    if (sub_R_.empty())
      sub_R_.insert(0.4);
    if (sub_R_match_.empty())
      sub_R_match_.insert(0.4);
  }
  
  std::vector<fastjet::JetDefinition> DijetMatrix::FillLeadJetDefinitions() {
    std::vector<fastjet::JetDefinition> ret;
    for (auto& alg : jet_algorithm_)
      for (auto& R : lead_R_)
        ret.push_back(fastjet::JetDefinition(alg, R, scheme_, strategy_));
    return ret;
  }
  
  std::vector<fastjet::JetDefinition> DijetMatrix::FillLeadMatchJetDefinitions() {
    std::vector<fastjet::JetDefinition> ret;
    for (auto& alg : jet_algorithm_)
      for (auto& R : lead_R_match_)
        ret.push_back(fastjet::JetDefinition(alg, R, scheme_, strategy_));
    return ret;
  }
  
  std::vector<fastjet::JetDefinition> DijetMatrix::FillSubJetDefinitions() {
    std::vector<fastjet::JetDefinition> ret;
    for (auto& alg : jet_algorithm_)
      for (auto& R : sub_R_)
        ret.push_back(fastjet::JetDefinition(alg, R, scheme_, strategy_));
    return ret;
  }
  
  std::vector<fastjet::JetDefinition> DijetMatrix::FillSubMatchJetDefinitions() {
    std::vector<fastjet::JetDefinition> ret;
    for (auto& alg : jet_algorithm_)
      for (auto& R : sub_R_match_)
        ret.push_back(fastjet::JetDefinition(alg, R, scheme_, strategy_));
    return ret;
  }
  
} // namespace dijetcore
