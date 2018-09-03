#ifndef DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_DIJET_MATRIX_H
#define DIJETCORE_WORKER_DIJET_WORKER_DIJET_MATRIX_DIJET_MATRIX_H

// the DijetMatrix defines a set of parameters
// that defines a dijet definition (as used by
// the dijet imbalance measurement). The user
// can specify a set of values for any/each parameter,
// and the DijetMatrix will create a DijetDefinition
// for each valid set of the parameters.

// .... obviously, having multiple values in, say, 10
// parameters is a bad idea. (even if only two values are
// specified for each parameter, the total number of
// DijetDefinitions would be quite large for clustering).
// But the matrix will attempt to find sets of dijet_definitions
// which can use the clustersequence & input particles, so
// that repetition can be reduced

#include <set>
#include <unordered_map>
#include <vector>

#include "dijetcore/lib/memory.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/Selector.hh"

namespace dijetcore {
  
  class DijetMatrix {
    
  public:
    
    // default constructor populates all fields that aren't filled
    // via InitializeEmptyFields()
    DijetMatrix();
    
    // single entry construction
    DijetMatrix(fastjet::JetAlgorithm jet_alg_in,
                double lead_jet_pt_in = 20.0,
                double lead_jet_R_in = 0.4,
                double lead_jet_R_match_in = 0.4,
                double sub_jet_pt_in = 10.0,
                double sub_jet_R_in = 0.4,
                double sub_jet_R_match_in = 0.4,
                double const_lead_pt_init_in = 2.0,
                double const_lead_pt_match_in = 0.2,
                double const_sub_pt_init_in = 2.0,
                double const_sub_pt_match_in = 0.2,
                double eta_in = 1.0);
    
    
    // set construction
    DijetMatrix(std::set<fastjet::JetAlgorithm> jet_alg_in,
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
                std::set<double> eta_in);
    
    DijetMatrix(const DijetMatrix& rhs);
    
    virtual ~DijetMatrix() { };
    
    // removes all parameters
    void Clear();
    // removes only the DijetDefinitions, if Initialize()
    // has been called
    void ClearDijetDefs();
    
    // used to initialize the matrix.
    void Initialize();
    
    // get the number of dijet definitions currently stored
    // (will be zero if no initialization has occured)
    inline std::size_t Size() {return dijet_defs_.size();}
    
    // get access to the dijet definitions
    std::unordered_map<string, unique_ptr<DijetDefinition>>& DijetDefinitions()
        {return dijet_defs_;}
    
    // dijet definition sets that have been selected such that
    // each set can use a single cluster sequence - only jet
    // selectors can change
    std::vector<std::unordered_map<string, DijetDefinition*>>& SortedDefinitions()
        {return ordered_defs_;}
    
    
    // force leading and subleading jets to have equal constituent
    // pt cuts or eta cuts (default true)
    void ForceConstituentPtEquality(bool flag = true);
    void ForceConstituentEtaEquality(bool flag = true);
    
    // force jets to have same R for leading/subleading (default true)
    void ForceJetResolutionEquality(bool flag = true);
    
    // force matched jets to have same R (default false)
    void ForceMatchJetResolutionEquality(bool flag = true);
    
    // get the keys to the map
    const std::set<string>& Keys() const {return keys_;}
    
    // keys sorted to match SortedDefinitions()
    const std::vector<string>& SortedKeys() const {return ordered_keys_;}
    
    // minimum pT needed for jet sets - returns {lead_hard_pt_min, sub_hard_pt_min}
    const std::vector<std::pair<double, double>>& SortedDefinitionsMinPt() {return ordered_min_pt_;}
    
    // The Matrix will build all logically consistent di-jet pairs
    // (i.e. if pt sub > pt lead, it will be ignored, that sort of
    // thing). Adding extra parameters will increase the effective
    // run time considerably, especially as the set of non-unique
    // parameters increases.
    void AddConstituentEta(double eta);
    void AddConstituentEta(std::set<double> eta);
    
    void AddConstituentLeadInitialPt(double pt);
    void AddConstituentLeadInitialPt(std::set<double> pt);
    
    void AddConstituentLeadMatchPt(double pt);
    void AddConstituentLeadMatchPt(std::set<double> pt);
    
    void AddConstituentSubInitialPt(double pt);
    void AddConstituentSubInitialPt(std::set<double> pt);
    
    void AddConstituentSubMatchPt(double pt);
    void AddConstituentSubMatchPt(std::set<double> pt);
    
    void AddJetAlgorithm(fastjet::JetAlgorithm alg);
    void AddJetAlgorithm(std::set<fastjet::JetAlgorithm> alg);
    
    void AddLeadJetPt(double pt);
    void AddLeadJetPt(std::set<double> pt);
    
    void AddLeadJetR(double R);
    void AddLeadJetR(std::set<double> R);
    
    void AddSubJetPt(double pt);
    void AddSubJetPt(std::set<double> pt);
    
    void AddSubJetR(double R);
    void AddSubJetR(std::set<double> R);
    
    // options to change the default fastjet settings
    // for area/bkg estimation
    inline void SetClusterStrategy(fastjet::Strategy strat)  {strategy_ = strat;}
    inline void SetRecombinationScheme(fastjet::RecombinationScheme schm)
        {scheme_ = schm;}
    inline void SetAreaType(fastjet::AreaType type)          {area_type_ = type;}
    inline void SetGhostRepeat(int repeat)                   {ghost_repeat_ = repeat;}
    inline void SetGhostArea(double area)                    {ghost_area_ = area;}
    inline void SetGridScatter(double scatter)               {grid_scatter_ = scatter;}
    inline void SetPtScatter(double scatter)                 {pt_scatter_ = scatter;}
    inline void SetMeanGhostPt(double mean)                  {mean_ghost_pt_ = mean;}
    inline void SetBkgDefinition(fastjet::JetDefinition def) {bkg_definition_ = def;}
    
    // access to the internally stored sets
    inline const std::set<double>& ConstituentEta() const              {return const_eta_;}
    inline const std::set<double>& LeadConstituentInitPt() const       {return const_lead_pt_init_;}
    inline const std::set<double>& LeadConstituentMatchPt() const      {return const_lead_pt_match_;}
    inline const std::set<double>& SubConstituentInitPt() const        {return const_sub_pt_init_;}
    inline const std::set<double>& SubConstituentMatchPt() const       {return const_sub_pt_match_;}
    inline const std::set<fastjet::JetAlgorithm>& JetAlgorithm() const {return jet_algorithm_;}
    inline const std::set<double>& LeadJetPt() const                   {return lead_pt_;}
    inline const std::set<double>& LeadJetR() const                    {return lead_R_;}
    inline const std::set<double>& SubJetPt() const                    {return sub_pt_;}
    inline const std::set<double>& SubJetR() const                     {return sub_R_;}
    
    inline fastjet::Strategy ClusterStrategy() const                {return strategy_;}
    inline fastjet::RecombinationScheme RecombinationScheme() const {return scheme_;}
    inline fastjet::AreaType AreaType() const                       {return area_type_;}
    
    
  protected:
    
    // used internally to update the dijet definitions, if there has
    // been a parameter updated after initialization
    void CheckToUpdate();
    
    // used internally to make sure there can be at least one valid
    // dijet definition when Initialize is called. Otherwise, it sets
    // default values for missing fields
    void InitializeEmptyFields();
    
    // used during initialization
    std::vector<fastjet::JetDefinition> FillLeadJetDefinitions();
    std::vector<fastjet::JetDefinition> FillSubJetDefinitions();
    std::vector<fastjet::JetDefinition> FillLeadMatchJetDefinitions();
    std::vector<fastjet::JetDefinition> FillSubMatchJetDefinitions();

    
    std::set<double> const_eta_;
    std::set<double> const_lead_pt_init_;
    std::set<double> const_lead_pt_match_;
    std::set<double> const_sub_pt_init_;
    std::set<double> const_sub_pt_match_;
    
    // forces constituent pt or eta cut to be equal in leading & subleading
    // jet definitions
    bool force_constituent_pt_equality_;
    bool force_constituent_eta_equality_;
    
    // forces jet R to be similar either between leading/subleading or between
    // initial/matched
    bool force_jet_resolution_equality_;
    bool force_match_jet_resolution_equality_;
    
    // allow different radii for leading & subleading
    // jets, but do not allow different jet algorithms
    std::set<fastjet::JetAlgorithm> jet_algorithm_;
    std::set<double> lead_pt_;
    std::set<double> lead_R_;
    std::set<double> lead_R_match_;
    std::set<double> sub_pt_;
    std::set<double> sub_R_;
    std::set<double> sub_R_match_;
    
    // details of fastjet that one might want to change from default,
    // for one reason or another
    
    // recombination strategy, probably shouldn't change. Default is set
    // to 'best' which tells fastjet to choose the most efficient
    fastjet::Strategy strategy_;
    
    // recombination scheme
    fastjet::RecombinationScheme scheme_;
    
    // the type of area estimation done, default is active area with
    // explicit ghosts
    fastjet::AreaType area_type_;
    
    // ghost repeating, default = 1
    int ghost_repeat_;
    
    // ghost area (smaller = more ghosts per unit area)
    double ghost_area_;
    
    // for completeness, these are included - can look up in fastjet
    // manual for details.
    double grid_scatter_;
    double pt_scatter_;
    double mean_ghost_pt_;
    
    // background jet definition
    fastjet::JetDefinition bkg_definition_;
    
    // populated by the DijetMatrix
    std::unordered_map<string, unique_ptr<DijetDefinition>> dijet_defs_;
    std::set<unique_ptr<MatchDef>> lead_matchdefs_;
    std::set<unique_ptr<MatchDef>> sub_matchdefs_;
    
    // set of dijet_definitions that can be run with a single
    // cluster sequence for more efficiency
    std::vector<std::unordered_map<string, DijetDefinition*>> ordered_defs_;
    
    // keys for the dictionary of DijetDefintions are identifier
    // strings describing the detailed settings of each definition.
    // the keys are the minimum length needed to be unique to each
    // dijet definition, not including background subtraction settings
    std::set<string> keys_;
    
    // organized keys
    std::vector<string> ordered_keys_;
    std::vector<std::pair<double, double>> ordered_min_pt_;
  };
  
} // namespace dijetcore

#endif // DIJET_MATRIX_HH
