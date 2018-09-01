#ifndef DIJETCORE_UTIL_DATA_EFFICIENCY_RUN7_EFF_H
#define DIJETCORE_UTIL_DATA_EFFICIENCY_RUN7_EFF_H

#include <vector>
#include <string>

#include "TH2.h"
#include "TF2.h"
#include "TFile.h"

namespace dijetcore {
  
  enum class TrackingUncY7 {
    NONE = 0,
    POSITIVE = 1,
    NEGATIVE = -1
  };
  
  class Run7Eff {
    public:
    
    Run7Eff(std::string filename = "resources/efficiencies/y7_effic.root");
    
    ~Run7Eff();
    
    void loadFile(std::string filename);
    
    double AuAuEff(double pt, double eta, int cent);
    double AuAuEff020Avg(double pt, double eta);
    double ppEff(double pt, double eta);
    double ratio(double pt, double eta);
    double ratioUncertainty(double pt, double eta);
    
    void setSystematicUncertainty(TrackingUncY7 sys = TrackingUncY7::NONE) {sys_ = sys;}
    
    
    private:
    
    void loadCurves();
    TF2* GetEffY06();
    TF2* GetEffY04(int cb);
    
    TFile* file;
    
    //centrality bins: 0 is 0-5%, 1 is 5-10%, 2 is 10-20%
    TF2* effY04[3]; // Run 4 parameterization
    TH2D* effY07pteta[3]; // Run 7 HT / Run 4 MB pt-eta map, used for pt <= 1.5 GeV/c
    TH1D* effY07eta[3]; // Run 7 HT / Run 4 MB eta map, used for pt > 1.5 GeV/c
    TF2* effY06; // Run 6 parameterization
    
    double maxPt;
    double maxPtpp;
    TrackingUncY7 sys_;
  };
  
} // namespace dijetcore
  
#endif // DIJETCORE_UTIL_DATA_EFFICIENCY_RUN7_EFF_H
