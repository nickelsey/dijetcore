#ifndef DIJETCORE_UTIL_DATA_EFFICIENCY_RUN6_EFF_H
#define DIJETCORE_UTIL_DATA_EFFICIENCY_RUN6_EFF_H

#include <vector>
#include <string>

#include "TF2.h"
#include "TF1.h"
#include "TFile.h"

namespace dijetcore {
  
  class Run6Eff {
    public:
    
    Run6Eff();
    
    ~Run6Eff();
    
    // get the efficiency of the TPC at a given track
    // (pT, eta)
    double ppEff(double pt, double eta);

    // returns a scaling factor for pT 
    // for instance, if smearPt() returns 1.3, then
    // smearedPt = 1.3 * pT
    double smearPt(double pt);

    TF2* getEffCurve() { return effY06; }
    TF1* getResCurve() { return resY06; }

    private:

    void LoadCurves();
    
    TF2* GetEffY06();
    TF1* GetResY06();
    
    TF2* effY06; // Run 6 parameterization

    TF1* resY06; // Run 6 pT resolution

    TF1* gaus; // Used for pT smearing
    
    double maxPtpp;
  };
  
} // namespace dijetcore
  
#endif // DIJETCORE_UTIL_DATA_EFFICIENCY_RUN6_EFF_H
