#ifndef DIJETCORE_UTIL_DATA_EFFICIENCY_RUN4_EFF_H
#define DIJETCORE_UTIL_DATA_EFFICIENCY_RUN4_EFF_H

#include <vector>
#include <string>

#include "TF2.h"

namespace dijetcore {

class Run4Eff {
public:
  Run4Eff();
  
  ~Run4Eff();
  
  double AuAuEff(double pt, double eta, int cent);
  double PPEff(double pt, double eta);
  double AuAuPPRatio(double pt, double eta, int cent=0);
  double CentRatio(double pt, double eta, int cent, int reference_cent=0);
  
private:
  std::vector<double> parset0;
  std::vector<double> parset1;
  std::vector<double> parset2;
  std::vector<double> parset3;
  std::vector<double> parset4;
  std::vector<double> parset5;
  std::vector<double> parset6;
  std::vector<double> parsetpp;
  
  double maxPt;
  double maxPtpp;
  
  TF2* y4_functions[9];
  TF2* y6_function;
  
};

} // namespace dijetcore

#endif // DIJETCORE_UTIL_DATA_EFFICIENCY_RUN4_EFF_H
