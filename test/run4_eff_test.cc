// run4_eff_test.cc

#include "jet_analysis/efficiency/run4_eff.hh"

#include <random>

int main() {
  
  // testing to make sure all efficiencies are
  // logically valid (i.e, >= 0, <= 1)
  
  Run4Eff efficiency;
  
  // define distributions for pt, eta, phi
  double eta_max = 1.0;
  double eta_min = -eta_max;
  std::uniform_real_distribution<> eta_dist(eta_min, eta_max);
  
  double pt_auau_max = 5.0;
  double pt_auau_min = 0.2;
  std::uniform_real_distribution<> pt_auau_dist(pt_auau_min, pt_auau_max);
  
  double pt_pp_max = 10.0;
  double pt_pp_min = 0.2;
  std::uniform_real_distribution<> pt_pp_dist(pt_pp_min, pt_pp_max);
  
  int cent_max = 8;
  double cent_min = 0;
  std::uniform_int_distribution<> cent_dist(cent_min, cent_max);
  // create the RNG
  // standard mersenne_twister_engine seeded with a constant,
  // so that it is reproducible
  std::mt19937 gen(14342);
  
  for (int i = 0; i < 1000; ++i) {
    double eta = eta_dist(gen);
    double pt_auau = pt_auau_dist(gen);
    double pt_pp = pt_pp_dist(gen);
    int cent = cent_dist(gen);

    if (efficiency.AuAuEff(pt_auau, eta, cent) > 1.0 ||
        efficiency.AuAuEff(pt_auau, eta, cent) < 0.0)
      return 1;
    if (efficiency.PPEff(pt_pp, eta) > 1.0 ||
        efficiency.PPEff(pt_pp, eta) < 0.0)
      return 1;
    
    
  }
  
  return 0;
}
