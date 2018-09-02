// pt_distribution_test.cc

#include "jet_analysis/util/pt_distribution.hh"

#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"

#include <iostream>
#include <chrono>

int main() {
  
  // create the RNG
  // standard mersenne_twister_engine seeded with a constant,
  // so that it is reproducible
  std::mt19937 gen(14342);
  
  // my pt distribution
  pt_distribution<double> dis;
  
  TH1D* h = new TH1D("tmp", "h", 100, 0, 5);
  TH1D* h2 = new TH1D("tmp1", "h2", 100, 0, 5);
  TF1* f = new TF1("f", "[0]*x*exp(-x/[1])", 0, 5);
  f->SetParameters(1, 0.291);
  
  for (int i = 0; i < 1000000; ++i) {
    h->Fill(dis(gen));
  }
  
  h->Scale(1.0/h->Integral());
  h->Fit("f");

  if (f->GetChisquare()/f->GetNDF() > 2.0)
    return 1;
    
  return 0;
}
