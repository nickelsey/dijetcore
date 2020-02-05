#include "dijetcore/util/data/efficiency/run6_eff.h"

namespace dijetcore {
  
  Run6Eff::Run6Eff() : maxPtpp(9.0) {
    LoadCurves();
  }
  
  Run6Eff::~Run6Eff() {
    delete effY06;
  }
  
  double Run6Eff::ppEff(double pt, double eta) {
    return  effY06->Eval(pt, eta);
  }

  double Run6Eff::smearPt(double pt) {
    double width = resY06->Eval(pt);
    gaus->SetParameter(2, width);
    double smear_val = gaus->GetRandom();
    return (1.0 + smear_val);
  }
  
  void Run6Eff::LoadCurves() {
    effY06 = GetEffY06();

    resY06 = GetResY06();

    gaus = new TF1("tmpgausrun6eff", "gaus(0)", -5, 5);
    gaus->SetParameters(1.0, 0.0, 1.0);
  }
  
  TF2* Run6Eff::GetEffY06()
  {
    
    TF2* funcpp = new TF2("effY06","[0]-0.06-[1]*exp([2]/x)+[3]*exp(-0.5*((x-[4])/[5])**2)/sqrt(2*pi*[5]*[5])-[6]*exp(-0.5*((x-[7])/[8])**2)/sqrt(2*pi*[8]*[8])+([9]-[10]*(y-[11])^2-[12]*(y-[11])^4-[13]*(y-[11])^6-[14]*(y-[11])^8)*exp(-[15]*x)",0.,10.,-1.,1.);
    
    Double_t parset[] = {0.869233,0.0223402,0.44061,0.558762,0.145162,0.508033,110.008,-4.63659,1.73765,0.0452674,-0.101279,0.0081551,0.945287,-2.00949,1.61746,1.39352};
    
    ((TF2*)funcpp)->SetParameters(parset);
    
    return funcpp;
  }

  TF1* Run6Eff::GetResY06() {
    TF1* funcpp = new TF1("resY06", "[0]+[1]*x", 0.0, maxPtpp);
    
    funcpp->SetParameters(0.016, 0.013);

    return funcpp;
  }
  
} // namespace dijetcore
