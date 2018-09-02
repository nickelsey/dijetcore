#include "dijetcore/util/data/efficiency/run7_eff.h"

namespace dijetcore {
  
  Run7Eff::Run7Eff(std::string filename) : file(nullptr), maxPt(5.0), maxPtpp(9.0), sys_(TrackingUncY7::NONE) {
    loadFile(filename);
  }
  
  Run7Eff::~Run7Eff() {
    delete effY06;
    
    delete effY04[0];
    delete effY04[1];
    delete effY04[2];
    
    delete effY07pteta[0];
    delete effY07pteta[1];
    delete effY07pteta[2];
    
    delete effY07eta[0];
    delete effY07eta[1];
    delete effY07eta[2];
  }
  
  void Run7Eff::loadFile(std::string filename) {
    if (file)
    file->Close();
    file = new TFile(filename.c_str(), "READ");
    
    loadCurves();
  }
  
  double Run7Eff::AuAuEff(double pt, double eta, int cent) {
    Double_t effWeight=1.0;
    if(pt < 5.)
    effWeight = effY04[cent]->Eval(eta, pt);
    else
    effWeight = effY04[cent]->Eval(eta, 5.0);
    if(pt > 1.5)
    effWeight *= effY07eta[cent]->GetBinContent(effY07eta[cent]->GetXaxis()->FindBin(eta));
    else
    effWeight *= effY07pteta[cent]->GetBinContent(effY07pteta[cent]->GetXaxis()->FindBin(pt),effY07pteta[cent]->GetYaxis()->FindBin(eta));
    
    return effWeight;
  }
  
  double Run7Eff::AuAuEff020Avg(double pt, double eta) {
    double eff05 = AuAuEff(pt, eta, 0);
    double eff10 = AuAuEff(pt, eta, 1);
    double eff20 = AuAuEff(pt, eta, 2);
    return (eff05 + eff10 + 2 * eff20) / 4.0;
  }
  
  double Run7Eff::ppEff(double pt, double eta) {
    Double_t effWeight=1.0;
    
    effWeight = effY06->Eval(pt, eta);
    
    return effWeight;
  }
  
  void Run7Eff::loadCurves() {
    effY07pteta[0] = (TH2D*)file->Get("ptetaScale_0");
    effY07pteta[0]->SetName("effY07pteta_0");
    effY07pteta[1] = (TH2D*)file->Get("ptetaScale_1");
    effY07pteta[1]->SetName("effY07pteta_1");
    effY07pteta[2] = (TH2D*)file->Get("ptetaScale_2");
    effY07pteta[2]->SetName("effY07pteta_2");
    effY07eta[0] = (TH1D*)file->Get("etaScale_0");
    effY07eta[0]->SetName("effY07eta_0");
    effY07eta[1] = (TH1D*)file->Get("etaScale_1");
    effY07eta[1]->SetName("effY07eta_1");
    effY07eta[2] = (TH1D*)file->Get("etaScale_2");
    effY07eta[2]->SetName("effY07eta_2");
    effY04[0] = GetEffY04(0);
    effY04[0]->SetName("effY04_0");
    effY04[1] = GetEffY04(1);
    effY04[1]->SetName("effY04_1");
    effY04[2] = GetEffY04(2);
    effY04[2]->SetName("effY04_2");
    
    effY06 = GetEffY06();
    effY06->SetName("effY06");
  }
  
  TF2* Run7Eff::GetEffY06()
  {
    
    TF2* funcpp = new TF2("ppEfficiency","[0]-0.06-[1]*exp([2]/x)+[3]*exp(-0.5*((x-[4])/[5])**2)/sqrt(2*pi*[5]*[5])-[6]*exp(-0.5*((x-[7])/[8])**2)/sqrt(2*pi*[8]*[8])+([9]-[10]*(y-[11])^2-[12]*(y-[11])^4-[13]*(y-[11])^6-[14]*(y-[11])^8)*exp(-[15]*x)",0.,10.,-1.,1.);
    
    Double_t parset[] = {0.869233,0.0223402,0.44061,0.558762,0.145162,0.508033,110.008,-4.63659,1.73765,0.0452674,-0.101279,0.0081551,0.945287,-2.00949,1.61746,1.39352};
    
    ((TF2*)funcpp)->SetParameters(parset);
    
    return funcpp;
  }
  
  
  
  TF2* Run7Eff::GetEffY04(Int_t cb)
  {
    Double_t fMaxPtPara = 5;
    
    char name[30];
    sprintf(name,"EfficiencyFunction%i",cb);
    
    TF2* func = new TF2(name,"[0]+[1]*x^2+[2]*x^4+[3]*x^6+[4]*x^8+[5]*exp([6]*y)+[7]*y+[8]*y*y +  [9]*exp(-((abs(x)-[10])^2)  /[11] - ((  abs(y)-[12]  )^2) /[13]) ",-1.,1.,0.,fMaxPtPara);
    
    Double_t parset4[]={ 0.698809, 0.0347652, -0.00819333, 0.112736, -0.356907, -1.62506, -7.26695, 0.0436162, -0.00453185, 0.249514, 0.308879, 0.133046, 0.295414, 0.0019349 };
    Double_t parset5[]={ 0.652242, 0.0760859, -0.0784171, 0.0393619, -0.247293, -2.04786, -8.96039, 0.0603416, -0.006971, -0.0435101, 0.131005, 0.00053132, 0.74369, 0.00576589 };
    Double_t parset6[]={ 0.631911, 0.117639, -0.29002, 0.522928, -0.569609, -1.3921, -6.73044, 0.0588101, -0.00686795, 0.110982, 0.2951, 0.14493, 0.295612, 0.00290843 };
    
    if(cb == 0)
      ((TF2*)func)->SetParameters(parset6);
    else if(cb == 1)
      ((TF2*)func)->SetParameters(parset5);
    else if(cb == 2)
      ((TF2*)func)->SetParameters(parset4);
    
    if(cb < 0 || cb > 2)
      std::cout << "Error: Nonsensical Centrality Class!" << std::endl;
    
    return func;
  }
  
  double Run7Eff::ratio(double pt, double eta, int cent) {
    Double_t pp_eff = ppEff(pt, eta);
    Double_t auau_eff = AuAuEff(pt, eta, cent);
    
    double ratio_ = auau_eff / pp_eff;
    
    int sign_ = static_cast<int>(sys_);
    if (sign_ == 0)
      return ratio_;
    
    double systematic_ = sign_ * ratioUncertainty(pt, eta);
    
    ratio_ = ratio_ + ratio_ * systematic_;
    
    if (ratio_ > 1.0)
      return 1.0;
    
    return ratio_;
  }
  
  double Run7Eff::ratio020Avg(double pt, double eta) {
    Double_t pp_eff = ppEff(pt, eta);
    Double_t auau_eff = AuAuEff020Avg(pt, eta);
    
    double ratio_ = auau_eff / pp_eff;
    
    int sign_ = static_cast<int>(sys_);
    if (sign_ == 0)
      return ratio_;
    
    double systematic_ = sign_ * ratioUncertainty(pt, eta);
    
    ratio_ = ratio_ + ratio_ * systematic_;
    
    if (ratio_ > 1.0)
      return 1.0;
    
    return ratio_;
  }
  
  
  double Run7Eff::ratioUncertainty(double pt, double eta) {
    
    double cent_u_ = 0.04;
    double auau_u_  = 0.03;
    double pp_u_  = 0.03;
    
    double auau_ = AuAuEff020Avg(pt, eta);
    double pp_ = ppEff(pt, eta);
    double ratio_ = auau_ / pp_;
    
    double au_term_ = pow(auau_u_ / auau_, 2.0);
    double pp_term_ = pow(pp_u_ / pp_, 2.0);
    double cent_term_ = pow(cent_u_ / auau_, 2.0);
    return sqrt(au_term_ + pp_term_ + cent_term_);
  }
  
} // namespace dijetcore
