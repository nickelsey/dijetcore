#include "dijetcore/util/data/efficiency/run14_eff_scaled.h"

namespace dijetcore {
  
  Run14EffScaled::Run14EffScaled(std::string filename) : file(nullptr), max_pt_(4.5), sys_(TrackingUnc::NONE),
  auau_u_(0.05), pp_u_(0.03), cent_u_(0.00) {
    if (!filename.empty())
    loadFile(filename);
  }
  
  Run14EffScaled::~Run14EffScaled() {
    
  }
  
  void Run14EffScaled::loadFile(std::string filename, int nBinsZDC, int nBinsCent) {
    if (file.get())
    file->Close();
    
    file = std::make_shared<TFile>(filename.c_str(), "READ");
    loadCurves(nBinsZDC, nBinsCent);
    
    // load year 6
    Double_t parset[] = {0.869233,0.0223402,0.44061,0.558762,0.145162,0.508033,110.008,-4.63659,1.73765,0.0452674,-0.101279,0.0081551,0.945287,-2.00949,1.61746,1.39352};
    effY06 = std::make_shared<TF2>("ppEfficiency","[0]-0.06-[1]*exp([2]/x)+[3]*exp(-0.5*((x-[4])/[5])**2)/sqrt(2*pi*[5]*[5])-[6]*exp(-0.5*((x-[7])/[8])**2)/sqrt(2*pi*[8]*[8])+([9]-[10]*(y-[11])^2-[12]*(y-[11])^4-[13]*(y-[11])^6-[14]*(y-[11])^8)*exp(-[15]*x)",0.,10.,-1.,1.);
    effY06->SetParameters(parset);
  }
  
  double Run14EffScaled::AuAuEff(double pt, double eta, int cent, double zdcrate) {
    if (pt > max_pt_)
      pt = max_pt_;
    
    int zdcBin = luminosityBin(zdcrate);
    
    int bin = curves.at(zdcBin).at(cent)->FindBin(pt, eta);
    double eff = curves.at(zdcBin).at(cent)->GetBinContent(bin);
    double scalar = GetScalar(pt, cent);
    return eff / scalar;
  }
  
  double Run14EffScaled::pp6Eff(double pt, double eta) {
    if (abs(eta) > 1.0) 
    std::cerr << "warning: efficiency curves only valid for |eta| < 1.0" << std::endl;
    if (pt > 10.0)
    pt = 10.0;
    return effY06->Eval(pt, eta);
  }
  
  double Run14EffScaled::ratio(double pt, double eta, int cent, double zdcrate) {
    double ratio_ = AuAuEff(pt, eta, cent, zdcrate) / pp6Eff(pt, eta);
    
    int sign_ = static_cast<int>(sys_);
    if (sign_ == 0)
    return ratio_;
    
    double systematic_ = sign_ * ratioUncertainty(pt, eta, cent, zdcrate);
    
    ratio_ = ratio_ + ratio_ * systematic_;
    
    if (ratio_ > 1.0)
    return 1.0;
    
    return ratio_;
  }
  
  double Run14EffScaled::ratioUncertainty(double pt, double eta, int cent, double zdcrate) {
    double auau_ = AuAuEff(pt, eta, cent, zdcrate);
    double pp_ = pp6Eff(pt, eta);
    double ratio_ = auau_ / pp_;
    
    double au_term_ = pow(auau_u_ / auau_, 2.0);
    double pp_term_ = pow(pp_u_ / pp_, 2.0);
    double cent_term_ = pow(cent_u_ / auau_, 2.0);
    return sqrt(au_term_ + pp_term_ + cent_term_);
  }
  
  int Run14EffScaled::luminosityBin(double zdcrate) {
    int zdcBin;
    if (zdcrate/1000.0 < 33.0)
    zdcBin = 0;
    else if (zdcrate/1000.0 < 66.0)
    zdcBin = 1;
    else
    zdcBin = 2;
    
    return zdcBin;
  }
  
  double Run14EffScaled::GetScalar(double pt, int cent) {
    if (pt > max_pt_)
      pt = max_pt_;
    
    int scalar_cent_bin = -1;
    double scalar = 1.0;
    if (curves.size() > 0) {
      if (curves[0].size() == 16) {
        if (cent == 0 || cent == 1)
          scalar_cent_bin = cent;
        else if (cent == 2 || cent == 3)
          scalar_cent_bin = 2;
      }
      else {
        scalar_cent_bin = cent;
      }
    }
    if (scalar_cent_bin >= scalars.size())
      scalar = 1.0;
    else
      scalar = scalars[scalar_cent_bin]->GetBinContent(scalars[scalar_cent_bin]->FindBin(pt));
    
    if (scalar > 1.0 || scalar < 0.0)
      scalar = 1.0;
    
    return scalar;
  }
  
  void Run14EffScaled::loadCurves(int nBinsZDC, int nBinsCent) {
    curves.clear();
    for (int i = 0; i < nBinsZDC; ++i) {
      curves.push_back(std::vector<TH2D*>());
      for (int j = 0; j < nBinsCent; ++j) {
        std::string name = "efficiency_lumi_" + std::to_string(i)
        + "_cent_" + std::to_string(j);
        curves[i].push_back((TH2D*) file->Get(name.c_str()));
      }
    }
  }
  
  void Run14EffScaled::loadScalars(int nBinsScalar) {
    scalars.clear();
    for (int i = 0; i < nBinsScalar; ++i) {
      std::string name = "eff_ratio_" + std::to_string(i);
      scalars.push_back((TH1D*) file->Get(name.c_str()));
    }
  }
  
} // namespace dijetcore
