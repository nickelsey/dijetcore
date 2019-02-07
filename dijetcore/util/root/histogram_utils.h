#ifndef DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H
#define DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

namespace dijetcore {

  template<typename T>
  void Mean1D(T* h, double& mean, double& counts, double minx= -1e10, double maxx = 1e10) {
    int binlow = h->GetXaxis()->FindBin(minx);
    int binhigh = h->GetXaxis()->FindBin(maxx);
    if (binlow < 1) binlow = 1;
    if (binhigh > h->GetXaxis()->GetNbins()) binhigh = h->GetXaxis()->GetNbins();
    
    double sum = 0;
    
    for (int i = binlow; i <= binhigh; ++i) {
      sum += h->GetBinContent(i)*h->GetBinCenter(i);
      counts += h->GetBinContent(i);
    }
    mean = sum / counts;
  }
  
  template<typename T>
  void MeanStdY(T* h, double& mean, double& std, double minx= -1e10, double maxx = 1e10) {
    int binlow = h->GetXaxis()->FindBin(minx);
    int binhigh = h->GetXaxis()->FindBin(maxx);
    if (binlow < 1) binlow = 1;
    if (binhigh > h->GetXaxis()->GetNbins()) binhigh = h->GetXaxis()->GetNbins();
    
    double sum = 0;
    double counts = 0;;
    
    for (int i = binlow; i <= binhigh; ++i) {
      // we dont want dead towers to skew the results
      if ( h->GetBinContent(i) == 0) continue;
      sum += h->GetBinContent(i);
      counts++;
    }
    mean = sum / counts;
    
    sum = 0;
    for (int i = binlow; i <= binhigh; ++i) {
      // we dont want dead towers to skew the results
      if ( h->GetBinContent(i) == 0) continue;
      sum += pow(h->GetBinContent(i) - mean,2);
    }
    std = sqrt(sum/counts);
    
  }
  
  void MeanStdDev2D(TH2F* h, double& mean, double& std, double miny= -1e10, double maxy = 1e10) {
    
    double sum = 0;
    double counts = 0;
    std::vector<std::pair<double,double> > binMeans(h->GetNbinsX(), {0.0, 0.0});
    
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      h->GetXaxis()->SetRange(i,i);
      TH1D* tmp = h->ProjectionY();
      if (tmp->GetEntries() == 0) continue;
      double tmpMean = 0;
      double tmpCounts  = 0;
      Mean1D(tmp, tmpMean, tmpCounts, miny, maxy);
      binMeans[i-1] = std::pair<double,double>{tmpMean, tmpCounts};
      
      sum += tmpMean*tmpCounts;
      counts += tmpCounts;
      
      delete tmp;
    }
    h->GetXaxis()->SetRange();
    mean = sum / counts;
    
    sum = 0;
    counts = 0;
    
    // now calculate std dev
    for (int i = 0; i < binMeans.size(); ++i){
      if (binMeans[i].second <= 0) continue;
      sum += (binMeans[i].first - mean) * (binMeans[i].first - mean) * binMeans[i].second;
      counts += binMeans[i].second;
    }
    std = sqrt(sum / counts);
  }
  
  template<typename T>
  void PrintWithBounds(T* p, double mean, double std, double nsigma, std::string outname) {
    
    TCanvas c1;
    p->Draw();
    
    TF1* tmp1 = new TF1("tmp1", "[0]", 0, 20000);
    tmp1->SetParameter(0, mean - std * nsigma);
    TF1* tmp2 = new TF1("tmp2", "[0]", 0, 20000);
    tmp2->SetParameter(0, mean + std * nsigma);
    
    tmp1->Draw("same");
    tmp2->Draw("same");
    
    c1.SaveAs(outname.c_str());
    
    delete tmp1;
    delete tmp2;
  }

} // namespace dijetcore
  
#endif // DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H
