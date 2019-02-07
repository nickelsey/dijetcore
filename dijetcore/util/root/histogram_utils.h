#ifndef DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H
#define DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H

#include "dijetcore/lib/logging.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

namespace dijetcore {

template <class H>
H* Norm(H* h) {
  h->Scale(1.0 / h->Integral());
  return h;
}

template <class H1, class H2>
TH1D* Diff(H1* h1, H2* h2) {
  if (h1->GetNbinsX() != h2->GetNbinsX() ||
      h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax() ||
      h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin()) {
    LOG(ERROR) << "can't subtract histograms with different bins";
  }
  string name = MakeString(h1->GetName(), "minus", h2->GetName());
  TH1D* diff = new TH1D(name.c_str(), "", h1->GetNbinsX(),
                        h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());

  for (int i = 1; i <= h1->GetNbinsX(); ++i) {
    double binContent = h1->GetBinContent(i) - h2->GetBinContent(i);
    double error1 = h1->GetBinError(i);
    double error2 = h2->GetBinError(i);
    double binError = sqrt(pow(error1, 2.0) + pow(error2, 2.0));

    diff->SetBinContent(i, binContent);
    diff->SetBinError(i, binError);
  }
  return diff;
}

template <class H1, class H2>
TH1D* RelativeDiff(H1* h1, H2* h2, double ymin = -1.0, double ymax = 1.0) {
  TH1D* diff = Diff(h1, h2);
  diff->Divide(h1);
  diff->GetYaxis()->SetRangeUser(ymin, ymax);
  return diff;
}

}  // namespace dijetcore

#endif  // DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H
