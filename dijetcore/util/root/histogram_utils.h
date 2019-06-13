#ifndef DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H
#define DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H

#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/types.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace dijetcore {

// normalizes the histogram such that h->Integral() == 1.0, modifies h in place
// and returns ptr to h
template <class H> H *Norm(H *h) {
  h->Scale(1.0 / h->Integral());
  return h;
}

// returns the bin-by-bin difference between two 1D histograms
template <class H1, class H2> TH1D *Diff(H1 *h1, H2 *h2) {
  if (h1->GetNbinsX() != h2->GetNbinsX() ||
      h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax() ||
      h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin()) {
    LOG(ERROR) << "can't subtract histograms with different bins";
  }
  string name = MakeString(h1->GetName(), "minus", h2->GetName());
  TH1D *diff = new TH1D(name.c_str(), "", h1->GetNbinsX(),
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

// returns the bin-by-bin relative difference between two 1D histograms:
// (h1 - h2)/ h1, and gives the option of setting the y axis range (does not
// rescale y axis if ymin == ymax == 0)
template <class H1, class H2>
TH1D *RelativeDiff(H1 *h1, H2 *h2, double ymin = 0.0, double ymax = 0.0) {
  TH1D *diff = Diff(h1, h2);
  diff->Divide(h1);
  if (ymin != 0.0 || ymax != 0.0)
    diff->GetYaxis()->SetRangeUser(ymin, ymax);
  return diff;
}

// splits a 2D histogram into TH1D histograms along the requested axis.
// name_prefix is an optional argument to append to the name of each individual
// histogram to stop ROOT problems with equal histogram names
template <class H1>
std::vector<TH1D *> Split2DByBin(H1 *h, int split_on_x = true,
                                 string name_suffix = "") {
  std::vector<TH1D *> ret;
  if (split_on_x) {
    for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
      string name = string(h->GetName()) + std::to_string(i) + name_suffix;
      TH1D *tmp = h->ProjectionY(name.c_str(), i, i);
      ret.push_back(tmp);
    }
    h->GetXaxis()->SetRange();
  } else {
    for (int i = 1; i <= h->GetYaxis()->GetNbins(); ++i) {
      string name = string(h->GetName()) + std::to_string(i) + name_suffix;
      TH1D *tmp = h->ProjectionX(name.c_str(), i, i);
      ret.push_back(tmp);
    }
    h->GetYaxis()->SetRange();
  }
  return ret;
}

// splits a 3D histogram into TH2D histograms along the requested axis (x=1,
// y=2, z=3). name_prefix is an optional argument to append to the name of each
// individual histogram to stop ROOT problems with equal histogram names
template <class H1>
std::vector<TH2D *> Split3DByBin(H1 *h, int split_on_axis = 1,
                                 string name_suffix = "") {
  std::vector<TH2D *> ret;
  if (split_on_axis == 1) {
    for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
      string name = string(h->GetName()) + std::to_string(i) + name_suffix;
      h->GetXaxis()->SetRange(i, i);
      TH2D *tmp = (TH2D *)h->Project3D("zy");
      tmp->SetName(name.c_str());
      ret.push_back(tmp);
    }
    h->GetXaxis()->SetRange();
  } else if (split_on_axis == 2) {
    for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
      string name = string(h->GetName()) + std::to_string(i) + name_suffix;
      h->GetYaxis()->SetRange(i, i);
      TH2D *tmp = (TH2D *)h->Project3D("zx");
      tmp->SetName(name.c_str());
      ret.push_back(tmp);
    }
    h->GetYaxis()->SetRange();
  } else if (split_on_axis == 3) {
    for (int i = 1; i <= h->GetXaxis()->GetNbins(); ++i) {
      string name = string(h->GetName()) + std::to_string(i) + name_suffix;
      h->GetZaxis()->SetRange(i, i);
      TH2D *tmp = (TH2D *)h->Project3D("yx");
      tmp->SetName(name.c_str());
      ret.push_back(tmp);
    }
    h->GetZaxis()->SetRange();
  } else
    LOG(ERROR) << "requested invalid axis: " << split_on_axis;
  return ret;
}

// splits a 2D histogram into 1D (see Split2DByBin() above) and normalizes each
// of the resulting TH1Ds so that their integral == 1.0 (see Norm() above).
// h should be a 2D histogram, split_on_x selects either x or y axis to split
// on, and name_suffix is appended to each individual TH1D's name to make it
// easier to avoid ROOT object name clashes
template <class H1>
std::vector<TH1D *> Split2DByBinNormalized(H1 *h, int split_on_x = true,
                                           string name_suffix = "") {
  auto ret = Split2DByBin(h, split_on_x, name_suffix);
  for (auto hist : ret)
    Norm(hist);
  return ret;
}

} // namespace dijetcore

#endif // DIJETCORE_UTIL_ROOT_HISTOGRAM_UTILS_H
