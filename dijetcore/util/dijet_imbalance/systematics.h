#ifndef DIJETCORE_UTIL_DIJET_IMBALANCE_SYSTEMATICS_H
#define DIJETCORE_UTIL_DIJET_IMBALANCE_SYSTEMATICS_H

#include "TGraphErrors.h"

namespace dijetcore {

template<class T1>
TGraphErrors *GetSystematic(T1 *nom, T1 *var1_a, T1 *var1_b, T1 *var2_a,
                            T1 *var2_b) {
  int nBins = nom->GetNbinsX();
  double x_[nBins];
  double y_[nBins];
  double x_err_[nBins];
  double y_err_[nBins];

  for (int i = 0; i < nBins; ++i) {
    x_[i] = nom->GetBinCenter(i + 1);
    y_[i] = nom->GetBinContent(i + 1);
    x_err_[i] = nom->GetXaxis()->GetBinWidth(1) / 2.0;
    double diff_var_1_a =
        fabs(nom->GetBinContent(i + 1) - var1_a->GetBinContent(i + 1));
    double diff_var_1_b =
        fabs(nom->GetBinContent(i + 1) - var1_b->GetBinContent(i + 1));
    double diff_var_2_a =
        fabs(nom->GetBinContent(i + 1) - var2_a->GetBinContent(i + 1));
    double diff_var_2_b =
        fabs(nom->GetBinContent(i + 1) - var2_b->GetBinContent(i + 1));
    double max_var_1 =
        (diff_var_1_a > diff_var_1_b ? diff_var_1_a : diff_var_1_b);
    double max_var_2 =
        (diff_var_2_a > diff_var_2_b ? diff_var_2_a : diff_var_2_b);
    y_err_[i] = sqrt(max_var_1 * max_var_1 + max_var_2 * max_var_2);
  }
  TGraphErrors *ret = new TGraphErrors(nBins, x_, y_, x_err_, y_err_);
  return ret;
}

template <class T1, class T2>
TH1D *GetFractionalError(T1 *nom, T2 *variation, std::string ext = "errfrac") {
  if (nom->GetXaxis()->GetNbins() != variation->GetXaxis()->GetNbins()) {
    LOG(ERROR) << "bin mismatch: exiting";
    return nullptr;
  }

  string name = nom->GetName() + ext;
  TH1D *ret = new TH1D(name.c_str(), "", nom->GetXaxis()->GetNbins(),
                       nom->GetXaxis()->GetXmin(), nom->GetXaxis()->GetXmax());
  ret->GetYaxis()->SetRangeUser(0.000001, 0.599999);
  for (int i = 1; i < nom->GetXaxis()->GetNbins(); ++i) {
    double bin_nom = nom->GetBinContent(i);
    double bin_sys = variation->GetBinContent(i);
    if (bin_nom != 0)
      ret->SetBinContent(i, fabs(bin_nom - bin_sys) / bin_nom);
    else
      ret->SetBinContent(i, 0.0);
  }
  return ret;
}

template <class T1>
TH1D *GetFractionalError(T1 *nom, TGraphErrors *errors,
                         std::string ext = "errfrac") {
  if (nom->GetXaxis()->GetNbins() != errors->GetN()) {
    LOG(ERROR) << "bin mismatch: exiting";
    return nullptr;
  }

  string name = nom->GetName() + ext;
  TH1D *ret = new TH1D(name.c_str(), "", nom->GetXaxis()->GetNbins(),
                       nom->GetXaxis()->GetXmin(), nom->GetXaxis()->GetXmax());
  ret->GetYaxis()->SetRangeUser(0.000001, 0.599999);
  for (int i = 1; i < nom->GetXaxis()->GetNbins(); ++i) {
    double bin_content = nom->GetBinContent(i);
    double bin_error = errors->GetErrorY(i - 1);
    if (bin_content != 0)
      ret->SetBinContent(i, bin_error / bin_content);
    else
      ret->SetBinContent(i, 0.0);
  }
  return ret;
}

} // namespace dijetcore

#endif // DIJETCORE_UTIL_DIJET_IMBALANCE_SYSTEMATICS_H