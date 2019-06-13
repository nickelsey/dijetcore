#ifndef DIJETCORE_UTIL_DIJET_IMBALANCE_PRINT_ROUTINES_H
#define DIJETCORE_UTIL_DIJET_IMBALANCE_PRINT_ROUTINES_H

#include "dijetcore/lib/types.h"
#include "dijetcore/util/root/root_print_utils.h"

#include "TLatex.h"
#include "TPaveText.h"

namespace dijetcore {

template <typename H>
void AjPrintout(H *h1, H *h2, TGraphErrors *sys, double y_min, double y_max,
                double x_min, double x_max, TPaveText &text, string h1_title,
                string h2_title, dijetcore::histogramOpts hopts,
                dijetcore::canvasOpts copts, string output_loc,
                string output_name, string canvas_title, string x_axis_label,
                string y_axis_label, string legend_title = "") {
  // we assume the output location exists, so create
  // the final output string that will be used for pdf creation
  string canvas_name = output_loc + "/" + output_name + ".pdf";

  // axis labels, and title
  h1->SetTitle(canvas_title.c_str());
  h1->GetXaxis()->SetTitle(x_axis_label.c_str());
  h1->GetYaxis()->SetTitle(y_axis_label.c_str());
  h1->GetYaxis()->SetRangeUser(y_min, y_max);
  h1->GetXaxis()->SetRangeUser(x_min, x_max);

  // set draw options
  hopts.SetHistogram(h1);
  hopts.SetHistogram(h2);

  if (sys == nullptr) {
    h2->SetLineColor(kBlue);
    h2->SetMarkerColor(kBlue);
  }

  TCanvas c;
  copts.SetMargins(&c);
  copts.SetLogScale(&c);

  h1->Draw("9");
  h2->Draw("9SAME");

  if (sys != nullptr) {
    sys->SetFillColorAlpha(h2->GetLineColor(), 0.45);
    sys->SetFillStyle(1001);
    sys->SetLineWidth(0);
    sys->SetMarkerSize(0);
    sys->Draw("9e2SAME");
  }

  TLegend *leg = copts.Legend();
  if (leg != nullptr) {
    leg->SetHeader(legend_title.c_str());
    leg->SetTextFont(62);
    leg->AddEntry(h1, h1_title.c_str(), "p")->SetTextSize(.035);
    leg->AddEntry(h2, h2_title.c_str(), "p")->SetTextSize(.035);
    leg->Draw();
  }

  // and STAR preliminary message
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.065);
  latex.SetTextColor(kRed + 3);
  // latex.DrawLatex(0.19, 0.8, "STAR Preliminary");

  text.Draw();

  c.SaveAs(canvas_name.c_str());
}

// gives us a vector of TPaveText axis labels for the top and right side of the
// grid (jet radius and pT const, respectively)
std::vector<TPaveText> GetGridLabelsForAj(
    unsigned nbins_x, unsigned nbins_y, double min_x, double max_x,
    double min_y, double max_y, double x_low_axis_margin = 0.1,
    double x_high_axis_margin = 0.1, double y_low_axis_margin = 0.1,
    double y_high_axis_margin = 0.1, double text_size_x = 0.03,
    double text_size_y = 0.03);

} // namespace dijetcore

#endif // DIJETCORE_UTIL_DIJET_IMBALANCE_PRINT_ROUTINES_H