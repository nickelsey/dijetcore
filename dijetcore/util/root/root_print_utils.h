
#ifndef DIJETCORE_UTIL_ROOT_ROOT_PRINT_UTILS_H
#define DIJETCORE_UTIL_ROOT_ROOT_PRINT_UTILS_H

#include <iostream>
#include <vector>
#include <string>

#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/logging.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TROOT.h"

#include "boost/filesystem.hpp"

namespace dijetcore {
  
  struct histogramOpts {
    double label_size_x;
    double title_size_x;
    double title_offset_x;
    bool   center_title_x;
    
    double label_size_y;
    double title_size_y;
    double title_offset_y;
    bool   center_title_y;
    
    int line_width;
    double marker_size;
    
    // pick a set of colors and markers
    int current;
    std::vector<int> colors;
    std::vector<int> markers;
    
    histogramOpts() :
    colors({kBlack, kRed, kBlue, kGreen, kCyan, kMagenta, kOrange, kYellow, kRed+2, kGreen+3, kBlue-7}),
    markers({kFullCircle, kFullSquare, kFullDiamond, kFullCrossX}){
      label_size_x = 0.060;
      title_size_x = 0.075;
      title_offset_x = 0.800;
      center_title_x = false;
      label_size_y = 0.055;
      title_size_y = 0.075;
      title_offset_y = 0.800;
      center_title_y = true;
      current = 0;
      line_width = 2;
      marker_size = 1.5;
    }
    
    template <typename H>
    void SetHistogram(H* h) {
      h->GetXaxis()->SetLabelSize(label_size_x);
      h->GetXaxis()->SetTitleSize(title_size_x);
      h->GetXaxis()->SetTitleOffset(title_offset_x);
      h->GetXaxis()->CenterTitle(center_title_x);
      h->GetYaxis()->SetLabelSize(label_size_y);
      h->GetYaxis()->SetTitleSize(title_size_y);
      h->GetYaxis()->SetTitleOffset(title_offset_y);
      h->GetYaxis()->CenterTitle(center_title_y);
      
      h->SetLineWidth(line_width);
      h->SetMarkerSize(marker_size);
      h->SetMarkerStyle(markers[current % markers.size()]);
      h->SetMarkerColor(colors[current % colors.size()]);
      h->SetLineColor(colors[current % colors.size()]);
      ++current;
    }
  };
  
  struct canvasOpts {
    double left_margin;
    double bottom_margin;
    double right_margin;
    double upper_margin;
    
    bool do_legend;
    double leg_upper_bound;
    double leg_left_bound;
    double leg_lower_bound;
    double leg_right_bound;
    
    bool log_x;
    bool log_y;
    bool log_z;
    
    canvasOpts() {
      left_margin = 0.13;
      bottom_margin = 0.15;
      right_margin = 0.10;
      upper_margin = 0.10;
      
      do_legend = true;
      leg_upper_bound = 0.88;
      leg_right_bound = 0.88;
      leg_lower_bound = 0.65;
      leg_left_bound = 0.65;
      
      log_x = false;
      log_y = false;
      log_z = false;
    }
    
    template <typename C>
    void SetMargins(C* c) {
      c->SetLeftMargin(left_margin);
      c->SetRightMargin(right_margin);
      c->SetBottomMargin(bottom_margin);
      c->SetTopMargin(upper_margin);
    }
    
    template <typename C>
    void SetLogScale(C* c) {
      c->SetLogx(log_x);
      c->SetLogy(log_y);
      c->SetLogz(log_z);
    }
    
    TLegend* Legend() {
      if (do_legend) {
        return new TLegend(leg_left_bound, leg_lower_bound,
                           leg_right_bound, leg_upper_bound);
      }
      else
        return nullptr;
    }
    
  };
  
  template<typename H>
  void PrettyPrint1D(H* h,
                     histogramOpts hopts,
                     canvasOpts copts,
                     std::string hist_title,
                     std::string output_loc,
                     std::string output_name,
                     std::string canvas_title,
                     std::string x_axis_label,
                     std::string y_axis_label,
                     std::string legend_title = "") {
    // we assume the output location exists, so create
    // the final output string that will be used for pdf creation
    std::string canvas_name = output_loc + "/" + output_name + ".pdf";
    
    // and axis labels, and title
    hopts.SetHistogram(h);
    h->GetXaxis()->SetTitle(x_axis_label.c_str());
    h->GetYaxis()->SetTitle(y_axis_label.c_str());
    
    // generate a canvas
    TCanvas c;
    copts.SetMargins(&c);
    copts.SetLogScale(&c);
    
    h->Draw();
    
    TLegend* leg = copts.Legend();
    if (leg != nullptr) {
      leg->SetHeader(legend_title.c_str());
      leg->AddEntry(h, hist_title.c_str(), "lep");
      leg->Draw();
    }
    
    c.SaveAs(canvas_name.c_str());
  }
  
  template<typename H>
  void Overlay1D(const std::vector<H*>& h,
                 std::vector<std::string> hist_titles,
                 histogramOpts hopts,
                 canvasOpts copts,
                 std::string output_loc,
                 std::string output_name,
                 std::string canvas_title,
                 std::string x_axis_label,
                 std::string y_axis_label,
                 std::string legend_title = "",
                 bool find_good_range = false) {
    
    // first, check that there is a name for each histogram
    if (h.size() != hist_titles.size()) {
      std::cerr << "incorrect number of histogram names for given set of histograms" << std::endl;
      std::cerr << "for canvas: " << canvas_title << ", exiting" << std::endl;
      return;
    }
    
    // we assume the output location exists, so create
    // the final output string that will be used for pdf creation
    std::string canvas_name = output_loc + "/" + output_name + ".pdf";
    
    // set axis labels on zeroth histogram
    h[0]->GetXaxis()->SetTitle(x_axis_label.c_str());
    h[0]->GetYaxis()->SetTitle(y_axis_label.c_str());
    
    // generate a canvas
    TCanvas c;
    copts.SetMargins(&c);
    copts.SetLogScale(&c);
    
    if (find_good_range) {
      double lowest, highest;
      h[0]->GetMinimumAndMaximum(lowest, highest);
      for (int i = 1; i < h.size(); ++i) {
        if (h[i]->GetMinimum() < lowest)
          lowest = h[i]->GetMinimum();
        if (h[i]->GetMaximum() > highest)
          highest = h[i]->GetMaximum();
      }
      double delta = highest - lowest;
      double range_min = lowest - delta * 0.1;
      double range_max = highest + delta * 0.1;
      if (range_min <= 0 && copts.log_y)
        range_min = lowest;
      h[0]->GetYaxis()->SetRangeUser(range_min, range_max);
    }
    
    // print histograms, giving them some nominal settings to differentiate them
    for (int i = 0; i < h.size(); ++i) {
      hopts.SetHistogram(h[i]);
      
      if (i == 0)
        h[i]->Draw();
      else
        h[i]->Draw("SAME");
    }
    
    TLegend* leg = copts.Legend();
    if (leg != nullptr) {
      leg->SetTextSize(0.04);
      leg->SetHeader(legend_title.c_str());
      for (int i = 0; i < h.size(); ++i) {
        leg->AddEntry(h[i], hist_titles[i].c_str(), "lep");
      }
      leg->Draw();
    }
    
    c.SaveAs(canvas_name.c_str());
  }
  
  template<typename H>
  void Overlay1D(H* h1,
                 H* h2,
                 std::string h1_title,
                 std::string h2_title,
                 histogramOpts hopts,
                 canvasOpts copts,
                 std::string output_loc,
                 std::string output_name,
                 std::string canvas_title,
                 std::string x_axis_label,
                 std::string y_axis_label,
                 std::string legend_title = "") {
    // we assume the output location exists, so create
    // the final output string that will be used for pdf creation
    std::string canvas_name = output_loc + "/" + output_name + ".pdf";
    
    // axis labels, and title
    h1->SetTitle(canvas_title.c_str());
    h1->GetXaxis()->SetTitle(x_axis_label.c_str());
    h1->GetYaxis()->SetTitle(y_axis_label.c_str());
    
    // set draw options
    hopts.SetHistogram(h1);
    hopts.SetHistogram(h2);
    
    TCanvas c;
    copts.SetMargins(&c);
    copts.SetLogScale(&c);
    
    h1->Draw();
    h2->Draw("SAME");
    
    TLegend* leg = copts.Legend();
    if (leg != nullptr) {
      leg->SetHeader(legend_title.c_str());
      leg->AddEntry(h1, h1_title.c_str(), "lep");
      leg->AddEntry(h2, h2_title.c_str(), "lep");
      leg->Draw();
    }
    
    c.SaveAs(canvas_name.c_str());
  }
  
  template<typename H>
  void Overlay1D(H* h1,
                 H* h2,
                 double y_min,
                 double y_max,
                 double x_min,
                 double x_max,
                 std::string h1_title,
                 std::string h2_title,
                 histogramOpts hopts,
                 canvasOpts copts,
                 std::string output_loc,
                 std::string output_name,
                 std::string canvas_title,
                 std::string x_axis_label,
                 std::string y_axis_label,
                 std::string legend_title = "") {
    // we assume the output location exists, so create
    // the final output string that will be used for pdf creation
    std::string canvas_name = output_loc + "/" + output_name + ".pdf";
    
    // axis labels, and title
    h1->SetTitle(canvas_title.c_str());
    h1->GetXaxis()->SetTitle(x_axis_label.c_str());
    h1->GetYaxis()->SetTitle(y_axis_label.c_str());
    
    // set draw options
    hopts.SetHistogram(h1);
    hopts.SetHistogram(h2);
    
    TCanvas c;
    copts.SetMargins(&c);
    copts.SetLogScale(&c);
    
    h1->Draw();
    h2->Draw("SAME");
    
    TLegend* leg = copts.Legend();
    if (leg != nullptr) {
      leg->SetHeader(legend_title.c_str());
      leg->AddEntry(h1, h1_title.c_str(), "lep");
      leg->AddEntry(h2, h2_title.c_str(), "lep");
      leg->Draw();
    }
    
    c.SaveAs(canvas_name.c_str());
  }
  
  template<typename H>
  void Overlay1D(H* h1,
                 H* h2,
                 TGraphErrors* sys,
                 double y_min,
                 double y_max,
                 double x_min,
                 double x_max,
                 std::string h1_title,
                 std::string h2_title,
                 std::string sys_title,
                 histogramOpts hopts,
                 canvasOpts copts,
                 std::string output_loc,
                 std::string output_name,
                 std::string canvas_title,
                 std::string x_axis_label,
                 std::string y_axis_label,
                 std::string legend_title = "") {
    // we assume the output location exists, so create
    // the final output string that will be used for pdf creation
    std::string canvas_name = output_loc + "/" + output_name + ".pdf";
    
    // axis labels, and title
    h1->SetTitle(canvas_title.c_str());
    h1->GetXaxis()->SetTitle(x_axis_label.c_str());
    h1->GetYaxis()->SetTitle(y_axis_label.c_str());
    h1->GetYaxis()->SetRangeUser(y_min, y_max);
    h1->GetXaxis()->SetRangeUser(x_min, x_max);
    
    // set draw options
    hopts.SetHistogram(h1);
    hopts.SetHistogram(h2);
    
    TCanvas c;
    copts.SetMargins(&c);
    copts.SetLogScale(&c);
    
    sys->SetFillColorAlpha(h2->GetLineColor(), 0.6);
    sys->SetFillStyle(1001);
    sys->SetLineWidth(0);
    sys->SetMarkerSize(0);
    
    h1->Draw("9");
    h2->Draw("9SAME");
    sys->Draw("9e2SAME");
    
    TLegend* leg = copts.Legend();
    if (leg != nullptr) {
      leg->SetHeader(legend_title.c_str());
      leg->AddEntry(h1, h1_title.c_str(), "lep");
      leg->AddEntry(h2, h2_title.c_str(), "lep");
      leg->AddEntry(sys,sys_title.c_str(), "f");
      leg->Draw();
    }
    
    c.SaveAs(canvas_name.c_str());
  }
  
  template<typename H>
  void Print2DSimple(H* h,
                     histogramOpts hopts,
                     canvasOpts copts,
                     std::string output_loc,
                     std::string output_name,
                     std::string canvas_title,
                     std::string x_axis_label,
                     std::string y_axis_label,
                     std::string opt = "COLZ") {
    // we assume the output location exists, so create
    // the final output string that will be used for pdf creation
    std::string canvas_name = output_loc + "/" + output_name + ".pdf";
    
    // and axis labels, and title
    h->SetTitle(canvas_title.c_str());
    h->GetXaxis()->SetTitle(x_axis_label.c_str());
    h->GetYaxis()->SetTitle(y_axis_label.c_str());
    hopts.SetHistogram(h);
    
    TCanvas c;
    copts.SetLogScale(&c);
    copts.SetMargins(&c);
    
    h->Draw(opt.c_str());
    c.SaveAs(canvas_name.c_str());
  }
  
  // print histograms & their ratios
  template <typename H>
  void PrintWithRatio(H* h1,
                      H* h2,
                      std::string h1_title,
                      std::string h2_title,
                      histogramOpts hopts,
                      canvasOpts copts,
                      std::string output_loc,
                      std::string output_name,
                      std::string canvas_title,
                      std::string x_axis_label,
                      std::string y_axis_label,
                      std::string legend_title = "") {
    
    // we assume the output location exists, so create
    // the final output string that will be used for pdf creation
    std::string canvas_name = output_loc + "/" + output_name + ".pdf";
    
    TCanvas c;
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.35, 1.0, 1.0);
    copts.SetMargins(pad1);
    pad1->SetBottomMargin(0.0);
    copts.SetLogScale(pad1);
    pad1->Draw();
    pad1->cd();
    
    // set histograms
    hopts.SetHistogram(h1);
    h1->GetXaxis()->SetTitle("");
    h1->GetYaxis()->SetTitle(y_axis_label.c_str());
    h1->SetTitle(canvas_title.c_str());
    h1->Draw();
    
    hopts.SetHistogram(h2);
    h2->Draw("SAME");
    
    TLegend* leg = copts.Legend();
    if (leg != nullptr) {
      leg->SetTextSize(0.04);
      leg->SetHeader(legend_title.c_str());
      leg->AddEntry(h1, h1_title.c_str(), "lep");
      leg->AddEntry(h2, h2_title.c_str(), "lep");
      leg->Draw();
    }
    
    // lower pad
    c.cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);
    copts.SetMargins(pad2);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.35);
    copts.SetLogScale(pad2);
    pad2->SetLogy(false);
    pad2->Draw();
    pad2->cd();
    
    TH1D* tmp = (TH1D*) h1->Clone();
    tmp->Divide(h2);
    hopts.SetHistogram(tmp);
    
    tmp->GetYaxis()->SetRangeUser(0.7, 1.3);
    tmp->GetYaxis()->SetNdivisions(8);
    tmp->GetXaxis()->SetTitle(x_axis_label.c_str());
    tmp->GetYaxis()->SetTitle("Ratio");
    tmp->SetLineColor(h1->GetLineColor());
    tmp->SetMarkerColor(h1->GetMarkerColor());
    tmp->SetMarkerStyle(h1->GetMarkerStyle());
    
    tmp->GetXaxis()->SetLabelSize(hopts.label_size_x*2);
    tmp->GetXaxis()->SetTitleSize(hopts.title_size_x*2);
    tmp->GetXaxis()->SetTitleOffset(hopts.title_offset_x);
    tmp->GetXaxis()->CenterTitle(hopts.center_title_x);
    tmp->GetYaxis()->SetLabelSize(hopts.label_size_y*2);
    tmp->GetYaxis()->SetTitleSize(hopts.title_size_y*2);
    tmp->GetYaxis()->SetTitleOffset(hopts.title_offset_y/2);
    tmp->GetYaxis()->CenterTitle(hopts.center_title_y);
    
    tmp->Draw("ep");
    
    c.SaveAs(canvas_name.c_str());
  }

  // convenience wrapper to PrintWithRatio above
  template<class H>
  void PrintWithRatio(std::vector<H*>& histograms, std::vector<string>& legend_entries,
                      dijetcore::histogramOpts hopts, dijetcore::canvasOpts copts,
                      string loc, string name, string x_axis, string y_axis, 
                      string leg_title = "") {
    if (histograms.size() != 2 ) {
      LOG(ERROR) << "wrong number of histograms, PrintRatios requires 2";
      return;
    }
    if (legend_entries.size() != 2) {
      LOG(ERROR) << "wrong number of legend entries, PrintRatios requires 2";
      return ;
    }
    dijetcore::PrintWithRatio(histograms[0], histograms[1], legend_entries[0],
                              legend_entries[1], hopts, copts, loc, name, "",
                              x_axis, y_axis, leg_title);
};
  
  
  // print histograms & their ratios
  template <typename H>
  void PrintWithRatios(H* ref,
                       std::vector<H*> h,
                       std::string ref_name,
                       std::vector<std::string> h_titles,
                       histogramOpts hopts,
                       canvasOpts copts,
                       std::string output_loc,
                       std::string output_name,
                       std::string canvas_title,
                       std::string x_axis_label,
                       std::string y_axis_label,
                       std::string legend_title = "") {
    if (h.size() == 0 ) {
      std::cerr << "warning: only one histogram, cant take ratios, exiting" << std::endl;
      return;
    }
    if (h.size() + 1 != h_titles.size()) {
      std::cerr << "warning, mismatched number of histograms & histogram names, exiting" << std::endl;
      return;
    }
    
    // we assume the output location exists, so create
    // the final output string that will be used for pdf creation
    std::string canvas_name = output_loc + "/" + output_name + ".pdf";
    
    TCanvas c;
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1.0, 1.0);
    copts.SetMargins(pad1);
    pad1->SetBottomMargin(0.0);
    copts.SetLogScale(pad1);
    pad1->Draw();
    pad1->cd();
    
    // set histograms
    hopts.SetHistogram(ref);
    ref->GetXaxis()->SetTitle("");
    ref->GetYaxis()->SetTitle(y_axis_label.c_str());
    ref->SetTitle(canvas_title.c_str());
    ref->Draw();
    
    for (auto hist : h) {
      hopts.SetHistogram(hist);
      hist->Draw("SAME");
    }
    
    TLegend* leg = copts.Legend();
    if (leg != nullptr) {
      leg->SetTextSize(0.04);
      leg->SetHeader(legend_title.c_str());
      leg->AddEntry(ref, ref_name.c_str(), "lep");
      for (int i = 0; i < h.size(); ++i) {
        leg->AddEntry(h[i], h_titles[i].c_str(), "lep");
      }
      leg->Draw();
    }
    
    // lower pad
    c.cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    copts.SetMargins(pad2);
    pad2->SetTopMargin(0.0);
    copts.SetLogScale(pad2);
    pad2->SetLogx(false);
    pad2->Draw();
    pad2->cd();
    
    for (int i = 0; i < h.size(); ++i) {
      std::string tmp_name = "tmp_" + std::to_string(i);
      TH1D* tmp = (TH1D*) ref->Clone(tmp_name.c_str());
      tmp->Divide(h[i]);
      tmp->GetXaxis()->SetTitle(x_axis_label.c_str());
      tmp->GetYaxis()->SetTitle("Ratio");
      tmp->SetLineColor(h[i]->GetLineColor());
      tmp->SetMarkerColor(h[i]->GetMarkerColor());
      tmp->SetMarkerStyle(h[i]->GetMarkerStyle());
      if (i == 0)
        tmp->Draw();
      else
        tmp->Draw("SAME");
    }
    c.SaveAs(canvas_name.c_str());
  }  

} // namespace dijetcore

#endif // DIJETCORE_UTIL_ROOT_ROOT_PRINT_UTILS_H
