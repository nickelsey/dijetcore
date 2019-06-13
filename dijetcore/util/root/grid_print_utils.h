#ifndef DIJETCORE_UTIL_ROOT_GRID_PRINT_UTILS_H
#define DIJETCORE_UTIL_ROOT_GRID_PRINT_UTILS_H

#include <string>
#include <vector>

#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/memory.h"
#include "dijetcore/lib/string/string_utils.h"

#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TROOT.h"

#include "boost/filesystem.hpp"

namespace dijetcore {

// object that can be passed to grid print functions to print text (in the form
// of a TPaveText or TLegend) in any grid entry - location 0 x 0 is the bottom
// left corner
struct GridTextObject {
  unsigned location_x = 0;
  unsigned location_y = 0;
  unique_ptr<TPaveText> text;
  unique_ptr<TLegend> legend;

  GridTextObject(){};

  GridTextObject(const GridTextObject &rhs)
      : location_x(rhs.location_x), location_y(rhs.location_y) {
    if (rhs.text.get())
      text = make_unique<TPaveText>(*rhs.text.get());
    if (rhs.legend.get())
      legend = make_unique<TLegend>(*rhs.legend.get());
  }
};

// structure with layer specific grid print options
struct GridPrintOptions {
  bool layer_1_active = true;
  bool layer_2_active = true;
  bool layer_3_active = true;
  std::string layer_1_print = "";
  std::string layer_2_print = "SAME";
  std::string layer_3_print = "SAME";
};

// partitions a TCanvas into an equal sized grid of Nx x Ny  with the specified
// margins
std::vector<std::vector<TPad *>>
CanvasPartition(TCanvas *C, const Int_t Nx, const Int_t Ny, Float_t lMargin,
                Float_t rMargin, Float_t bMargin, Float_t tMargin);

template <class H1, class H2, class H3>
void PrintGrid(TCanvas *C, const std::vector<std::vector<TPad *>> &pads,
               const std::vector<std::vector<H1 *>> &layer_1,
               const std::vector<std::vector<H2 *>> &layer_2,
               const std::vector<std::vector<H3 *>> &layer_3,
               const GridPrintOptions &opts,
               const std::vector<GridTextObject> &text, TPad *extra_pad,
               const std::vector<TPaveText> &extra_pad_text,
               const std::string &output_file) {

  // loop over every canvas and print the histograms/graphs
  for (int x = 0; x < pads.size(); ++x) {
    for (int y = 0; y < pads[x].size(); ++y) {
      pads[x][y]->cd();
      // also have to scale  the tick lengths so each pad is uniform
      double x_factor = pads[0][0]->GetAbsWNDC() / pads[x][y]->GetAbsWNDC();
      double y_factor = pads[0][0]->GetAbsHNDC() / pads[x][y]->GetAbsHNDC();
      bool axis_drawn = false;

      if (opts.layer_1_active) {
        if (!axis_drawn) {
          axis_drawn = true;
          layer_1[x][y]->GetYaxis()->SetTickLength(x_factor * 0.04 / y_factor);
          layer_1[x][y]->GetXaxis()->SetTickLength(y_factor * 0.06 / x_factor);
        }
        layer_1[x][y]->Draw(opts.layer_1_print.c_str());
      }
      if (opts.layer_2_active) {
        if (!axis_drawn) {
          axis_drawn = true;
          layer_2[x][y]->GetYaxis()->SetTickLength(x_factor * 0.04 / y_factor);
          layer_2[x][y]->GetXaxis()->SetTickLength(y_factor * 0.06 / x_factor);
        }
        layer_2[x][y]->Draw(opts.layer_2_print.c_str());
      }
      if (opts.layer_3_active) {
        if (!axis_drawn) {
          axis_drawn = true;
          layer_3[x][y]->GetYaxis()->SetTickLength(x_factor * 0.04 / y_factor);
          layer_3[x][y]->GetXaxis()->SetTickLength(y_factor * 0.06 / x_factor);
        }
        layer_3[x][y]->Draw(opts.layer_3_print.c_str());
      }
    }
  }
  
  // now print TPaveText in individual canvases
  for (auto &t : text) {
    pads[t.location_x][t.location_y]->cd();
    if (t.legend.get() != nullptr)
      t.legend->Draw();
    if (t.text.get() != nullptr)
      t.text->Draw();
  }

  // print larger, cross-pad text on an invisible pad
  if (extra_pad != nullptr) {
    C->cd();
    extra_pad->Draw();
    extra_pad->cd();
    for (auto &t : extra_pad_text)
      t.DrawClone();
  }

  C->cd();
  C->SaveAs(output_file.c_str());
}

} // namespace dijetcore

#endif // DIJETCORE_UTIL_ROOT_GRID_PRINT_UTILS_H