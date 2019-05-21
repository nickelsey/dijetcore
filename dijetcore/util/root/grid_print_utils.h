#ifndef DIJETCORE_UTIL_ROOT_GRID_PRINT_UTILS_H
#define DIJETCORE_UTIL_ROOT_GRID_PRINT_UTILS_H

#include <string>
#include <vector>

#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TROOT.h"

#include "boost/filesystem.hpp"

namespace dijetcore {

// partitions a TCanvas into an equal sized grid of Nx x Ny  with the specified
// margins
std::vector<std::vector<TPad *>>
CanvasPartition(TCanvas *C, const Int_t Nx, const Int_t Ny, Float_t lMargin,
                Float_t rMargin, Float_t bMargin, Float_t tMargin) {
  if (!C)
    return std::vector<std::vector<TPad *>>();
  std::vector<std::vector<TPad *>> ret(Nx, std::vector<TPad *>(Ny));
  // Setup Pad layout:
  Float_t vSpacing = 0.0;
  Float_t vStep = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny;

  Float_t hSpacing = 0.0;
  Float_t hStep = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx;

  Float_t vposd, vposu, vmard, vmaru, vfactor;
  Float_t hposl, hposr, hmarl, hmarr, hfactor;

  for (Int_t i = 0; i < Nx; i++) {

    if (i == 0) {
      hposl = 0.0;
      hposr = lMargin + hStep;
      hfactor = hposr - hposl;
      hmarl = lMargin / hfactor;
      hmarr = 0.0;
    } else if (i == Nx - 1) {
      hposl = hposr + hSpacing;
      hposr = hposl + hStep + rMargin;
      hfactor = hposr - hposl;
      hmarl = 0.0;
      hmarr = rMargin / (hposr - hposl);
    } else {
      hposl = hposr + hSpacing;
      hposr = hposl + hStep;
      hfactor = hposr - hposl;
      hmarl = 0.0;
      hmarr = 0.0;
    }

    for (Int_t j = 0; j < Ny; j++) {

      if (j == 0) {
        vposd = 0.0;
        vposu = bMargin + vStep;
        vfactor = vposu - vposd;
        vmard = bMargin / vfactor;
        vmaru = 0.0;
      } else if (j == Ny - 1) {
        vposd = vposu + vSpacing;
        vposu = vposd + vStep + tMargin;
        vfactor = vposu - vposd;
        vmard = 0.0;
        vmaru = tMargin / (vposu - vposd);
      } else {
        vposd = vposu + vSpacing;
        vposu = vposd + vStep;
        vfactor = vposu - vposd;
        vmard = 0.0;
        vmaru = 0.0;
      }

      C->cd(0);

      std::string name = MakeString(C->GetName(), "pad_", i, "_", j).c_str();
      TPad *pad = (TPad *)gROOT->FindObject(name.c_str());
      if (pad)
        delete pad;
      pad = new TPad(name.c_str(), "", hposl, vposd, hposr, vposu);
      pad->SetLeftMargin(hmarl);
      pad->SetRightMargin(hmarr);
      pad->SetBottomMargin(vmard);
      pad->SetTopMargin(vmaru);

      pad->SetFrameBorderMode(0);
      pad->SetBorderMode(0);
      pad->SetBorderSize(0);

      pad->Draw();

      ret[i][j] = pad;
    }
  }
  return ret;
}

} // namespace dijetcore

#endif // DIJETCORE_UTIL_ROOT_GRID_PRINT_UTILS_H