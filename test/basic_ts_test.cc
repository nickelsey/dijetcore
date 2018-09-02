#include <iostream>

#include "TLorentzVector.h"

#include "TStarJetVector.h"

// see if TStarJetPico library was found
// and being linked against properly

int main() {
  
  TLorentzVector v(1, 2, 3, 4);
  TStarJetVector a(v);

  int ret = 0;
  if (a.Px() != 1) ret = 1;
  if (a.Py() != 2) ret = 1;
  if (a.Pz() != 3) ret = 1;
  if (a.E() != 4) ret = 1;

  return ret;
}
