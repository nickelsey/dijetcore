#include <iostream>

#include "TLorentzVector.h"

// see if root was located and linked
// against properly

int main() {
  
  TLorentzVector a( 1, 2, 3, 4 );
  
  int ret = 0;
  if (a.Px() != 1) ret = 1;
  if (a.Py() != 2) ret = 1;
  if (a.Pz() != 3) ret = 1;
  if (a.E() != 4) ret = 1;
  
  return ret;
}
