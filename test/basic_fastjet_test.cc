#include <iostream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

// see if fastjet was located and linked
// against properly

int main() {
  
  fastjet::PseudoJet a( 1, 2, 3, 4 );
  
  int ret = 0;
  if (a.px() != 1) ret = 1;
  if (a.py() != 2) ret = 1;
  if (a.pz() != 3) ret = 1;
  if (a.E() != 4) ret = 1;
  
  return ret;
}
