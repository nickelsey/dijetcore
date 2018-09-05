#ifndef DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H
#define DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H

// vector_conversion functions to convert
// between TLorentzVector, fastjet::PseudoJet,
// etc

#include "fastjet/PseudoJet.hh"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "Pythia8/Pythia.h"

namespace dijetcore {
  
  // converts a TStarJetVectorContainer of TStarJetVectors
  // into an std::vector of fastjet::PseudoJets
  void ConvertTStarJetVector(TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet>& particles) {
    
    // resize container
    size_t offset = particles.size();
    particles.resize(particles.size() + container->GetEntries());
    
    // Transform TStarJetVectors into (FastJet) PseudoJets
    // ---------------------------------------------------
    for(int i = 0; i < container->GetEntries() ; ++i) {
      TStarJetVector* sv = container->Get(i);
      
      fastjet::PseudoJet tmpPJ = fastjet::PseudoJet(*sv);
      tmpPJ.set_user_index(sv->GetCharge());
      particles[offset + i] = tmpPJ;
    }
  }
 
} // dijetcore

#endif // DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H
