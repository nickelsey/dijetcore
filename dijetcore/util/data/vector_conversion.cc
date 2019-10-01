#include "dijetcore/util/data/vector_conversion.h"

#include "fastjet/PseudoJet.hh"
#include "TStarJetVector.h"

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