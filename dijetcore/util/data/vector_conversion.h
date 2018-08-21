#ifndef DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H
#define DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H

// vector_conversion functions to convert
// between TLorentzVector, fastjet::PseudoJet,
// etc

#include <random>

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
    TStarJetVector* sv;
    for(int i = 0; i < container->GetEntries() ; ++i) {
      sv = container->Get(i);
      
      fastjet::PseudoJet tmpPJ = fastjet::PseudoJet(*sv);
      tmpPJ.set_user_index(sv->GetCharge());
      particles[offset + i] = tmpPJ;
    }
  }
  
  void ConvertPythiaToPseudoJet(const Pythia8::Pythia& gen, std::vector<fastjet::PseudoJet>& particles) {
    
    size_t offset = particles.size();
    size_t final_state = 0;
    particles.resize(offset + gen.event.size());
    
    for (size_t i = 0; i < gen.event.size(); ++i) {
      if ( gen.event[i].isFinal() && gen.event[i].isVisible() ) {
        fastjet::PseudoJet tmp(gen.event[i].px(), gen.event[i].py(), gen.event[i].pz(), gen.event[i].e());
        tmp.set_user_index(gen.event[i].charge());
        particles[offset + final_state] = tmp;
        ++final_state;
      }
    }
    particles.resize(offset + final_state);
  }
  
} // dijetcore

#endif // DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H
