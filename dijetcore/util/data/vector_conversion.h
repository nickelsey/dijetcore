#ifndef DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H
#define DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H

// vector_conversion functions to convert
// between TLorentzVector, fastjet::PseudoJet,
// etc

#include <vector>

#include "TStarJetVectorContainer.h"

class TStarJetVector;
namespace fastjet {
class PseudoJet;
}

namespace dijetcore {

// converts a TStarJetVectorContainer of TStarJetVectors
// into an std::vector of fastjet::PseudoJets
void ConvertTStarJetVector(TStarJetVectorContainer<TStarJetVector> *container,
                           std::vector<fastjet::PseudoJet> &particles);

} // namespace dijetcore

#endif // DIJETCORE_UTIL_DATA_VECTOR_CONVERSION_H
