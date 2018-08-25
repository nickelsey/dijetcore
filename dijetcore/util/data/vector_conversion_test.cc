#include "gtest/gtest.h"

#include "dijetcore/util/data/vector_conversion.h"

#include "fastjet/PseudoJet.hh"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"

TEST(VectorConversion, ConvertToPseudoJet) {
    TStarJetVectorContainer<TStarJetVector> container;
    TLorentzVector A(1.0, 2.0, 3.0, 4.0);
    TStarJetVector vecA(A);
    vecA.SetCharge(1);
    TLorentzVector B(2.5, 3.5, 1.0, 3);
  TStarJetVector vecB(B);
    vecB.SetCharge(-1);
    container.Add(&vecA);
    container.Add(&vecB);

    std::vector<fastjet::PseudoJet> container_out;
    dijetcore::ConvertTStarJetVector(&container, container_out);

    EXPECT_NEAR(vecA.E(), container_out[0].E(), 1e-4);
    EXPECT_NEAR(vecA.Pt(), container_out[0].pt(), 1e-4);
    EXPECT_EQ(vecA.GetCharge(), container_out[0].user_index());
    EXPECT_NEAR(vecB.E(), container_out[1].E(), 1e-4);
    EXPECT_NEAR(vecB.Pt(), container_out[1].pt(), 1e-4);
    EXPECT_EQ(vecB.GetCharge(), container_out[1].user_index());
    
}
