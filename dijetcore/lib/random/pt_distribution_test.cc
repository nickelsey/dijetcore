#include "gtest/gtest.h"

#include <random>

#include "TF1.h"
#include "TH1.h"

#include "dijetcore/lib/random/pt_distribution.h"

TEST(PtDistribution, BoundsTest) {
    std::mt19937 gen(14342);

    dijetcore::pt_distribution<double> dis;

    for(int i = 0; i < 100000; ++i) {
        double generated = dis(gen);
        EXPECT_LE(generated, 5.0);
        EXPECT_GE(generated, 0.0);
    }
}

TEST(PtDistribution, FitTest) {
    std::mt19937 gen(14342);

    dijetcore::pt_distribution<double> dis;
  
    TH1D* h = new TH1D("h", "", 100, 0, 5);

    TF1* f = new TF1("f", "[0]*x*exp(-x/[1])", 0, 5);
    f->SetParameters(1, 0.291);
    
    for (int i = 0; i < 100000; ++i) {
        h->Fill(dis(gen));
    }
    
    h->Scale(1.0 / h->Integral());
    h->Fit("f");

    EXPECT_LE(f->GetChisquare()/f->GetNDF(), 2.0);
}
