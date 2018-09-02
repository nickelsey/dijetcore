#ifndef DIJETCORE_UTIL_DATA_CENTRALITY_CENTRALITY_RUN7_H
#define DIJETCORE_UTIL_DATA_CENTRALITY_CENTRALITY_RUN7_H

// Run 7 centrality

#include <vector>
#include <random>

namespace dijetcore {
  
  class CentralityRun7 {
    public:
    CentralityRun7();
    ~CentralityRun7();

    
    // given a corrected refmult, calculate the 16 or 9 bin centrality
    int Centrality9(unsigned grefmult);
    
    std::vector<unsigned> CentralityBounds9Bin() const {return cent_bin_9_;}
    
    private:

    std::vector<unsigned> cent_bin_9_;
    
  };
  
} // namespace dijetcore

#endif // DIJETCORE_UTIL_DATA_CENTRALITY_CENTRALITY_RUN7_HH

