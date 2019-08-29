#ifndef DIJETCORE_UTIL_MC_JEWEL_READER_H
#define DIJETCORE_UTIL_MC_JEWEL_READER_H

#include "dijetcore/lib/types.h"

#include <vector>

#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "dijetcore/lib/root_ext/TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

namespace dijetcore {

class JewelReader {
public:
  JewelReader(std::string &filename);

  bool next();
  bool read(unsigned idx);

  std::vector<fastjet::PseudoJet> &event() { return processed_; }
  float vx() { return *vx_; }
  float vy() { return *vy_; }
  float weight() { return *w_; }

  float totalXsec() { return *xsec_; }
  float totalWeight() {return *total_weight_;}
  unsigned nEvents() { return event_.GetTree()->GetEntries(); }

  fastjet::PseudoJet leadingParton() {return lead_p_;}
  fastjet::PseudoJet subParton() {return sub_p_;}

private:
  // clears the container of pseudojets. Called by read() and next()
  void clear() { processed_.clear(); }

  // processes the current event
  void create_pseudojets();

  TFile input_;
  TTreeReader event_;
  TTreeReader header_;

  // branches for the trees
  TTreeReaderValue<double> xsec_;
  TTreeReaderValue<double> total_weight_;
  TTreeReaderValue<int> nevents_;

  TTreeReaderArray<short> pid_;
  TTreeReaderArray<char> status_;
  TTreeReaderArray<float> px_;
  TTreeReaderArray<float> py_;
  TTreeReaderArray<float> pz_;
  TTreeReaderArray<float> e_;
  TTreeReaderArray<float> m_;

  TTreeReaderValue<float> vx_;
  TTreeReaderValue<float> vy_;
  TTreeReaderValue<double> w_;

  // event container
  std::vector<fastjet::PseudoJet> processed_;
  fastjet::PseudoJet lead_p_;
  fastjet::PseudoJet sub_p_;
};

} // namespace dijetcore

#endif // DIJETCORE_UTIL_MC_JEWEL_READER_H