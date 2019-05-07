#include "dijetcore/util/mc/jewel_reader.h"

#include "dijetcore/lib/assert.h"
#include "dijetcore/lib/logging.h"

namespace dijetcore {

JewelReader::JewelReader(std::string &filename)
    : input_(filename.c_str(), "read"), event_("T", &input_),
      header_("runInfo", &input_), xsec_(header_, "crossSection"),
      total_weight_(header_, "totalWeight"), nevents_(header_, "nEvents"),
      pid_(event_, "particleID"), status_(event_, "status"), px_(event_, "px"),
      py_(event_, "py"), pz_(event_, "pz"), e_(event_, "energy"), m_(event_, "mass"),
      vx_(event_, "vertex_x"), vy_(event_, "vertex_y"), w_(event_, "weight") {
        header_.SetEntry(0);
      }

bool JewelReader::next() {
  clear();
  if (!event_.Next())
    return false;

  create_pseudojets();
  return true;
}

bool JewelReader::read(unsigned idx) {
  clear();
  if (event_.SetEntry(idx) != TTreeReader::kEntryValid)
    return false;
  
  create_pseudojets();
  return true;
}

void JewelReader::create_pseudojets() {
  clear();
  LOG(INFO) << "new event";
  for (int i = 0; i < px_.GetSize(); ++i) {
    LOG(INFO) << "particle: " << i;
    LOG(INFO) << "particle status: " << status_[i];
    LOG(INFO) << "particle pid: " << pid_[i];
    LOG(INFO) << "particle e: " << e_[i];
    if (status_[i] != 0)
      continue;
    fastjet::PseudoJet tmp(px_[i], py_[i], pz_[i], e_[i]);
    tmp.set_user_index(pid_[i]);
    processed_.push_back(tmp);
  }
}

} // namespace dijetcore