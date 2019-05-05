#ifndef DIJETCORE_UTIL_FASTJET_DIJET_KEY_H
#define DIJETCORE_UTIL_FASTJET_DIJET_KEY_H

#include <iomanip>
#include <sstream>
#include <string>

#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"

#include "dijetcore/lib/assert.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/fastjet/selector_compare.h"

namespace dijetcore {

struct DijetKey {
  double lead_init_r;
  int lead_init_alg;
  double lead_init_pt;
  double lead_match_r;
  int lead_match_alg;
  double lead_match_pt;
  double lead_init_const_pt;
  double lead_init_const_eta;
  double lead_match_const_pt;
  double lead_match_const_eta;
  double sub_init_r;
  int sub_init_alg;
  double sub_init_pt;
  double sub_match_r;
  int sub_match_alg;
  double sub_match_pt;
  double sub_init_const_pt;
  double sub_init_const_eta;
  double sub_match_const_pt;
  double sub_match_const_eta;
};

DijetKey ParseStringToDijetKey(string key);

string ParseDijetKeyToString(DijetKey &key);

string MakeKeyFromDijetDefinition(const DijetDefinition &def);

string MakeSortedKeyFromDijetDefinition(const DijetDefinition &def);

} // namespace dijetcore

#endif // DIJETCORE_LIB_FASTJET_DIJET_KEY_H
