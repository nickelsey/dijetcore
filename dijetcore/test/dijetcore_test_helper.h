#ifndef DIJETCORE_TEST_DIJETCORE_TEST_HELPER_H
#define DIJETCORE_TEST_DIJETCORE_TEST_HELPER_H

#include "dijetcore/lib/logging.h"
#include "dijetcore/util/fastjet/selector_compare.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_definition.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/dijet_matrix.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/jet_def.h"
#include "dijetcore/worker/dijet_worker/dijet_matrix/match_def.h"
#include "dijetcore/worker/dijet_worker/dijet_worker.h"

#include "gtest/gtest.h"

namespace dijetcore {
namespace testing {
void CheckJetDef(dijetcore::JetDef def, fastjet::JetAlgorithm alg,
                 double radius, double jet_pt, double const_pt,
                 double const_eta, fastjet::RecombinationScheme scheme,
                 fastjet::Strategy strategy, fastjet::AreaType area_type,
                 int ghost_repeat, double ghost_area, double grid_scatter,
                 double pt_scatter, double mean_ghost_pt,
                 fastjet::JetDefinition bkg_def,
                 fastjet::AreaDefinition bkg_area);

void CheckMatchDef(dijetcore::MatchDef *def, fastjet::JetAlgorithm alg,
                   double radius, double match_radius, double jet_pt,
                   double const_init_pt, double const_match_pt,
                   double const_eta, fastjet::RecombinationScheme scheme,
                   fastjet::Strategy strategy, fastjet::AreaType area_type,
                   int ghost_repeat, double ghost_area, double grid_scatter,
                   double pt_scatter, double mean_ghost_pt,
                   fastjet::JetDefinition bkg_def,
                   fastjet::AreaDefinition bkg_area);

void CheckDijetDefinition(
    dijetcore::DijetDefinition *def, fastjet::JetAlgorithm lead_alg,
    fastjet::JetAlgorithm sub_alg, double lead_R, double lead_match_R,
    double sub_R, double sub_match_R, double lead_pt, double sub_pt,
    double lead_const_init_pt, double lead_const_match_pt,
    double sub_const_init_pt, double sub_const_match_pt, double const_eta,
    fastjet::RecombinationScheme scheme, fastjet::Strategy strategy,
    fastjet::AreaType area_type, int ghost_repeat, double ghost_area,
    double grid_scatter, double pt_scatter, double mean_ghost_pt,
    fastjet::JetDefinition bkg_def, fastjet::AreaDefinition bkg_area_lead,
    fastjet::AreaDefinition bkg_area_sub);

} // namespace testing
} // namespace dijetcore

#endif // DIJETCORE_TEST_DIJETCORE_TEST_HELPER_H