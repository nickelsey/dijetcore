#include "dijetcore/worker/dijet_worker/eta_cone_worker.h"

#include "dijetcore/lib/assert.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/util/data/vector_conversion.h"
#include "dijetcore/util/fastjet/selector_compare.h"

#include <algorithm>
#include <vector>

#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

namespace dijetcore {

EtaConeWorker::EtaConeWorker() {}

EtaConeWorker::~EtaConeWorker() {}

std::unordered_map<std::string, unique_ptr<EtaConeOutput>> &
EtaConeWorker::run(DijetWorker &worker,
                   std::vector<fastjet::PseudoJet> &input) {

  // clear our output
  eta_cone_result_.clear();

  // see if we need to do anything at all -
  // if there was no match in the Dijetworker, we exit out
  bool found_match = false;
  for (auto &result : worker.dijets())
    if (result.second->found_match)
      found_match = true;

  if (found_match == false)
    return eta_cone_result_;

  // get the input container
  auto &dijets = worker.dijets();

  for (auto &dijet : dijets) {

    if (!dijet.second->found_match)
      continue;

    auto &key = dijet.first;
    auto &def = dijet.second->dijet_def;
    auto &lead_hard_jet = dijet.second->lead_hard;
    auto &sub_hard_jet = dijet.second->sublead_hard;
    auto &lead_match_jet = dijet.second->lead_match;
    auto &sub_match_jet = dijet.second->sublead_match;

    // only select constituents that were used during jetfinding for the matched
    // jets.
    std::vector<fastjet::PseudoJet> lead_constituents =
        def->lead->matchedJetDef().constituentSelector()(input);
    std::vector<fastjet::PseudoJet> sub_constituents =
        def->sub->matchedJetDef().constituentSelector()(input);

    // and only select particles that have a pt less than the particles used for
    // hard-core constituent selection. This way we don't "double count" random
    // hard-core fluctuations
    double max_pt_lead = ExtractDoubleFromSelector(
        def->lead->initialJetDef().constituentSelector(), "pt >=");
    double max_pt_sub = ExtractDoubleFromSelector(
        def->sub->initialJetDef().constituentSelector(), "pt >=");
    fastjet::Selector max_pt_sel_lead = fastjet::SelectorPtMax(max_pt_lead);
    fastjet::Selector max_pt_sel_sub = fastjet::SelectorPtMax(max_pt_sub);

    lead_constituents = max_pt_sel_lead(lead_constituents);
    sub_constituents = max_pt_sel_sub(sub_constituents);

    // find an acceptable eta based on jet radius and position within the event
    double max_const_eta_lead = ExtractDoubleFromSelector(
        def->lead->matchedJetDef().constituentSelector(), "|eta| <=");
    std::pair<bool, double> lead_cone_eta = findEta(
        lead_match_jet, def->lead->matchedJetDef().R(), max_const_eta_lead);

    double max_const_eta_sub = ExtractDoubleFromSelector(
        def->sub->matchedJetDef().constituentSelector(), "|eta| <=");
    std::pair<bool, double> sub_cone_eta = findEta(
        sub_match_jet, def->sub->matchedJetDef().R(), max_const_eta_sub);

    // if either jet could not fit an eta cone, we continue
    if (lead_cone_eta.first == false || sub_cone_eta.first == false) {
      unique_ptr<EtaConeOutput> res = make_unique<EtaConeOutput>();
      eta_cone_result_[key] = std::move(res);
      continue;
    }

    // valid locations exist for both the lead and the sublead jet - now we will
    // do the cone "jetfinding"

    unique_ptr<EtaConeOutput> res =
        make_unique<EtaConeOutput>(lead_match_jet, sub_match_jet);

    eta_cone_result_[key] = std::move(res);
  }

  return eta_cone_result_;
}

std::pair<bool, double> EtaConeWorker::findEta(const fastjet::PseudoJet &reference, double r,
                                double const_max_eta) {
  // check that it is possible to fit another cone into the constituent region
  if (2.0 * r >= const_max_eta) {
    LOG(ERROR) << "jet radius " << r
               << " is too large to fit an eta offset cone into a the "
                  "constituent region. (max constituent eta="
               << const_max_eta << ")";
    return {false, 0.0};
  }

  double jet_max_eta = const_max_eta - r;
  double ref_eta = reference.eta();

  // build the two possible eta regions where it can be placed
  double negative_window_low = -jet_max_eta;
  double negative_window_high = ref_eta - 2.0 * r;

  double positive_window_low = ref_eta + 2.0 * r;
  double positive_window_high = jet_max_eta;

  // randomly shuffle so that if both windows are valid, we don't always choose
  // the same one
  std::vector<std::pair<double, double>> windows{
      {negative_window_low, negative_window_high},
      {positive_window_low, positive_window_high}};
  std::shuffle(windows.begin(), windows.end(), gen_);

  // now loop over both windows: if the range [low, high] is ordered correctly -
  // that is, if high > low, then that means there is a valid range available
  // for the random cone. Otherwise, it means the two restrictions (from max
  // constituent eta & jet axis) overlap, and there is no available room for the
  // eta cone
  for (auto &window : windows) {
    if (window.second - window.first < 0.0)
      continue;
    std::uniform_real_distribution<> eta_range(window.first, window.second);
    double eta = eta_range(gen_);
    return {true, eta};
  }

  // neither window was valid - return failure
  return {false, 0.0};
}

} // namespace dijetcore
