// Aj of pythia in a detector simulation of the STAR detector
// for both generator and detector level jets

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/random.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/root/root_print_utils.h"

#include <algorithm>
#include <random>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include "RooUnfold/RooUnfold.h"
#include "RooUnfold/RooUnfoldBayes.h"

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "job", "name for output file");
DIJETCORE_DEFINE_int(id, 0, "job id (when running parallel jobs)");
DIJETCORE_DEFINE_string(config, "", "configuration file");

// define the set of parameters required in config
const std::set<string> required_params = {
    "output_directory",   // directory for output root file to written to
    "num_events",         // max number of events to read from tree
    "match_radius",       // matching radius for det/gen matching
    "lead_jet_pt",        // lead hard-core jet pt minimum det
    "lead_jet_pt_gen",    // lead hard-core jet pt minimum gen
    "sublead_jet_pt",     // sublead hard-core jet pt minimum det
    "sublead_jet_pt_gen", // sublead hard-core jet pt minimum gen
    "det_trig_threshold", // minimum detector neutral trigger threshold
    "gen_trig_threshold", // minimum generator neutral trigger threshold
    "seed",               // seed for rng, if <0 use random
    "pp_file",            // will plot a pp comparison if auau and pp is present
    "auau_file", // will plot an auau comparison if auau and pp is present
    "tree_name", // auau/pp tree name (must match)
    "force_gen_miss_with_det_cuts", // fill misses in "no_miss" responses when
                                    // generator level passes detector level
                                    // cuts
    "proj_leadpt_min",              // unfolded projection lead pt min
    "proj_leadpt_max",              // unfolded projection lead pt max
    "proj_subpt_min",               // unfolded projection sub pt min
    "proj_subpt_max",               // unfolded projection sub pt max
    "proj_aj_min",                  // unfolded projection aj min
    "proj_aj_max"                   // unfolded projection aj max
};

struct HistBinning {
  int lead_pt_bin = 10;
  int sub_pt_bin = 10;
  int aj_bin = 10;
  double lead_pt_min = 0.0;
  double lead_pt_max = 50.0;
  double sub_pt_min = 0.0;
  double sub_pt_max = 50.0;
  double aj_min = 0.0;
  double aj_max = 1.0;
};

struct ProjectionRange {
  double leadptmin = 0.0;
  double leadptmax = 50.0;
  double subptmin = 0.0;
  double subptmax = 50.0;
  double ajmin = 0.0;
  double ajmax = 1.0;

  void SetRanges(TH3D *h) {
    h->GetXaxis()->SetRangeUser(leadptmin, leadptmax);
    h->GetYaxis()->SetRangeUser(subptmin, subptmax);
    h->GetZaxis()->SetRangeUser(ajmin, ajmax);
  }
  void ResetRanges(TH3D *h) {
    h->GetXaxis()->SetRange();
    h->GetYaxis()->SetRange();
    h->GetZaxis()->SetRange();
  }
};

struct MatchEvent {
  fastjet::PseudoJet lead_hc_gen;
  fastjet::PseudoJet sub_hc_gen;
  fastjet::PseudoJet lead_match_gen;
  fastjet::PseudoJet sub_match_gen;
  fastjet::PseudoJet lead_hc_det;
  fastjet::PseudoJet sub_hc_det;
  fastjet::PseudoJet lead_match_det;
  fastjet::PseudoJet sub_match_det;
  fastjet::PseudoJet gen_trig;
  fastjet::PseudoJet det_trig;
  bool has_gen;
  bool has_det;
  bool gen_with_det_cuts;
  fastjet::PseudoJet gen1;
  fastjet::PseudoJet gen2;
  bool gen1_g;
  bool gen2_g;
};

struct ResponseSet {
  RooUnfoldResponse *full;
  RooUnfoldResponse *full_no_miss;
  RooUnfoldResponse *geom_match;
  RooUnfoldResponse *full_gluon;
  RooUnfoldResponse *full_no_miss_gluon;
  RooUnfoldResponse *geom_match_gluon;
  RooUnfoldResponse *full_quark;
  RooUnfoldResponse *full_no_miss_quark;
  RooUnfoldResponse *geom_match_quark;
};

struct HistSet {
  TH3D *full;
  TH3D *quark;
  TH3D *gluon;
};

ResponseSet BuildResponseSet(std::string name, TH3D *truth_template,
                             TH3D *det_template) {
  ResponseSet tmp;
  tmp.full = new RooUnfoldResponse(dijetcore::MakeString(name, "_full"),
                                   "3D Aj response");
  tmp.full->Setup(det_template, truth_template);
  tmp.full_no_miss = new RooUnfoldResponse(
      dijetcore::MakeString(name, "_full_no_miss"), "3D Aj response");
  tmp.full_no_miss->Setup(det_template, truth_template);
  tmp.geom_match = new RooUnfoldResponse(
      dijetcore::MakeString(name, "_geom_match"), "3D Aj response");
  tmp.geom_match->Setup(det_template, truth_template);
  tmp.full_gluon = new RooUnfoldResponse(
      dijetcore::MakeString(name, "_full_gluon"), "3D Aj response");
  tmp.full_gluon->Setup(det_template, truth_template);
  tmp.full_no_miss_gluon = new RooUnfoldResponse(
      dijetcore::MakeString(name, "_full_no_miss_gluon"), "3D Aj response");
  tmp.full_no_miss_gluon->Setup(det_template, truth_template);
  tmp.full_quark = new RooUnfoldResponse(
      dijetcore::MakeString(name, "_full_quark"), "3D Aj response");
  tmp.full_quark->Setup(det_template, truth_template);
  tmp.full_no_miss_quark = new RooUnfoldResponse(
      dijetcore::MakeString(name, "_full_no_miss_quark"), "3D Aj response");
  tmp.full_no_miss_quark->Setup(det_template, truth_template);
  tmp.geom_match_gluon = new RooUnfoldResponse(
      dijetcore::MakeString(name, "_gm_gluon"), "3D Aj response");
  tmp.geom_match_gluon->Setup(det_template, truth_template);
  tmp.geom_match_quark = new RooUnfoldResponse(
      dijetcore::MakeString(name, "_gm_quark"), "3D Aj response");
  tmp.geom_match_quark->Setup(det_template, truth_template);

  return tmp;
}

HistSet BuildHistSet(std::string name, HistBinning h) {
  HistSet tmp;

  tmp.full = new TH3D(dijetcore::MakeString(name, "_full").c_str(),
                      ";p_{T}^{lead};p_{T}^{sub};A_{J}", h.lead_pt_bin,
                      h.lead_pt_min, h.lead_pt_max, h.sub_pt_bin, h.sub_pt_min,
                      h.sub_pt_max, h.aj_bin, h.aj_min, h.aj_max);

  tmp.quark = new TH3D(dijetcore::MakeString(name, "_quark").c_str(),
                       ";p_{T}^{lead};p_{T}^{sub};A_{J}", h.lead_pt_bin,
                       h.lead_pt_min, h.lead_pt_max, h.sub_pt_bin, h.sub_pt_min,
                       h.sub_pt_max, h.aj_bin, h.aj_min, h.aj_max);
  tmp.gluon = new TH3D(dijetcore::MakeString(name, "_gluon").c_str(),
                       ";p_{T}^{lead};p_{T}^{sub};A_{J}", h.lead_pt_bin,
                       h.lead_pt_min, h.lead_pt_max, h.sub_pt_bin, h.sub_pt_min,
                       h.sub_pt_max, h.aj_bin, h.aj_min, h.aj_max);

  return tmp;
}

double Aj(double lead, double sub) { return fabs(lead - sub) / (lead + sub); }

using hist_map = std::unordered_map<std::string, TH1D *>;

hist_map Proj3D(TH3D *hist, ProjectionRange range, bool use_range = true) {
  std::unordered_map<std::string, TH1D *> ret;

  if (use_range)
    range.SetRanges(hist);
  ret["x"] = (TH1D*) hist->Project3D("x");
  ret["x"]->SetName(
      dijetcore::MakeString(hist->GetName(), "xproj",
                            dijetcore::Counter::instance().counter())
          .c_str());
  ret["x"]->Scale(1.0 / ret["x"]->Integral());
  ret["y"] = (TH1D*) hist->Project3D("y");
  ret["y"]->SetName(
      dijetcore::MakeString(hist->GetName(), "yproj",
                            dijetcore::Counter::instance().counter())
          .c_str());
  ret["y"]->Scale(1.0 / ret["y"]->Integral());
  ret["z"] = (TH1D*) hist->Project3D("z");
  ret["z"]->SetName(
      dijetcore::MakeString(hist->GetName(), "zproj",
                            dijetcore::Counter::instance().counter())
          .c_str());
  ret["z"]->Scale(1.0 / ret["z"]->Integral());
  if (use_range)
    range.ResetRanges(hist);

  return ret;
}

using result_map = std::unordered_map<std::string, std::vector<TH1D *>>;

result_map Proj3D(std::vector<RooUnfold *> &responses, ProjectionRange range,
                  bool use_range = true) {
  std::unordered_map<std::string, std::vector<TH1D *>> ret;
  for (auto h : responses) {
    RooUnfoldBayes *b = (RooUnfoldBayes *)h;

    if (use_range)
      range.SetRanges((TH3D *)h->Hreco());
    ret["x"].push_back((TH1D*) ((TH3D *)h->Hreco())->Project3D("x"));
    ret["x"].back()->SetName(
        dijetcore::MakeString(h->Hreco()->GetName(), "xproj",
                              dijetcore::Counter::instance().counter())
            .c_str());
    ret["x"].back()->Scale(1.0 / ret["x"].back()->Integral());
    ret["y"].push_back((TH1D*) ((TH3D *)h->Hreco())->Project3D("y"));
    ret["y"].back()->SetName(
        dijetcore::MakeString(h->Hreco()->GetName(), "yproj",
                              dijetcore::Counter::instance().counter())
            .c_str());
    ret["y"].back()->Scale(1.0 / ret["y"].back()->Integral());
    ret["z"].push_back((TH1D*) ((TH3D *)h->Hreco())->Project3D("z"));
    ret["z"].back()->SetName(
        dijetcore::MakeString(h->Hreco()->GetName(), "zproj",
                              dijetcore::Counter::instance().counter())
            .c_str());
    ret["z"].back()->Scale(1.0 / ret["z"].back()->Integral());

    if (use_range)
      range.ResetRanges((TH3D *)h->Hreco());
  }
  return ret;
}

std::vector<RooUnfold *> Unfold(RooUnfoldResponse *response,
                                RooUnfold::Algorithm alg, TH1 *measured,
                                std::vector<int> &iters) {
  std::vector<RooUnfold *> ret;

  for (auto &val : iters) {
    RooUnfoldBayes *b = new RooUnfoldBayes(response, measured, val, false);
    ret.push_back(b);
  }

  return ret;
}

bool matched_dijet(fastjet::PseudoJet &d1, fastjet::PseudoJet &d2,
                   fastjet::PseudoJet &t1, fastjet::PseudoJet &t2,
                   double match_rad) {
  fastjet::Selector matcher1 = fastjet::SelectorCircle(match_rad);
  matcher1.set_reference(d1);
  fastjet::Selector matcher2 = fastjet::SelectorCircle(match_rad);
  matcher2.set_reference(d2);

  std::vector<fastjet::PseudoJet> match1 =
      matcher1(std::vector<fastjet::PseudoJet>{t1});
  std::vector<fastjet::PseudoJet> match2 =
      matcher2(std::vector<fastjet::PseudoJet>{t2});
  return match1.size() == 1 && match2.size() == 1;
}

void FillResponseQG(ResponseSet &set, MatchEvent &e, bool match, bool g) {
  if (g) {
    if (match) {
      set.geom_match_gluon->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                                 Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()),
                                 e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                                 Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    } else {
      set.geom_match_gluon->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                                 Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
      set.geom_match_gluon->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                                 Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    }
    set.full_gluon->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                         Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()),
                         e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                         Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    set.full_no_miss_gluon->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                                 Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()),
                                 e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                                 Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  } else {
    if (match) {
      set.geom_match_quark->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                                 Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()),
                                 e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                                 Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    } else {
      set.geom_match_quark->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                                 Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
      set.geom_match_quark->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                                 Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    }
    set.full_quark->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                         Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()),
                         e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                         Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    set.full_no_miss_quark->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                                 Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()),
                                 e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                                 Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  }
}

void FillResponse(ResponseSet &set, MatchEvent &e, bool match) {
  if (match) {
    set.geom_match->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                         Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()),
                         e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                         Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  } else {
    set.geom_match->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                         Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    set.geom_match->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                         Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
  }
  set.full->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                 Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()), e.lead_hc_gen.pt(),
                 e.sub_hc_gen.pt(), Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  set.full_no_miss->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                         Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()),
                         e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                         Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  if (e.gen1_g && e.gen2_g) {
    FillResponseQG(set, e, match, true);
  } else if (!e.gen1_g && !e.gen2_g) {
    FillResponseQG(set, e, match, false);
  }
}

void FillResponseFake(ResponseSet &s, MatchEvent &e) {
  s.full->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
               Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
  s.full_no_miss->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                       Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
  s.geom_match->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                     Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
  if (e.gen1_g && e.gen2_g) {
    s.full_gluon->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                       Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    s.full_no_miss_gluon->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                               Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    s.geom_match_gluon->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                             Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
  } else if (!e.gen1_g && !e.gen2_g) {
    s.full_quark->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                       Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    s.full_no_miss_quark->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                               Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    s.geom_match_quark->Fake(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                             Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
  }
}

void FillResponseMiss(ResponseSet &s, MatchEvent &e) {
  s.full->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
               Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  s.geom_match->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                     Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  if (e.gen1_g && e.gen2_g) {
    s.full_gluon->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                       Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    s.geom_match_gluon->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                             Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  } else if (!e.gen1_g && !e.gen2_g) {
    s.full_quark->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                       Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    s.geom_match_quark->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                             Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  }
}

void FillResponseNoMiss(ResponseSet &s, MatchEvent &e) {
  s.full_no_miss->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                       Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  if (e.gen1_g && e.gen2_g) {
    s.full_no_miss_gluon->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                               Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  } else if (!e.gen1_g && !e.gen2_g) {
    s.full_no_miss_quark->Miss(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                               Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
  }
}

void FillHists(HistSet &s, MatchEvent &e, bool gen) {
  if (gen) {
    s.full->Fill(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                 Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    if (e.gen1_g && e.gen2_g) {
      s.gluon->Fill(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                    Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    } else if (!e.gen1_g && !e.gen2_g) {
      s.quark->Fill(e.lead_hc_gen.pt(), e.sub_hc_gen.pt(),
                    Aj(e.lead_hc_gen.pt(), e.sub_hc_gen.pt()));
    }
  } else {
    s.full->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                 Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    if (e.gen1_g && e.gen2_g) {
      s.gluon->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                    Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    } else if (!e.gen1_g && !e.gen2_g) {
      s.quark->Fill(e.lead_hc_det.pt(), e.sub_hc_det.pt(),
                    Aj(e.lead_hc_det.pt(), e.sub_hc_det.pt()));
    }
  }
}

int main(int argc, char *argv[]) {

  string usage = "Unfolding of Pythia Aj with STAR detector sim";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // parse configuration file
  dijetcore::json config;
  try {
    config = dijetcore::LoadConfig(FLAGS_config, required_params);
  } catch (std::exception &e) {
    LOG(ERROR) << "error loading config: " << e.what();
    return 1;
  }

  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  // find output directory from configuration file
  string output_dir = config["output_directory"];

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(output_dir);
  boost::filesystem::create_directories(dir);

  // copy config file to output directory
  boost::filesystem::path input_file(FLAGS_config.c_str());
  boost::filesystem::path copy_path(dir);
  copy_path /= input_file.filename();
  boost::filesystem::copy_file(
      input_file, copy_path,
      boost::filesystem::copy_option::overwrite_if_exists);

  // load input file and setup TTreeReader
  TFile input(FLAGS_input.c_str(), "READ");
  TTreeReader input_reader("resultTree", &input);
  TTreeReaderValue<TLorentzVector> lead_hc_gen_r(input_reader, "leadhcgen");
  TTreeReaderValue<TLorentzVector> sub_hc_gen_r(input_reader, "subhcgen");
  TTreeReaderValue<TLorentzVector> lead_match_gen_r(input_reader,
                                                    "leadmatchgen");
  TTreeReaderValue<TLorentzVector> sub_match_gen_r(input_reader, "submatchgen");
  TTreeReaderValue<TLorentzVector> lead_hc_det_r(input_reader, "leadhcdet");
  TTreeReaderValue<TLorentzVector> sub_hc_det_r(input_reader, "subhcdet");
  TTreeReaderValue<TLorentzVector> lead_match_det_r(input_reader,
                                                    "leadmatchdet");
  TTreeReaderValue<TLorentzVector> sub_match_det_r(input_reader, "submatchdet");
  TTreeReaderValue<TLorentzVector> gen_trig_r(input_reader, "triggen");
  TTreeReaderValue<TLorentzVector> det_trig_r(input_reader, "trigdet");

  TTreeReaderValue<TLorentzVector> init_part_1(input_reader, "gen1");
  TTreeReaderValue<TLorentzVector> init_part_2(input_reader, "gen2");

  TTreeReaderValue<bool> init_part_1_g(input_reader, "gen1_g");
  TTreeReaderValue<bool> init_part_2_g(input_reader, "gen2_g");

  // get pt cuts and match radius
  double match_rad = config["match_radius"];
  double lead_gen_min_pt = config["lead_jet_pt_gen"];
  double lead_min_pt = config["lead_jet_pt"];
  double sub_gen_min_pt = config["sublead_jet_pt_gen"];
  double sub_min_pt = config["sublead_jet_pt"];
  double gen_trig_threshold = config["gen_trig_threshold"];
  double det_trig_threshold = config["det_trig_threshold"];
  int num_events = config["num_events"];
  int seed = config["seed"];
  bool force_gen_miss_with_det_cuts = config["force_gen_miss_with_det_cuts"];
  std::vector<MatchEvent> matched_events;
  unsigned event_counter = 0;
  matched_events.reserve(10000);
  while (input_reader.Next()) {
    fastjet::PseudoJet lead_hc_gen(*lead_hc_gen_r);
    fastjet::PseudoJet sub_hc_gen(*sub_hc_gen_r);
    fastjet::PseudoJet lead_match_gen(*lead_match_gen_r);
    fastjet::PseudoJet sub_match_gen(*sub_match_gen_r);
    fastjet::PseudoJet lead_hc_det(*lead_hc_det_r);
    fastjet::PseudoJet sub_hc_det(*sub_hc_det_r);
    fastjet::PseudoJet lead_match_det(*lead_match_det_r);
    fastjet::PseudoJet sub_match_det(*sub_match_det_r);
    fastjet::PseudoJet gen_trig(*gen_trig_r);
    fastjet::PseudoJet det_trig(*det_trig_r);
    fastjet::PseudoJet gen1(*init_part_1);
    fastjet::PseudoJet gen2(*init_part_2);

    // check pT and trigger - if we don't need this event for either,
    // continue
    bool found_gen = lead_hc_gen.pt() > lead_gen_min_pt &&
                     sub_hc_gen.pt() > sub_gen_min_pt &&
                     gen_trig.pt() > gen_trig_threshold;
    bool found_det = lead_hc_det.pt() > lead_min_pt &&
                     sub_hc_det.pt() > sub_min_pt &&
                     det_trig.pt() > det_trig_threshold;
    bool found_gen_with_det_cuts =
        lead_hc_gen.pt() > lead_min_pt && sub_hc_gen.pt() > sub_min_pt;

    if (!found_det && !found_gen)
      continue;

    MatchEvent tmp;
    tmp.lead_hc_gen = lead_hc_gen;
    tmp.lead_match_gen = lead_match_gen;
    tmp.sub_hc_gen = sub_hc_gen;
    tmp.sub_match_gen = sub_match_gen;
    tmp.lead_hc_det = lead_hc_det;
    tmp.lead_match_det = lead_match_det;
    tmp.sub_hc_det = sub_hc_det;
    tmp.sub_match_det = sub_match_det;
    tmp.gen_trig = gen_trig;
    tmp.det_trig = det_trig;
    tmp.has_gen = found_gen;
    tmp.has_det = found_det;
    tmp.gen_with_det_cuts = found_gen_with_det_cuts;
    tmp.gen1 = gen1;
    tmp.gen2 = gen2;
    tmp.gen1_g = *init_part_1_g;
    tmp.gen2_g = *init_part_2_g;
    matched_events.push_back(tmp);

    if (num_events > 0 && ++event_counter >= num_events)
      break;
  }

  LOG(INFO) << "accepted events " << matched_events.size();

  // shuffle the events just in case
  std::random_device rd;
  std::mt19937 g;
  if (seed < 0)
    g.seed(rd());
  else
    g.seed(seed);
  std::shuffle(matched_events.begin(), matched_events.end(), g);

  // split training data into training and testing data
  std::vector<MatchEvent> testing(matched_events.begin(),
                                  matched_events.begin() +
                                      matched_events.size() / 4);
  std::vector<MatchEvent> training(
      matched_events.begin() + matched_events.size() / 4, matched_events.end());

  // build response and histograms
  HistBinning binning;
  TH3D *truth_dist_template =
      new TH3D("truthdisttemp", ";p_{T}^{lead};p_{T}^{sub};A_{J}",
               binning.lead_pt_bin, binning.lead_pt_min, binning.lead_pt_max,
               binning.sub_pt_bin, binning.sub_pt_min, binning.sub_pt_max,
               binning.aj_bin, binning.aj_min, binning.aj_max);
  TH3D *det_dist_template =
      new TH3D("detdisttemp", ";p_{T}^{lead};p_{T}^{sub};A_{J}",
               binning.lead_pt_bin, binning.lead_pt_min, binning.lead_pt_max,
               binning.sub_pt_bin, binning.sub_pt_min, binning.sub_pt_max,
               binning.aj_bin, binning.aj_min, binning.aj_max);
  ResponseSet response_set_hc =
      BuildResponseSet("hcresponse", truth_dist_template, det_dist_template);
  HistSet det_hist_set_hc = BuildHistSet("dethchists", binning);
  HistSet gen_hist_set_hc = BuildHistSet("genhchists", binning);

  int n_response_entries = 0;
  int n_response_flipped = 0;
  int total_entries = 0;
  int miss_entries = 0;

  // for QA, build TH1s for generator pT
  TH1D *pt_dist_quark = new TH1D("quarkpt", ";p_{T};counts", 50, 0, 50);
  TH1D *pt_dist_qg = new TH1D("qgpt", ";p_{T};counts", 50, 0, 50);
  for (auto &e : training) {
    // first, plot generator hc pt
    if (e.has_gen) {
      if (!e.gen1_g && !e.gen1_g) {
        pt_dist_quark->Fill(e.lead_hc_gen.pt());
        pt_dist_qg->Fill(e.lead_hc_gen.pt());
      } else if (e.gen1_g && e.gen1_g) {
        pt_dist_qg->Fill(e.lead_hc_gen.pt());
      }
    }

    // handle the two cases separately - when detector is present and when only
    // generator is present
    if (e.has_det) {
      // split this into the case where a generator level di-jet is present or
      // not
      if (e.has_gen) {
        total_entries++;
        n_response_entries++;
        bool match = matched_dijet(e.lead_hc_det, e.sub_hc_det, e.lead_hc_gen,
                                   e.sub_hc_gen, match_rad);
        FillResponse(response_set_hc, e, match);
        if (!match)
          n_response_flipped++;

        FillHists(gen_hist_set_hc, e, true);
      } else {
        FillResponseFake(response_set_hc, e);
      }
      FillHists(det_hist_set_hc, e, false);

    } else if (e.has_gen) {
      total_entries++;
      miss_entries++;
      FillResponseMiss(response_set_hc, e);
      if (force_gen_miss_with_det_cuts && e.gen_with_det_cuts) {
        FillResponseNoMiss(response_set_hc, e);
      }
      FillHists(gen_hist_set_hc, e, true);
    }
  }

  LOG(INFO) << "fraction of response entries with leading/subleading flip: "
            << (double)n_response_entries / n_response_flipped;
  LOG(INFO) << "fraction of entries that are misses: "
            << (double)miss_entries / total_entries;

  HistSet det_hist_set_hc_test = BuildHistSet("dethchiststest", binning);
  HistSet gen_hist_set_hc_test = BuildHistSet("genhchiststest", binning);

  for (auto &e : testing) {
    if (e.has_det) {
      FillHists(det_hist_set_hc_test, e, false);
    }
    if (e.has_gen) {
      FillHists(gen_hist_set_hc_test, e, true);
    }
  }

  // run unfolding on our two samples: train_hc, test_hc
  // and for quark/gluon distributions
  // test_match for {1, 2, 3, 4, 5} iters
  // std::vector<int> iters{1, 2, 3, 4, 5};
  std::vector<int> iters{4};

  std::vector<RooUnfold *> train_hc_unf = Unfold(
      response_set_hc.full, RooUnfold::kBayes, det_hist_set_hc.full, iters);
  std::vector<RooUnfold *> test_hc_unf =
      Unfold(response_set_hc.full, RooUnfold::kBayes, det_hist_set_hc_test.full,
             iters);
  std::vector<RooUnfold *> hc_unf_quark = Unfold(
      response_set_hc.full, RooUnfold::kBayes, det_hist_set_hc.quark, iters);
  std::vector<RooUnfold *> hc_unf_quark_no_miss =
      Unfold(response_set_hc.full_no_miss, RooUnfold::kBayes,
             det_hist_set_hc.quark, iters);
  std::vector<RooUnfold *> hc_unf_gluon = Unfold(
      response_set_hc.full, RooUnfold::kBayes, det_hist_set_hc.gluon, iters);
  std::vector<RooUnfold *> hc_unf_gluon_no_miss =
      Unfold(response_set_hc.full_no_miss, RooUnfold::kBayes,
             det_hist_set_hc.gluon, iters);
  std::vector<RooUnfold *> hc_unf_quark_g_res =
      Unfold(response_set_hc.full_gluon, RooUnfold::kBayes,
             det_hist_set_hc.quark, iters);
  std::vector<RooUnfold *> hc_unf_gluon_q_res =
      Unfold(response_set_hc.full_quark, RooUnfold::kBayes,
             det_hist_set_hc.gluon, iters);

  // now for response where lead jet gen/det are geometrically matched
  std::vector<RooUnfold *> train_hc_unf_geom_match =
      Unfold(response_set_hc.geom_match, RooUnfold::kBayes,
             det_hist_set_hc.full, iters);
  std::vector<RooUnfold *> test_hc_unf_geom_match =
      Unfold(response_set_hc.geom_match, RooUnfold::kBayes,
             det_hist_set_hc_test.full, iters);

  // build the projection ranges (one full range, one restricted range)
  ProjectionRange full_range;
  ProjectionRange user_range;
  user_range.leadptmin = config["proj_leadpt_min"];
  user_range.leadptmax = config["proj_leadpt_max"];
  user_range.subptmin = config["proj_subpt_min"];
  user_range.subptmax = config["proj_subpt_max"];
  user_range.ajmin = config["proj_aj_min"];
  user_range.ajmax = config["proj_aj_max"];


  result_map train_hc_unf_proj = Proj3D(train_hc_unf, user_range);
  result_map test_hc_unf_proj = Proj3D(test_hc_unf, user_range);
  result_map hc_unf_gluon_proj = Proj3D(hc_unf_gluon, user_range);
  result_map hc_unf_gluon_no_miss_proj = Proj3D(hc_unf_gluon_no_miss, user_range);
  result_map hc_unf_quark_proj = Proj3D(hc_unf_quark, user_range);
  result_map hc_unf_quark_no_miss_proj = Proj3D(hc_unf_quark_no_miss, user_range);
  result_map hc_unf_gluon_q_res_proj = Proj3D(hc_unf_gluon_q_res, user_range);
  result_map hc_unf_quark_g_res_proj = Proj3D(hc_unf_quark_g_res, user_range);
  result_map train_hc_unf_geom_match_proj = Proj3D(train_hc_unf_geom_match, user_range);
  result_map test_hc_unf_geom_match_proj = Proj3D(test_hc_unf_geom_match, user_range);

  // get truth level quantities
  hist_map train_hc_truth_proj = Proj3D(gen_hist_set_hc.full, user_range);
  hist_map test_hc_truth_proj = Proj3D(gen_hist_set_hc_test.full, user_range);

  // get detector level quantities
  hist_map train_hc_det_proj = Proj3D(det_hist_set_hc.full, user_range);
  hist_map test_hc_det_proj = Proj3D(det_hist_set_hc_test.full, user_range);

  // quark/gluon
  hist_map hc_truth_proj_gluon = Proj3D(gen_hist_set_hc.gluon, user_range);
  hist_map hc_truth_proj_quark = Proj3D(gen_hist_set_hc.quark, user_range);
  hist_map hc_det_proj_gluon = Proj3D(det_hist_set_hc.gluon, user_range);
  hist_map hc_det_proj_quark = Proj3D(det_hist_set_hc.quark, user_range);

  // setup printing utilities
  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;

  // first print QA plots
  pt_dist_quark->Divide(pt_dist_qg);
  dijetcore::canvasOpts coptsnoleg;
  coptsnoleg.do_legend = false;
  dijetcore::PrettyPrint1D(pt_dist_quark, hopts, coptsnoleg, "",
                           config["output_directory"], "quark_fraction", "",
                           "p_{T}", "quark fraction", "");

  // print Aj unfolding results
  // select the 4th iteration as our baseline (4th index)
  int print_iter = 0;

  // print response for hard-core and matched
  dijetcore::Print2DSimple(
      (TH2D *)train_hc_unf[print_iter]->response()->Hresponse()->Clone(), hopts,
      copts, config["output_directory"], "hc_response_no_match", "",
      "detector p_{T}^{lead}:p_{T}^{sub}:A_{J}",
      "generator p_{T}^{lead}:p_{T}^{sub}:A_{J}", "COLZ");
  dijetcore::Print2DSimple((TH2D *)train_hc_unf_geom_match[print_iter]
                               ->response()
                               ->Hresponse()
                               ->Clone(),
                           hopts, copts, config["output_directory"],
                           "hc_response_geom_match", "",
                           "detector p_{T}^{lead}:p_{T}^{sub}:A_{J}",
                           "generator p_{T}^{lead}:p_{T}^{sub}:A_{J}", "COLZ");

  // hard-core
  std::vector<TH1D *> aj_hc_test_results{test_hc_det_proj["z"],
                                         test_hc_unf_proj["z"][print_iter]};
  std::vector<std::string> aj_hc_test_results_labels{"detector", "unfolded"};
  test_hc_truth_proj["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(
      test_hc_truth_proj["z"], aj_hc_test_results, "generator",
      aj_hc_test_results_labels, hopts, copts, config["output_directory"],
      "hc_unfolded_test_no_match", "", "A_{J}", "event fraction");

  std::vector<TH1D *> aj_hc_train_results{train_hc_det_proj["z"],
                                          train_hc_unf_proj["z"][print_iter]};
  std::vector<std::string> aj_hc_train_results_labels{"detector", "unfolded"};
  train_hc_truth_proj["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(
      train_hc_truth_proj["z"], aj_hc_train_results, "generator",
      aj_hc_test_results_labels, hopts, copts, config["output_directory"],
      "hc_unfolded_train_no_match", "", "A_{J}", "event fraction");

  // hard-core with geometric matching
  std::vector<TH1D *> aj_hc_geom_match_test_results{
      test_hc_det_proj["z"], test_hc_unf_geom_match_proj["z"][print_iter]};
  std::vector<std::string> aj_hc_test_results_geom_match_labels{"detector",
                                                                "unfolded"};
  test_hc_truth_proj["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(
      test_hc_truth_proj["z"], aj_hc_geom_match_test_results, "generator",
      aj_hc_test_results_geom_match_labels, hopts, copts,
      config["output_directory"], "hc_unfolded_test_geom_match", "", "A_{J}",
      "event fraction");

  std::vector<TH1D *> aj_hc_geom_match_train_results{
      train_hc_det_proj["z"], train_hc_unf_geom_match_proj["z"][print_iter]};
  std::vector<std::string> aj_hc_train_results_geom_match_labels{"detector",
                                                                 "unfolded"};
  train_hc_truth_proj["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(
      train_hc_truth_proj["z"], aj_hc_geom_match_train_results, "generator",
      aj_hc_test_results_geom_match_labels, hopts, copts,
      config["output_directory"], "hc_unfolded_train_geom_match", "", "A_{J}",
      "event fraction");

  // overlay the unfolded distribution for matched and unmatched
  test_hc_unf_geom_match_proj["z"][print_iter]->GetYaxis()->SetRangeUser(0.0,
                                                                         0.4);
  dijetcore::PrintWithRatio(test_hc_unf_geom_match_proj["z"][print_iter],
                            test_hc_unf_proj["z"][print_iter],
                            "geometric match", "no match", hopts, copts,
                            config["output_directory"], "compare_matching", "",
                            "A_{J}", "event_fraction");

  // show quark and gluon differences
  // overlay the response for matched and unmatched
  std::vector<TH1D *> qg_comp_gen{hc_truth_proj_quark["z"],
                                  hc_truth_proj_gluon["z"]};
  std::vector<std::string> qg_comp_gen_labels{"quark", "gluon"};
  test_hc_truth_proj["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(test_hc_truth_proj["z"], qg_comp_gen, "q+g",
                             qg_comp_gen_labels, hopts, copts,
                             config["output_directory"], "quark_gluon_gen", "",
                             "A_{J}", "event_fraction");

  std::vector<TH1D *> qg_comp_det{hc_det_proj_quark["z"],
                                  hc_det_proj_gluon["z"]};
  std::vector<std::string> qg_comp_det_labels{"quark", "gluon"};
  train_hc_det_proj["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(train_hc_det_proj["z"], qg_comp_det, "q+g",
                             qg_comp_det_labels, hopts, copts,
                             config["output_directory"], "quark_gluon_det", "",
                             "A_{J}", "event_fraction");

  std::vector<TH1D *> qg_comp_unf{
      hc_unf_quark_proj["z"][print_iter],
      hc_unf_gluon_proj["z"][print_iter],
  };
  std::vector<std::string> qg_comp_unf_labels{"quark", "gluon"};
  test_hc_unf_proj["z"][print_iter]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(test_hc_unf_proj["z"][print_iter], qg_comp_unf,
                             "q+g", qg_comp_unf_labels, hopts, copts,
                             config["output_directory"], "quark_gluon_unf", "",
                             "A_{J}", "event_fraction");

  std::vector<TH1D *> quark_compare{hc_det_proj_quark["z"],
                                    hc_unf_quark_proj["z"][print_iter]};
  std::vector<std::string> quark_labels{"det quark", "unf quark"};
  hc_truth_proj_quark["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(hc_truth_proj_quark["z"], quark_compare,
                             "gen quark", quark_labels, hopts, copts,
                             config["output_directory"], "hc_quark_unf_compare",
                             "", "A_{J}", "event fraction");

  std::vector<TH1D *> gluon_compare{hc_det_proj_gluon["z"],
                                    hc_unf_gluon_proj["z"][print_iter]};
  std::vector<std::string> gluon_labels{"det gluon", "unf gluon"};
  hc_truth_proj_gluon["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(hc_truth_proj_gluon["z"], gluon_compare,
                             "gen gluon", gluon_labels, hopts, copts,
                             config["output_directory"], "hc_gluon_unf_compare",
                             "", "A_{J}", "event fraction");

  std::vector<TH1D *> quark_g_res_compare{
      hc_det_proj_quark["z"], hc_unf_quark_g_res_proj["z"][print_iter]};
  std::vector<std::string> quark_g_res_labels{"det quark", "unf quark"};
  hc_truth_proj_quark["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(
      hc_truth_proj_quark["z"], quark_g_res_compare, "gen quark",
      quark_g_res_labels, hopts, copts, config["output_directory"],
      "hc_quark_g_res_unf_compare", "", "A_{J}", "event fraction");

  std::vector<TH1D *> gluon_q_res_compare{
      hc_det_proj_gluon["z"], hc_unf_gluon_q_res_proj["z"][print_iter]};
  std::vector<std::string> gluon_q_res_labels{"det gluon", "unf gluon"};
  hc_truth_proj_gluon["z"]->GetYaxis()->SetRangeUser(0.0, 0.4);
  dijetcore::PrintWithRatios(
      hc_truth_proj_gluon["z"], gluon_q_res_compare, "gen gluon",
      gluon_q_res_labels, hopts, copts, config["output_directory"],
      "hc_gluon_g_res_unf_compare", "", "A_{J}", "event fraction");

  // compare the unfolding result with and without misses for gluon
  std::vector<TH1D *> gluon_miss_compare{
      hc_unf_gluon_proj["z"][print_iter],
      hc_unf_gluon_no_miss_proj["z"][print_iter]};
  std::vector<std::string> gluon_miss_compare_labels{"full response",
                                                     "no misses"};
  dijetcore::PrintWithRatios(
      hc_truth_proj_gluon["z"], gluon_miss_compare, "gen gluon",
      gluon_miss_compare_labels, hopts, copts, config["output_directory"],
      "hc_gluon_unf_compare_no_miss", "", "A_{J}", "event fraction");

  // auau & pp comparison if it is present
  if (!config["pp_file"].empty() && !config["auau_file"].empty()) {
    std::string pp_filename = config["pp_file"];
    std::string auau_filename = config["auau_file"];
    std::string treename = config["tree_name"];

    TFile pp_file(pp_filename.c_str(), "READ");
    TTreeReader pp_reader(treename.c_str(), &pp_file);
    TTreeReaderValue<TLorentzVector> pp_lead_hc(pp_reader, "jl");
    TTreeReaderValue<TLorentzVector> pp_sub_hc(pp_reader, "js");

    TH1D *pp_aj = new TH1D("pp_aj", ";A_{J};event fraction", 10, 0, 1.0);
    TH3D *pp_det =
        new TH3D("pp_det", ";p_{T}^{lead};p_{T}^{sub};A_{J}",
                 binning.lead_pt_bin, binning.lead_pt_min, binning.lead_pt_max,
                 binning.sub_pt_bin, binning.sub_pt_min, binning.sub_pt_max,
                 binning.aj_bin, binning.aj_min, binning.aj_max);
    while (pp_reader.Next()) {
      pp_aj->Fill(Aj((*pp_lead_hc).Pt(), (*pp_sub_hc).Pt()));
      pp_det->Fill((*pp_lead_hc).Pt(), (*pp_sub_hc).Pt(),
                   Aj((*pp_lead_hc).Pt(), (*pp_sub_hc).Pt()));
    }

    TFile auau_file(auau_filename.c_str(), "READ");
    TTreeReader auau_reader(treename.c_str(), &auau_file);
    TTreeReaderValue<int> auau_cent(auau_reader, "cent");
    TTreeReaderValue<TLorentzVector> auau_lead_hc(auau_reader, "jl");
    TTreeReaderValue<TLorentzVector> auau_sub_hc(auau_reader, "js");

    TH1D *auau_aj = new TH1D("auau_aj", ";A_{J};event fraction", 10, 0, 1.0);
    TH3D *auau_det =
        new TH3D("auau_det", ";p_{T}^{lead};p_{T}^{sub};A_{J}",
                 binning.lead_pt_bin, binning.lead_pt_min, binning.lead_pt_max,
                 binning.sub_pt_bin, binning.sub_pt_min, binning.sub_pt_max,
                 binning.aj_bin, binning.aj_min, binning.aj_max);
    while (auau_reader.Next()) {
      if ((*auau_cent) >= 0 && (*auau_cent) <= 2)
        auau_aj->Fill(Aj((*auau_lead_hc).Pt(), (*auau_sub_hc).Pt()));
      auau_det->Fill((*auau_lead_hc).Pt(), (*auau_sub_hc).Pt(),
                     Aj((*auau_lead_hc).Pt(), (*auau_sub_hc).Pt()));
    }

    pp_aj->Scale(1.0 / pp_aj->Integral());
    auau_aj->Scale(1.0 / auau_aj->Integral());
    pp_aj->GetYaxis()->SetRangeUser(0.0, 0.4);
    auau_aj->GetYaxis()->SetRangeUser(0.0, 0.4);
    std::vector<TH1D *> auau_pp_comp{auau_aj, pp_aj};
    std::vector<std::string> auau_pp_comp_labels{"AuAu", "pp"};
    dijetcore::PrintWithRatios(
        train_hc_det_proj["z"], auau_pp_comp, "q+g detector",
        auau_pp_comp_labels, hopts, copts, config["output_directory"],
        "auau_pp_comp_det", "", "A_{J}", "event_fraction");

    std::vector<RooUnfold *> pp_unfold =
        Unfold(response_set_hc.full, RooUnfold::kBayes, pp_det, iters);
    std::vector<RooUnfold *> auau_unfold =
        Unfold(response_set_hc.full, RooUnfold::kBayes, auau_det, iters);

    result_map pp_unf_proj = Proj3D(pp_unfold, user_range);
    result_map auau_unf_proj = Proj3D(auau_unfold, user_range);

    std::vector<TH1D *> auau_pp_unf_comp{auau_unf_proj["z"][print_iter],
                                         pp_unf_proj["z"][print_iter]};
    std::vector<std::string> auau_pp_unf_comp_labels{"AuAu", "pp"};
    dijetcore::PrintWithRatios(
        train_hc_truth_proj["z"], auau_pp_unf_comp, "q+g generator",
        auau_pp_unf_comp_labels, hopts, copts, config["output_directory"],
        "auau_pp_comp_unf", "", "A_{J}", "event_fraction");
  }

  gflags::ShutDownCommandLineFlags();
  return 0;
}
