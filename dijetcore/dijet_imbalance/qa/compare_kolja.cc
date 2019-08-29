// compares my 20-10 aj trees to Koljas to find differences in jet selection

#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/json.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/map.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/types.h"
#include "dijetcore/util/data/centrality/centrality_run7.h"
#include "dijetcore/util/data/reader_util.h"
#include "dijetcore/util/data/trigger_lookup.h"
#include "dijetcore/util/data/vector_conversion.h"
#include "dijetcore/util/root/root_print_utils.h"

#include <set>

#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "sct/lib/root_ext/TTreeReader.h"
#include "TTreeReaderValue.h"

// define a structure to store event-by-event results
struct event_info {
  unsigned refmult = 0;
  double rho = 0.0;
  TLorentzVector jl;
  TLorentzVector js;
  TLorentzVector jlm;
  TLorentzVector jsm;
};

double hardAj(event_info &info) {
  return (info.jl.Pt() - info.js.Pt()) / (info.jl.Pt() + info.js.Pt());
}

double matchAj(event_info &info) {
  return (info.jlm.Pt() - info.jsm.Pt()) / (info.jlm.Pt() + info.jsm.Pt());
}

using dijetcore::string;

DIJETCORE_DEFINE_string(kolja_input, "", "kolja's aj file");
DIJETCORE_DEFINE_string(my_input, "", "my aj file");
DIJETCORE_DEFINE_string(outputDir, "tmp", "directory for output");

int main(int argc, char *argv[]) {
  string usage = "compares my 20_10 0-20% Aj results to Kolja's event-by-event";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outputDir);
  boost::filesystem::create_directories(dir);

  // set drawing preferences for histograms and graphs
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  TFile kolja_file(FLAGS_kolja_input.c_str(), "READ");
  TFile input_file(FLAGS_my_input.c_str(), "READ");

  // initialize a reader for both of our trees - we both use different
  // structures, so the tree reader setup is different
  TTreeReader kolja("ResultTree", &kolja_file);
  TTreeReaderValue<unsigned> k_run(kolja, "runid");
  TTreeReaderValue<unsigned> k_event(kolja, "eventid");
  TTreeReaderValue<double> k_refmult(kolja, "refmult");
  TTreeReaderValue<TLorentzVector> k_jl(kolja, "j1");
  TTreeReaderValue<TLorentzVector> k_js(kolja, "j2");
  TTreeReaderValue<TLorentzVector> k_jlm(kolja, "jm1");
  TTreeReaderValue<TLorentzVector> k_jsm(kolja, "jm2");
  TTreeReaderValue<float> k_rho(kolja, "rho");

  std::string dijet_key =
      "LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_"
      "0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_"
      "pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2;1LEAD_INIT_R_0.4_"
      "alg_2_pt_20_const_eta_1_const_pt_2_MATCH_R_0.4_alg_2_pt_0_const_eta_1_"
      "const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_pt_2_MATCH_R_"
      "0.4_alg_2_pt_0_const_eta_1_const_pt_0.2";
  TTreeReader mytree(dijet_key.c_str(), &input_file);
  TTreeReaderValue<int> m_run(mytree, "runid");
  TTreeReaderValue<int> m_event(mytree, "eventid");
  TTreeReaderValue<int> m_refmult(mytree, "grefmult");
  TTreeReaderValue<TLorentzVector> m_jl(mytree, "jl");
  TTreeReaderValue<TLorentzVector> m_js(mytree, "js");
  TTreeReaderValue<TLorentzVector> m_jlm(mytree, "jlm");
  TTreeReaderValue<TLorentzVector> m_jsm(mytree, "jsm");
  TTreeReaderValue<double> m_rho(mytree, "jlmrho");

  // define a map to hold results
  using result_map = dijetcore::dijetcore_map<std::pair<unsigned, unsigned>,
                                              event_info, dijetcore::PairHash>;

  result_map kolja_events;
  result_map my_events;

  std::set<std::pair<unsigned, unsigned>> all_event_keys;
  while (kolja.Next()) {
    std::pair<unsigned, unsigned> key;
    key.first = *k_run;
    key.second = *k_event;
    event_info tmp;
    tmp.refmult = *k_refmult;
    tmp.rho = *k_rho;
    tmp.jl = *k_jl;
    tmp.js = *k_js;
    tmp.jlm = *k_jlm;
    tmp.jsm = *k_jsm;
    if (tmp.refmult < 269)
      continue;
    all_event_keys.insert(key);
    kolja_events[key] = tmp;
  }

  while (mytree.Next()) {
    std::pair<unsigned, unsigned> key;
    key.first = *m_run;
    key.second = *m_event;
    event_info tmp;
    tmp.refmult = *m_refmult;
    tmp.rho = *m_rho;
    tmp.jl = *m_jl;
    tmp.js = *m_js;
    tmp.jlm = *m_jlm;
    tmp.jsm = *m_jsm;
    if (tmp.refmult < 269)
      continue;
    all_event_keys.insert(key);
    my_events[key] = tmp;
  }

  // calculate all unmatched quantities
  TH1D *unmatched_rho = new TH1D("unmatch_rho", ";#rho", 50, 30, 100);
  TH1D *unmatched_jlpt = new TH1D("unmatchedleadpt", ";p_{T}", 50, 0, 50);
  TH1D *unmatched_jspt = new TH1D("unmatchedsubpt", ";p_{T}", 50, 0, 50);
  TH1D *unmatched_jlmpt = new TH1D("unmatchedleadmatchpt", ";p_{T}", 50, 0, 50);
  TH1D *unmatched_jsmpt = new TH1D("unmatchedsubmatchpt", ";p_{T}", 50, 0, 50);
  TH1D *unmatched_ref = new TH1D("unmatched_refmult", ";remfult", 100, 0, 600);
  TH1D *unmatched_jleta = new TH1D("unmatchedleadeta", "#eta", 50, -1.0, 1.0);
  TH1D *unmatched_jseta = new TH1D("unmatchedsubeta", "#eta", 50, -1.0, 1.0);
  TH1D *unmatched_jlphi =
      new TH1D("unmatchedleadphi", "#phi", 50, -TMath::Pi(), TMath::Pi());
  TH1D *unmatched_jsphi =
      new TH1D("unmatchedsubphi", "#phi", 50, -TMath::Pi(), TMath::Pi());
  TH1D *unmatched_aj_hard = new TH1D("unmatchedajhard", "A_{J}", 20, -0.3, 0.9);
  TH1D *unmatched_aj_matched =
      new TH1D("unmatchedajmatched", "A_{J}", 20, -0.3, 0.9);

  for (auto &key : all_event_keys) {
    if (!my_events.count(key)) {
      LOG(INFO) << "kolja has event that I dont: " << key;
      LOG(INFO) << "refmult: " << kolja_events[key].refmult;
      unmatched_rho->Fill(kolja_events[key].rho);
      unmatched_jlpt->Fill(kolja_events[key].jl.Pt());
      unmatched_jspt->Fill(kolja_events[key].js.Pt());
      unmatched_jlmpt->Fill(kolja_events[key].jlm.Pt());
      unmatched_jsmpt->Fill(kolja_events[key].jsm.Pt());
      unmatched_ref->Fill(kolja_events[key].refmult);
      unmatched_jleta->Fill(kolja_events[key].jl.Eta());
      unmatched_jseta->Fill(kolja_events[key].js.Eta());
      unmatched_jlphi->Fill(kolja_events[key].jl.Phi());
      unmatched_jsphi->Fill(kolja_events[key].js.Phi());
      unmatched_aj_hard->Fill(hardAj(kolja_events[key]));
      unmatched_aj_matched->Fill(matchAj(kolja_events[key]));
    } else if (!kolja_events.count(key)) {
      LOG(INFO) << "I have event that Kolja doesnt: " << key;
      LOG(INFO) << "refmult: " << my_events[key].refmult;
      unmatched_rho->Fill(my_events[key].rho);
      unmatched_jlpt->Fill(my_events[key].jl.Pt());
      unmatched_jspt->Fill(my_events[key].js.Pt());
      unmatched_jlmpt->Fill(my_events[key].jlm.Pt());
      unmatched_jsmpt->Fill(my_events[key].jsm.Pt());
      unmatched_ref->Fill(my_events[key].refmult);
      unmatched_jleta->Fill(my_events[key].jl.Eta());
      unmatched_jseta->Fill(my_events[key].js.Eta());
      unmatched_jlphi->Fill(my_events[key].jl.Phi());
      unmatched_jsphi->Fill(my_events[key].js.Phi());
      unmatched_aj_hard->Fill(hardAj(my_events[key]));
      unmatched_aj_matched->Fill(matchAj(my_events[key]));
    }
  }
  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;

  dijetcore::PrettyPrint1D(unmatched_rho, hopts, copts, "unmatched rho",
                           FLAGS_outputDir, "unmatchedrho", "", "#rho",
                           "counts");
  dijetcore::PrettyPrint1D(unmatched_jlpt, hopts, copts, "unmatched lead pT",
                           FLAGS_outputDir, "unmatchedjlpt", "", "p_{T}",
                           "counts");
  dijetcore::PrettyPrint1D(unmatched_jlmpt, hopts, copts,
                           "unmatched lead match pT", FLAGS_outputDir,
                           "unmatchedjlmpt", "", "p_{T}", "counts");
  dijetcore::PrettyPrint1D(unmatched_jspt, hopts, copts, "unmatched sub pT",
                           FLAGS_outputDir, "unmatchedjspt", "", "p_{T}",
                           "counts");
  dijetcore::PrettyPrint1D(unmatched_jsmpt, hopts, copts,
                           "unmatched sub match pT", FLAGS_outputDir,
                           "unmatchedjsmpt", "", "p_{T}", "counts");
  dijetcore::PrettyPrint1D(unmatched_jleta, hopts, copts, "unmatched lead eta",
                           FLAGS_outputDir, "unmatchedjleta", "", "#eta",
                           "counts");
  dijetcore::PrettyPrint1D(unmatched_jseta, hopts, copts,
                           "unmatched sub lead eta", FLAGS_outputDir,
                           "unmatchedjseta", "", "#eta", "counts");
  dijetcore::PrettyPrint1D(unmatched_jlphi, hopts, copts, "unmatched lead phi",
                           FLAGS_outputDir, "unmatchedjlphi", "", "#phi",
                           "counts");
  dijetcore::PrettyPrint1D(unmatched_jsphi, hopts, copts,
                           "unmatched sub lead phi", FLAGS_outputDir,
                           "unmatchedjsphi", "", "#phi", "counts");
  dijetcore::PrettyPrint1D(unmatched_aj_hard, hopts, copts, "unmatched hard aj",
                           FLAGS_outputDir, "unmatchedhardaj", "", "A_{J}",
                           "counts");
  dijetcore::PrettyPrint1D(unmatched_aj_matched, hopts, copts,
                           "unmatched matched aj", FLAGS_outputDir,
                           "unmatchedmatchedaj", "", "A_{J}", "counts");

  // plot all matched quantities and see if there is any large difference
  TH2D *jlpt = new TH2D("jlpt", ";p_{T}", 50, 0, 50, 50, 0, 50);
  TH2D *jspt = new TH2D("jspt", ";p_{T}", 50, 0, 50, 50, 0, 50);
  TH2D *jlmpt = new TH2D("jlmpt", ";p_{T}", 50, 0, 50, 50, 0, 50);
  TH2D *jsmpt = new TH2D("jsmpt", ";p_{T}", 50, 0, 50, 50, 0, 50);
  TH1D *diff_pt_hard = new TH1D("diffptlead", ";p_{T}", 50, -0.1, 0.1);
  TH1D *diff_pt_match = new TH1D("diffptleadmatch", ";p_{T}", 50, -0.5, 0.5);
  TH1D *diff_pt_hard_sub = new TH1D("diffptsub", ";p_{T}", 50, -0.1, 0.1);
  TH1D *diff_pt_match_sub = new TH1D("diffptsubmatch", ";p_{T}", 50, -0.5, 0.5);
  TH1D *diff_rho = new TH1D("diffrho", ";d#rho", 50, -30, 30);
  TH1D *diff_rho_norm = new TH1D("diffrhonorm", ";d#rho / #rho", 50, -0.3, 0.3);
  TH2D *mismatched_jets = new TH2D("mismatch", ";#eta;#phi", 50, -1, 1, 50,
                                   -TMath::Pi(), TMath::Pi());

  for (auto &key : all_event_keys) {
    if (my_events.count(key) && kolja_events.count(key)) {
      if (my_events[key].refmult != kolja_events[key].refmult) {
        DIJETCORE_THROW("shouldn't be here ", my_events[key].refmult, " ",
                        kolja_events[key].refmult);
      }
      if (my_events[key].refmult < 269)
        continue;

      double my_hard_aj = hardAj(my_events[key]);
      double my_match_aj = matchAj(my_events[key]);
      double k_hard_aj = hardAj(kolja_events[key]);
      double k_match_aj = matchAj(kolja_events[key]);

      jlpt->Fill(my_events[key].jl.Pt(), kolja_events[key].jl.Pt());
      jspt->Fill(my_events[key].js.Pt(), kolja_events[key].js.Pt());
      jlmpt->Fill(my_events[key].jlm.Pt(), kolja_events[key].jlm.Pt());
      jsmpt->Fill(my_events[key].jsm.Pt(), kolja_events[key].jsm.Pt());

      diff_pt_hard->Fill((my_events[key].jl.Pt() - kolja_events[key].jl.Pt()) /
                         my_events[key].jl.Pt());
      diff_pt_match->Fill(
          (my_events[key].jlm.Pt() - kolja_events[key].jlm.Pt()) /
          my_events[key].jlm.Pt());

      diff_pt_hard_sub->Fill(
          (my_events[key].js.Pt() - kolja_events[key].js.Pt()) /
          my_events[key].js.Pt());
      diff_pt_match_sub->Fill(
          (my_events[key].jsm.Pt() - kolja_events[key].jsm.Pt()) /
          my_events[key].jsm.Pt());

      diff_rho->Fill((my_events[key].rho - kolja_events[key].rho));
      diff_rho_norm->Fill((my_events[key].rho - kolja_events[key].rho) /
                          my_events[key].rho);

      if (fabs(my_events[key].jl.Pt() - kolja_events[key].jl.Pt()) /
              my_events[key].jl.Pt() >
          0.05) {
        // LOG(INFO) << "found an event: " << key;
        // LOG(INFO) << "my pt: " << my_events[key].jl.Pt();
        // LOG(INFO) << "kolja pt: " << kolja_events[key].jl.Pt();
        // LOG(INFO) << "my eta: " << my_events[key].jl.Eta();
        // LOG(INFO) << "Kolja eta: " << kolja_events[key].jl.Eta();
        // LOG(INFO) << "my phi: " << my_events[key].jl.Phi();
        // LOG(INFO) << "Kolja phi: " << kolja_events[key].jl.Phi();
        mismatched_jets->Fill(my_events[key].jl.Eta(), my_events[key].jl.Phi());
      }
    }
  }

  dijetcore::Print2DSimple(jlpt, hopts, copts, FLAGS_outputDir, "leadhardpt2d",
                           "", "my jets p_{T}", "koljas jets p_{T}");
  dijetcore::Print2DSimple(jspt, hopts, copts, FLAGS_outputDir, "subhardpt2d",
                           "", "my jets p_{T}", "koljas jets p_{T}");
  dijetcore::Print2DSimple(jlmpt, hopts, copts, FLAGS_outputDir,
                           "leadmatchpt2d", "", "my jets p_{T}",
                           "koljas jets p_{T}");
  dijetcore::Print2DSimple(jsmpt, hopts, copts, FLAGS_outputDir, "submatchpt2d",
                           "", "my jets p_{T}", "koljas jets p_{T}");
  dijetcore::Print2DSimple(mismatched_jets, hopts, copts, FLAGS_outputDir,
                           "mismatched", "", "dp_{T}/p_{T}", "#phi");

  dijetcore::PrettyPrint1D(diff_pt_hard, hopts, copts, "diff hard pt lead",
                           FLAGS_outputDir, "diffhardptlead", "", "A_{J}",
                           "counts");
  dijetcore::PrettyPrint1D(diff_pt_match, hopts, copts, "diff hard pt match",
                           FLAGS_outputDir, "diffmatchptlead", "",
                           "dp_{T}/p_{T}", "counts");

  dijetcore::PrettyPrint1D(diff_pt_hard_sub, hopts, copts, "diff hard pt sub",
                           FLAGS_outputDir, "diffhardptsub", "", "A_{J}",
                           "counts");
  dijetcore::PrettyPrint1D(diff_pt_match_sub, hopts, copts,
                           "diff hard pt match sub", FLAGS_outputDir,
                           "diffmatchptsub", "", "dp_{T}/p_{T}", "counts");
  dijetcore::PrettyPrint1D(diff_rho, hopts, copts, "diff rho", FLAGS_outputDir,
                           "diffrho", "", "d#rho", "counts");
  dijetcore::PrettyPrint1D(diff_rho_norm, hopts, copts, "diff rho norm",
                           FLAGS_outputDir, "diffrhonorm", "", "d#rho / #rho",
                           "counts");

  LOG(INFO) << "mean relative matched pT: " << diff_pt_match->GetMean();
  LOG(INFO) << "mean relative matched pT sublead: "
            << diff_pt_match_sub->GetMean();
  LOG(INFO) << "mean delta rho: " << diff_rho->GetMean();
  LOG(INFO) << "mean relative rho: " << diff_rho_norm->GetMean();

  gflags::ShutDownCommandLineFlags();
}
