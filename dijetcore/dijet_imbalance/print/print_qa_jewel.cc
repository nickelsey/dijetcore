#include "dijetcore/lib/flags.h"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/util/fastjet/dijet_key.h"
#include "dijetcore/util/root/root_print_utils.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector2.h"
#include "TH2D.h"

#include "dijetcore/lib/map.h"
#include <string>

// include boost for filesystem manipulation
#include "boost/filesystem.hpp"

using std::string;

DIJETCORE_DEFINE_string(input, "",
                        "input root file with processed jewel di-jet trees");
DIJETCORE_DEFINE_string(outputDir, "tmp", "directory for output");
DIJETCORE_DEFINE_double(ptrigger, 5.4, "hadron trigger ET threshold");
DIJETCORE_DEFINE_double(jtrigger, 10, "jet trigger ET threshold");

int main(int argc, char *argv[]) {
  string usage = "JEWEL differential di-jet imbalance print routine";

  dijetcore::SetUsageMessage(usage);
  dijetcore::ParseCommandLineFlags(&argc, argv);
  dijetcore::InitLogging(argv[0]);

  // set drawing preferences for histograms and graphs
  // gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  // turn off print messages
  gErrorIgnoreLevel = kInfo + 1;

  // check to make sure we have valid inputs
  if (!boost::filesystem::exists(FLAGS_input)) {
    LOG(ERROR) << "input file " << FLAGS_input << "doesn't exist: exiting"
               << std::endl;
    return 1;
  }

  // build output directory if it doesn't exist, using boost::filesystem
  if (FLAGS_outputDir.empty())
    FLAGS_outputDir = "tmp";
  boost::filesystem::path dir(FLAGS_outputDir.c_str());
  boost::filesystem::create_directories(dir);

  // read in the file
  TFile input_file(FLAGS_input.c_str(), "READ");

  // save all histograms so we can do comparisons
  // between different keys if we want
  TH2D *jt_vtx = new TH2D("jtvtx", "", 80, -8, 8, 80, -8, 8);
  TH2D *pt_vtx = new TH2D("ptvtx", "", 80, -8, 8, 80, -8, 8);
  TH1D* p1_spectrum_jet = new TH1D("p1jet", "", 100, 0, 100);
  TH1D* p1_spectrum_had = new TH1D("p1had", "", 100, 0, 100);
  // create TTreeReader
  TTree *t = (TTree *)input_file.Get("T");
  TTreeReader reader(t);
  TTreeReaderValue<TLorentzVector> tj(reader, "tj");
  TTreeReaderValue<TLorentzVector> tp(reader, "tp");
  TTreeReaderValue<TLorentzVector> p1(reader, "p1");
  TTreeReaderValue<TLorentzVector> p2(reader, "p2");
  TTreeReaderValue<double> w(reader, "w");
  TTreeReaderValue<double> xsec(reader, "xsec");
  TTreeReaderValue<double> vx(reader, "vx");
  TTreeReaderValue<double> vy(reader, "vy");

  size_t total_events = reader.GetTree()->GetEntries();
  LOG(INFO) << "number of events: " << total_events;
  size_t event_counter = 0;
  while (reader.Next()) {
    if (event_counter % (total_events / 10) == 0)
      LOG(INFO) << "event: " << event_counter;
    event_counter++;

    if (*xsec > 4000000)
      continue;
    // if (*w > 0.01)
    //   continue;
    // to fill the rotated entries, we first have to find the angle of
    // rotation
    TVector2 event_vector(*vx, *vy);
    TVector2 jet_vector((*tj).Px(), (*tj).Py());
    TVector2 hadron_vector((*tp).Px(), (*tp).Py());
    TVector2 event_vector_rotated_jet = event_vector.Rotate(TMath::Pi() - jet_vector.Phi());
    TVector2 event_vector_rotated_hadron = event_vector.Rotate(TMath::Pi() - hadron_vector.Phi());

    if ((*tj).Pt() > FLAGS_jtrigger) {
      jt_vtx->Fill(event_vector_rotated_jet.X(), event_vector_rotated_jet.Y(), *xsec);
      p1_spectrum_jet->Fill(p1->Pt(), *xsec);
    }
    if ((*tp).Pt() > FLAGS_ptrigger) {
      pt_vtx->Fill(event_vector_rotated_hadron.X(), event_vector_rotated_hadron.Y(), *xsec);
      p1_spectrum_had->Fill(p1->Pt(), *xsec);
    }
    

  }

  // create our histogram and canvas options
  dijetcore::histogramOpts hopts;
  dijetcore::canvasOpts copts;
  copts.leg_left_bound = 0.6;
  dijetcore::canvasOpts coptslogz;
  coptslogz.leg_left_bound = 0.6;
  coptslogz.log_z = true;
  dijetcore::canvasOpts coptslogy;
  coptslogy.leg_left_bound = 0.6;
  coptslogy.log_y = true;
  dijetcore::canvasOpts cOptsBottomLeg;
  cOptsBottomLeg.leg_upper_bound = 0.4;
  cOptsBottomLeg.leg_lower_bound = 0.18;
  cOptsBottomLeg.leg_right_bound = 0.9;
  cOptsBottomLeg.leg_left_bound = 0.7;
  dijetcore::canvasOpts cOptsBottomLeftLeg;
  cOptsBottomLeftLeg.leg_upper_bound = 0.4;
  cOptsBottomLeftLeg.leg_lower_bound = 0.18;
  cOptsBottomLeftLeg.leg_right_bound = 0.18;
  cOptsBottomLeftLeg.leg_left_bound = 0.4;
  dijetcore::canvasOpts cOptsBottomLeftLegLogy;
  cOptsBottomLeftLegLogy.log_y = true;
  cOptsBottomLeftLegLogy.leg_upper_bound = 0.4;
  cOptsBottomLeftLegLogy.leg_lower_bound = 0.18;
  cOptsBottomLeftLegLogy.leg_right_bound = 0.18;
  cOptsBottomLeftLegLogy.leg_left_bound = 0.4;
  dijetcore::canvasOpts cOptsTopLeftLeg;
  cOptsTopLeftLeg.leg_right_bound = 0.18;
  cOptsTopLeftLeg.leg_left_bound = 0.4;

  //dijetcore::Print2DSimple
  std::string jet_string = dijetcore::MakeString("jet trigger E_{T} > ", FLAGS_jtrigger, " GeV");
  std::string hadron_string = dijetcore::MakeString("hadron trigger E_{T} > ", FLAGS_ptrigger, " GeV").c_str();

  dijetcore::Print2DSimple(jt_vtx, hopts, copts, FLAGS_outputDir, "jet_trig_2d", jet_string, "x [fm]", "y [fm]", "CONTZ");
  dijetcore::Print2DSimple(pt_vtx, hopts, copts, FLAGS_outputDir, "had_trig_2d", hadron_string, "x [fm]", "y [fm]", "CONTZ");

  TFile* output = new TFile(dijetcore::MakeString(FLAGS_outputDir, "out.root").c_str(), "RECREATE");
  output->cd();
  jt_vtx->Write();
  pt_vtx->Write();
  p1_spectrum_jet->Write();
  p1_spectrum_had->Write();
  output->Close();
  
  return 0;
}

// void Print2DSimple(H* h,
//                      histogramOpts hopts,
//                      canvasOpts copts,
//                      std::string output_loc,
//                      std::string output_name,
//                      std::string canvas_title,
//                      std::string x_axis_label,
//                      std::string y_axis_label,
//                      std::string opt = "COLZ")