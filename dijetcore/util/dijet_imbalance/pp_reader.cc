#define PPReader_cxx
#include "PPReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void PPReader::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PPReader.C
//      root> PPReader t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}


PPReader::PPReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("trees/y6/y6_grid_16_8/tow_0_track_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("trees/y6/y6_grid_16_8/tow_0_track_0.root");
      }
      f->GetObject("LEAD_INIT_R_0.4_alg_2_pt_16_const_eta_1_const_pt_1_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_8_const_eta_1_const_pt_1_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2",tree);

   }
   Init(tree);
}

PPReader::~PPReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PPReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PPReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PPReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jl = 0;
   js = 0;
   jlm = 0;
   jsm = 0;
   ppjl = 0;
   ppjs = 0;
   ppjlm = 0;
   ppjsm = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runid", &runid, &b_runid);
   fChain->SetBranchAddress("eventid", &eventid, &b_eventid);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("refmult", &refmult, &b_refmult);
   fChain->SetBranchAddress("grefmult", &grefmult, &b_grefmult);
   fChain->SetBranchAddress("refmultcorr", &refmultcorr, &b_refmultcorr);
   fChain->SetBranchAddress("grefmultcorr", &grefmultcorr, &b_grefmultcorr);
   fChain->SetBranchAddress("cent", &cent, &b_cent);
   fChain->SetBranchAddress("zdcrate", &zdcrate, &b_zdcrate);
   fChain->SetBranchAddress("rp", &rp, &b_rp);
   fChain->SetBranchAddress("nglobal", &nglobal, &b_nglobal);
   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("jl", &jl, &b_jl);
   fChain->SetBranchAddress("js", &js, &b_js);
   fChain->SetBranchAddress("jlm", &jlm, &b_jlm);
   fChain->SetBranchAddress("jsm", &jsm, &b_jsm);
   fChain->SetBranchAddress("jlconst", &jlconst, &b_jlconst);
   fChain->SetBranchAddress("jlrho", &jlrho, &b_jlrho);
   fChain->SetBranchAddress("jlsig", &jlsig, &b_jlsig);
   fChain->SetBranchAddress("jlarea", &jlarea, &b_jlarea);
   fChain->SetBranchAddress("jlmconst", &jlmconst, &b_jlmconst);
   fChain->SetBranchAddress("jlmrho", &jlmrho, &b_jlmrho);
   fChain->SetBranchAddress("jlmsig", &jlmsig, &b_jlmsig);
   fChain->SetBranchAddress("jlmarea", &jlmarea, &b_jlmarea);
   fChain->SetBranchAddress("jsconst", &jsconst, &b_jsconst);
   fChain->SetBranchAddress("jsrho", &jsrho, &b_jsrho);
   fChain->SetBranchAddress("jssig", &jssig, &b_jssig);
   fChain->SetBranchAddress("jsarea", &jsarea, &b_jsarea);
   fChain->SetBranchAddress("jsmconst", &jsmconst, &b_jsmconst);
   fChain->SetBranchAddress("jsmrho", &jsmrho, &b_jsmrho);
   fChain->SetBranchAddress("jsmsig", &jsmsig, &b_jsmsig);
   fChain->SetBranchAddress("jsmarea", &jsmarea, &b_jsmarea);
   fChain->SetBranchAddress("embed_eventid", &embed_eventid, &b_embed_eventid);
   fChain->SetBranchAddress("embed_runid", &embed_runid, &b_embed_runid);
   fChain->SetBranchAddress("embed_refmult", &embed_refmult, &b_embed_refmult);
   fChain->SetBranchAddress("embed_grefmult", &embed_grefmult, &b_embed_grefmult);
   fChain->SetBranchAddress("embed_refmultcorr", &embed_refmultcorr, &b_embed_refmultcorr);
   fChain->SetBranchAddress("embed_grefmultcorr", &embed_grefmultcorr, &b_embed_grefmultcorr);
   fChain->SetBranchAddress("embed_cent", &embed_cent, &b_embed_cent);
   fChain->SetBranchAddress("embed_npart", &embed_npart, &b_embed_npart);
   fChain->SetBranchAddress("embed_rp", &embed_rp, &b_embed_rp);
   fChain->SetBranchAddress("embed_zdcrate", &embed_zdcrate, &b_embed_zdcrate);
   fChain->SetBranchAddress("embed_vz", &embed_vz, &b_embed_vz);
   fChain->SetBranchAddress("foundpp", &foundpp, &b_foundpp);
   fChain->SetBranchAddress("ppjl", &ppjl, &b_ppjl);
   fChain->SetBranchAddress("ppjs", &ppjs, &b_ppjs);
   fChain->SetBranchAddress("ppjlm", &ppjlm, &b_ppjlm);
   fChain->SetBranchAddress("ppjsm", &ppjsm, &b_ppjsm);
   Notify();
}

Bool_t PPReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PPReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PPReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

