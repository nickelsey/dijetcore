//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 29 00:10:07 2019 by ROOT version 6.17/01
// from TTree LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_1.5_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_pt_1.5_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2/LEAD_INIT_R_0.4_alg_2_pt_20_const_eta_1_const_pt_1.5_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2_SUB_INIT_R_0.4_alg_2_pt_10_const_eta_1_const_pt_1.5_MATCH_R_0.4_alg_2_pt_0_const_eta_1_const_pt_0.2
// found on file: trees/y7/y7_grid_20_10.root
//////////////////////////////////////////////////////////

#ifndef DIJETCORE_UTIL_DIJET_IMBALANCE_AUAUREADER_H
#define DIJETCORE_UTIL_DIJET_IMBALANCE_AUAUREADER_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"

namespace dijetcore {

class AuAuReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runid;
   Int_t           eventid;
   Double_t        vz;
   Int_t           refmult;
   Int_t           grefmult;
   Double_t        refmultcorr;
   Double_t        grefmultcorr;
   Int_t           cent;
   Double_t        zdcrate;
   Double_t        rp;
   Int_t           nglobal;
   Int_t           npart;
   TLorentzVector  *jl;
   TLorentzVector  *js;
   TLorentzVector  *jlm;
   TLorentzVector  *jsm;
   TLorentzVector  *jloa;
   TLorentzVector  *jsoa;
   Int_t           jlconst;
   Double_t        jlrho;
   Double_t        jlsig;
   Double_t        jlarea;
   Int_t           jlmconst;
   Double_t        jlmrho;
   Double_t        jlmsig;
   Double_t        jlmarea;
   Int_t           jsconst;
   Double_t        jsrho;
   Double_t        jssig;
   Double_t        jsarea;
   Int_t           jsmconst;
   Double_t        jsmrho;
   Double_t        jsmsig;
   Double_t        jsmarea;
   Int_t           oacent;
   Int_t           jloaconst;
   Double_t        jloarho;
   Double_t        jloasig;
   Int_t           jsoaconst;
   Double_t        jsoarho;
   Double_t        jsoasig;

   // List of branches
   TBranch        *b_runid;   //!
   TBranch        *b_eventid;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_refmult;   //!
   TBranch        *b_grefmult;   //!
   TBranch        *b_refmultcorr;   //!
   TBranch        *b_grefmultcorr;   //!
   TBranch        *b_cent;   //!
   TBranch        *b_zdcrate;   //!
   TBranch        *b_rp;   //!
   TBranch        *b_nglobal;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_jl;   //!
   TBranch        *b_js;   //!
   TBranch        *b_jlm;   //!
   TBranch        *b_jsm;   //!
   TBranch        *b_jloa;   //!
   TBranch        *b_jsoa;   //!
   TBranch        *b_jlconst;   //!
   TBranch        *b_jlrho;   //!
   TBranch        *b_jlsig;   //!
   TBranch        *b_jlarea;   //!
   TBranch        *b_jlmconst;   //!
   TBranch        *b_jlmrho;   //!
   TBranch        *b_jlmsig;   //!
   TBranch        *b_jlmarea;   //!
   TBranch        *b_jsconst;   //!
   TBranch        *b_jsrho;   //!
   TBranch        *b_jssig;   //!
   TBranch        *b_jsarea;   //!
   TBranch        *b_jsmconst;   //!
   TBranch        *b_jsmrho;   //!
   TBranch        *b_jsmsig;   //!
   TBranch        *b_jsmarea;   //!
   TBranch        *b_oacent;   //!
   TBranch        *b_jloaconst;   //!
   TBranch        *b_jloarho;   //!
   TBranch        *b_jloasig;   //!
   TBranch        *b_jsoaconst;   //!
   TBranch        *b_jsoarho;   //!
   TBranch        *b_jsoasig;   //!

   AuAuReader(TTree *tree=0);
   virtual ~AuAuReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

} // namespace dijetcore 

#endif // DIJETCORE_UTIL_DIJET_IMBALANCE_AUAUREADER_H