#ifndef DIJETCORE_UTIL_DIJET_IMBALANCE_GENERIC_READER_H
#define DIJETCORE_UTIL_DIJET_IMBALANCE_GENERIC_READER_H

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"

namespace dijetcore {

class GenericReader {
public:
  TTree *fChain;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent; //! current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Int_t runid;
  Int_t eventid;
  Double_t vz;
  Int_t refmult;
  Int_t grefmult;
  Double_t refmultcorr;
  Double_t grefmultcorr;
  Int_t cent;
  Double_t zdcrate;
  Double_t rp;
  Int_t nglobal;
  Int_t npart;
  TLorentzVector *jl;
  TLorentzVector *js;
  TLorentzVector *jlm;
  TLorentzVector *jsm;
  Int_t jlconst;
  Double_t jlrho;
  Double_t jlsig;
  Double_t jlarea;
  Int_t jlmconst;
  Double_t jlmrho;
  Double_t jlmsig;
  Double_t jlmarea;
  Int_t jsconst;
  Double_t jsrho;
  Double_t jssig;
  Double_t jsarea;
  Int_t jsmconst;
  Double_t jsmrho;
  Double_t jsmsig;
  Double_t jsmarea;

  // auau only variables
  TLorentzVector *jloa;
  TLorentzVector *jsoa;
  Int_t oacent;
  Int_t jloaconst;
  Double_t jloarho;
  Double_t jloasig;
  Int_t jsoaconst;
  Double_t jsoarho;
  Double_t jsoasig;

  // pp only variables
  Int_t embed_eventid;
  Int_t embed_runid;
  Int_t embed_refmult;
  Int_t embed_grefmult;
  Double_t embed_refmultcorr;
  Double_t embed_grefmultcorr;
  Int_t embed_cent;
  Int_t embed_npart;
  Double_t embed_rp;
  Double_t embed_zdcrate;
  Double_t embed_vz;
  Bool_t foundpp;
  TLorentzVector *ppjl;
  TLorentzVector *ppjs;
  TLorentzVector *ppjlm;
  TLorentzVector *ppjsm;

  // List of branches
  TBranch *b_runid;        //!
  TBranch *b_eventid;      //!
  TBranch *b_vz;           //!
  TBranch *b_refmult;      //!
  TBranch *b_grefmult;     //!
  TBranch *b_refmultcorr;  //!
  TBranch *b_grefmultcorr; //!
  TBranch *b_cent;         //!
  TBranch *b_zdcrate;      //!
  TBranch *b_rp;           //!
  TBranch *b_nglobal;      //!
  TBranch *b_npart;        //!
  TBranch *b_jl;           //!
  TBranch *b_js;           //!
  TBranch *b_jlm;          //!
  TBranch *b_jsm;          //!
  TBranch *b_jlconst;      //!
  TBranch *b_jlrho;        //!
  TBranch *b_jlsig;        //!
  TBranch *b_jlarea;       //!
  TBranch *b_jlmconst;     //!
  TBranch *b_jlmrho;       //!
  TBranch *b_jlmsig;       //!
  TBranch *b_jlmarea;      //!
  TBranch *b_jsconst;      //!
  TBranch *b_jsrho;        //!
  TBranch *b_jssig;        //!
  TBranch *b_jsarea;       //!
  TBranch *b_jsmconst;     //!
  TBranch *b_jsmrho;       //!
  TBranch *b_jsmsig;       //!
  TBranch *b_jsmarea;      //!

  // auau only branches
  TBranch *b_jloa;      //!
  TBranch *b_jsoa;      //!
  TBranch *b_oacent;    //!
  TBranch *b_jloaconst; //!
  TBranch *b_jloarho;   //!
  TBranch *b_jloasig;   //!
  TBranch *b_jsoaconst; //!
  TBranch *b_jsoarho;   //!
  TBranch *b_jsoasig;   //!

  // pp only branches
  TBranch *b_embed_eventid;      //!
  TBranch *b_embed_runid;        //!
  TBranch *b_embed_refmult;      //!
  TBranch *b_embed_grefmult;     //!
  TBranch *b_embed_refmultcorr;  //!
  TBranch *b_embed_grefmultcorr; //!
  TBranch *b_embed_cent;         //!
  TBranch *b_embed_npart;        //!
  TBranch *b_embed_rp;           //!
  TBranch *b_embed_zdcrate;      //!
  TBranch *b_embed_vz;           //!
  TBranch *b_foundpp;            //!
  TBranch *b_ppjl;               //!
  TBranch *b_ppjs;               //!
  TBranch *b_ppjlm;              //!
  TBranch *b_ppjsm;              //!

  GenericReader(TTree *tree, bool off_axis = false, bool pp = false,
                bool embed = false);
  virtual ~GenericReader();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  virtual bool Next();
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree, bool off_axis = false, bool pp = false,
                bool embed = false);
  virtual void Loop();
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
};

} // namespace dijetcore

#endif // DIJETCORE_UTIL_DIJET_IMBALANCE_GENERIC_READER_H