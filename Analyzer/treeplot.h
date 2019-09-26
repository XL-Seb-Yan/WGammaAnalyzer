//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 24 01:54:13 2019 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: BackgroundMCCombined_Wwindow_presel_weightedTo41p54_fitData.root
//////////////////////////////////////////////////////////

#ifndef treeplot_h
#define treeplot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class treeplot : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Float_t> photon_pt = {fReader, "photon_pt"};
   TTreeReaderValue<Float_t> photon_eta = {fReader, "photon_eta"};
   TTreeReaderValue<Float_t> photon_phi = {fReader, "photon_phi"};
   TTreeReaderValue<Float_t> photon_e = {fReader, "photon_e"};
   TTreeReaderValue<Float_t> ak8puppijet_pt = {fReader, "ak8puppijet_pt"};
   TTreeReaderValue<Float_t> ak8puppijet_eta = {fReader, "ak8puppijet_eta"};
   TTreeReaderValue<Float_t> ak8puppijet_phi = {fReader, "ak8puppijet_phi"};
   TTreeReaderValue<Float_t> ak8puppijet_e = {fReader, "ak8puppijet_e"};
   TTreeReaderValue<Float_t> ak8puppijet_masssoftdropcorr = {fReader, "ak8puppijet_masssoftdropcorr"};
   TTreeReaderValue<Float_t> ak8puppijet_tau21 = {fReader, "ak8puppijet_tau21"};
   TTreeReaderValue<Float_t> sys_costhetastar = {fReader, "sys_costhetastar"};
   TTreeReaderValue<Float_t> sys_ptoverm = {fReader, "sys_ptoverm"};
   TTreeReaderValue<Float_t> sys_invmass = {fReader, "sys_invmass"};
   TTreeReaderValue<Float_t> xsec_weight = {fReader, "xsec_weight"};


   treeplot(TTree * /*tree*/ =0) { }
   virtual ~treeplot() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(treeplot,0);

};

#endif

#ifdef treeplot_cxx
void treeplot::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t treeplot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef treeplot_cxx
