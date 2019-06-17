//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 24 03:10:56 2019 by ROOT version 6.10/09
// from TTree tree/tree
// found on file: flatTuple_muon_johnC.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>

#include "string"

#include <map>



class test : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> ph_N = {fReader, "ph_N"};
   TTreeReaderArray<int> ph_pdgId = {fReader, "ph_pdgId"};
   TTreeReaderArray<float> ph_charge = {fReader, "ph_charge"};
   TTreeReaderArray<float> ph_e = {fReader, "ph_e"};
   TTreeReaderArray<float> ph_eta = {fReader, "ph_eta"};
   TTreeReaderArray<float> ph_phi = {fReader, "ph_phi"};
   TTreeReaderArray<float> ph_mass = {fReader, "ph_mass"};
   TTreeReaderArray<float> ph_pt = {fReader, "ph_pt"};
   TTreeReaderArray<float> ph_et = {fReader, "ph_et"};
   TTreeReaderArray<float> ph_rho = {fReader, "ph_rho"};
   TTreeReaderArray<float> ph_superCluster_eta = {fReader, "ph_superCluster_eta"};
   TTreeReaderArray<float> ph_superCluster_phi = {fReader, "ph_superCluster_phi"};
   TTreeReaderArray<float> ph_sigmaIetaIeta = {fReader, "ph_sigmaIetaIeta"};
   TTreeReaderArray<float> ph_hOverE = {fReader, "ph_hOverE"};
   TTreeReaderArray<float> ph_isoGamma = {fReader, "ph_isoGamma"};
   TTreeReaderArray<float> ph_isoCh = {fReader, "ph_isoCh"};
   TTreeReaderValue<vector<bool>> ph_passEleVeto = {fReader, "ph_passEleVeto"};
   TTreeReaderArray<int> ph_passLooseId = {fReader, "ph_passLooseId"};
   TTreeReaderArray<int> ph_passMediumId = {fReader, "ph_passMediumId"};
   TTreeReaderArray<int> ph_passTightId = {fReader, "ph_passTightId"};
   TTreeReaderArray<float> ph_mvaVal = {fReader, "ph_mvaVal"};
   TTreeReaderArray<float> ph_mvaCat = {fReader, "ph_mvaCat"};
   TTreeReaderValue<Float_t> rho = {fReader, "rho"};
   TTreeReaderValue<Int_t> jetAK4_N = {fReader, "jetAK4_N"};
   TTreeReaderArray<float> jetAK4_pt = {fReader, "jetAK4_pt"};
   TTreeReaderArray<float> jetAK4_eta = {fReader, "jetAK4_eta"};
   TTreeReaderArray<float> jetAK4_mass = {fReader, "jetAK4_mass"};
   TTreeReaderArray<float> jetAK4_phi = {fReader, "jetAK4_phi"};
   TTreeReaderArray<float> jetAK4_e = {fReader, "jetAK4_e"};
   TTreeReaderArray<float> jetAK4_jec = {fReader, "jetAK4_jec"};
   TTreeReaderArray<float> jetAK4_jecUp = {fReader, "jetAK4_jecUp"};
   TTreeReaderArray<float> jetAK4_jecDown = {fReader, "jetAK4_jecDown"};
   TTreeReaderValue<vector<bool>> jetAK4_IDLoose = {fReader, "jetAK4_IDLoose"};
   TTreeReaderValue<vector<bool>> jetAK4_IDTight = {fReader, "jetAK4_IDTight"};
   TTreeReaderValue<vector<bool>> jetAK4_IDTightLepVeto = {fReader, "jetAK4_IDTightLepVeto"};
   TTreeReaderArray<int> jetAK4_charge = {fReader, "jetAK4_charge"};
   TTreeReaderArray<float> jetAK4_csv = {fReader, "jetAK4_csv"};
   TTreeReaderArray<float> jetAK4_vtxMass = {fReader, "jetAK4_vtxMass"};
   TTreeReaderArray<float> jetAK4_vtxNtracks = {fReader, "jetAK4_vtxNtracks"};
   TTreeReaderArray<float> jetAK4_vtx3DVal = {fReader, "jetAK4_vtx3DVal"};
   TTreeReaderArray<float> jetAK4_vtx3DSig = {fReader, "jetAK4_vtx3DSig"};
   TTreeReaderValue<Int_t> jetAK8_N = {fReader, "jetAK8_N"};
   TTreeReaderArray<float> jetAK8_pt = {fReader, "jetAK8_pt"};
   TTreeReaderArray<float> jetAK8_eta = {fReader, "jetAK8_eta"};
   TTreeReaderArray<float> jetAK8_mass = {fReader, "jetAK8_mass"};
   TTreeReaderArray<float> jetAK8_phi = {fReader, "jetAK8_phi"};
   TTreeReaderArray<float> jetAK8_e = {fReader, "jetAK8_e"};
   TTreeReaderArray<float> jetAK8_jec = {fReader, "jetAK8_jec"};
   TTreeReaderArray<float> jetAK8_jecUp = {fReader, "jetAK8_jecUp"};
   TTreeReaderArray<float> jetAK8_jecDown = {fReader, "jetAK8_jecDown"};
   TTreeReaderValue<vector<bool>> jetAK8_IDLoose = {fReader, "jetAK8_IDLoose"};
   TTreeReaderValue<vector<bool>> jetAK8_IDTight = {fReader, "jetAK8_IDTight"};
   TTreeReaderValue<vector<bool>> jetAK8_IDTightLepVeto = {fReader, "jetAK8_IDTightLepVeto"};
   TTreeReaderArray<int> jetAK8_charge = {fReader, "jetAK8_charge"};
   TTreeReaderArray<float> jetAK8_Hbbtag = {fReader, "jetAK8_Hbbtag"};
   TTreeReaderArray<float> jetAK8_csv = {fReader, "jetAK8_csv"};
   TTreeReaderArray<float> jetAK8_tau1 = {fReader, "jetAK8_tau1"};
   TTreeReaderArray<float> jetAK8_tau2 = {fReader, "jetAK8_tau2"};
   TTreeReaderArray<float> jetAK8_tau3 = {fReader, "jetAK8_tau3"};
   TTreeReaderArray<float> jetAK8_pruned_mass = {fReader, "jetAK8_pruned_mass"};
   TTreeReaderArray<float> jetAK8_pruned_massCorr = {fReader, "jetAK8_pruned_massCorr"};
   TTreeReaderArray<float> jetAK8_pruned_jec = {fReader, "jetAK8_pruned_jec"};
   TTreeReaderArray<float> jetAK8_pruned_jecUp = {fReader, "jetAK8_pruned_jecUp"};
   TTreeReaderArray<float> jetAK8_pruned_jecDown = {fReader, "jetAK8_pruned_jecDown"};
   TTreeReaderArray<float> jetAK8_softdrop_mass = {fReader, "jetAK8_softdrop_mass"};
   TTreeReaderArray<float> jetAK8_softdrop_massCorr = {fReader, "jetAK8_softdrop_massCorr"};
   TTreeReaderArray<float> jetAK8_softdrop_jec = {fReader, "jetAK8_softdrop_jec"};
   TTreeReaderArray<float> jetAK8_softdrop_jecUp = {fReader, "jetAK8_softdrop_jecUp"};
   TTreeReaderArray<float> jetAK8_softdrop_jecDown = {fReader, "jetAK8_softdrop_jecDown"};
   TTreeReaderArray<int> jetAK8_subjet_softdrop_N = {fReader, "jetAK8_subjet_softdrop_N"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_softdrop_pt = {fReader, "jetAK8_subjet_softdrop_pt"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_softdrop_eta = {fReader, "jetAK8_subjet_softdrop_eta"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_softdrop_mass = {fReader, "jetAK8_subjet_softdrop_mass"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_softdrop_phi = {fReader, "jetAK8_subjet_softdrop_phi"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_softdrop_e = {fReader, "jetAK8_subjet_softdrop_e"};
   TTreeReaderArray<vector<int>> jetAK8_subjet_softdrop_charge = {fReader, "jetAK8_subjet_softdrop_charge"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_softdrop_csv = {fReader, "jetAK8_subjet_softdrop_csv"};
   TTreeReaderArray<float> jetAK8_puppi_pt = {fReader, "jetAK8_puppi_pt"};
   TTreeReaderArray<float> jetAK8_puppi_eta = {fReader, "jetAK8_puppi_eta"};
   TTreeReaderArray<float> jetAK8_puppi_mass = {fReader, "jetAK8_puppi_mass"};
   TTreeReaderArray<float> jetAK8_puppi_phi = {fReader, "jetAK8_puppi_phi"};
   TTreeReaderArray<float> jetAK8_puppi_e = {fReader, "jetAK8_puppi_e"};
   TTreeReaderArray<float> jetAK8_puppi_pruned_mass = {fReader, "jetAK8_puppi_pruned_mass"};
   TTreeReaderArray<float> jetAK8_puppi_pruned_massCorr = {fReader, "jetAK8_puppi_pruned_massCorr"};
   TTreeReaderArray<float> jetAK8_puppi_pruned_jec = {fReader, "jetAK8_puppi_pruned_jec"};
   TTreeReaderArray<float> jetAK8_puppi_softdrop_mass = {fReader, "jetAK8_puppi_softdrop_mass"};
   TTreeReaderArray<float> jetAK8_puppi_softdrop_massCorr = {fReader, "jetAK8_puppi_softdrop_massCorr"};
   TTreeReaderArray<float> jetAK8_puppi_softdrop_jec = {fReader, "jetAK8_puppi_softdrop_jec"};
   TTreeReaderArray<float> jetAK8_puppi_tau1 = {fReader, "jetAK8_puppi_tau1"};
   TTreeReaderArray<float> jetAK8_puppi_tau2 = {fReader, "jetAK8_puppi_tau2"};
   TTreeReaderArray<float> jetAK8_puppi_tau3 = {fReader, "jetAK8_puppi_tau3"};
   TTreeReaderArray<int> jetAK8_subjet_puppi_softdrop_N = {fReader, "jetAK8_subjet_puppi_softdrop_N"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_pt = {fReader, "jetAK8_subjet_puppi_softdrop_pt"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_eta = {fReader, "jetAK8_subjet_puppi_softdrop_eta"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_mass = {fReader, "jetAK8_subjet_puppi_softdrop_mass"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_phi = {fReader, "jetAK8_subjet_puppi_softdrop_phi"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_e = {fReader, "jetAK8_subjet_puppi_softdrop_e"};
   TTreeReaderArray<vector<int>> jetAK8_subjet_puppi_softdrop_charge = {fReader, "jetAK8_subjet_puppi_softdrop_charge"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_csv = {fReader, "jetAK8_subjet_puppi_softdrop_csv"};
   TTreeReaderArray<pair<string,bool>> HLT_isFired = {fReader, "HLT_isFired"};
   TTreeReaderValue<Bool_t> passFilter_HBHE_ = {fReader, "passFilter_HBHE"};
   TTreeReaderValue<Bool_t> passFilter_HBHELoose_ = {fReader, "passFilter_HBHELoose"};
   TTreeReaderValue<Bool_t> passFilter_HBHETight_ = {fReader, "passFilter_HBHETight"};
   TTreeReaderValue<Bool_t> passFilter_HBHEIso_ = {fReader, "passFilter_HBHEIso"};
   TTreeReaderValue<Bool_t> passFilter_CSCHalo_ = {fReader, "passFilter_CSCHalo"};
   TTreeReaderValue<Bool_t> passFilter_CSCTightHalo2015_ = {fReader, "passFilter_CSCTightHalo2015"};
   TTreeReaderValue<Bool_t> passFilter_HCALlaser_ = {fReader, "passFilter_HCALlaser"};
   TTreeReaderValue<Bool_t> passFilter_ECALDeadCell_ = {fReader, "passFilter_ECALDeadCell"};
   TTreeReaderValue<Bool_t> passFilter_GoodVtx_ = {fReader, "passFilter_GoodVtx"};
   TTreeReaderValue<Bool_t> passFilter_TrkFailure_ = {fReader, "passFilter_TrkFailure"};
   TTreeReaderValue<Bool_t> passFilter_EEBadSc_ = {fReader, "passFilter_EEBadSc"};
   TTreeReaderValue<Bool_t> passFilter_ECALlaser_ = {fReader, "passFilter_ECALlaser"};
   TTreeReaderValue<Bool_t> passFilter_TrkPOG_ = {fReader, "passFilter_TrkPOG"};
   TTreeReaderValue<Bool_t> passFilter_TrkPOG_manystrip_ = {fReader, "passFilter_TrkPOG_manystrip"};
   TTreeReaderValue<Bool_t> passFilter_TrkPOG_toomanystrip_ = {fReader, "passFilter_TrkPOG_toomanystrip"};
   TTreeReaderValue<Bool_t> passFilter_TrkPOG_logError_ = {fReader, "passFilter_TrkPOG_logError"};
   TTreeReaderValue<Bool_t> passFilter_METFilters_ = {fReader, "passFilter_METFilters"};
   TTreeReaderValue<Bool_t> passFilter_CSCTightHaloTrkMuUnvetoFilter_ = {fReader, "passFilter_CSCTightHaloTrkMuUnvetoFilter"};
   TTreeReaderValue<Bool_t> passFilter_globalTightHalo2016_ = {fReader, "passFilter_globalTightHalo2016"};
   TTreeReaderValue<Bool_t> passFilter_HcalStripHalo_ = {fReader, "passFilter_HcalStripHalo"};
   TTreeReaderValue<Bool_t> passFilter_chargedHadronTrackResolution_ = {fReader, "passFilter_chargedHadronTrackResolution"};
   TTreeReaderValue<Bool_t> passFilter_muonBadTrack_ = {fReader, "passFilter_muonBadTrack"};
   TTreeReaderValue<Int_t> EVENT_event = {fReader, "EVENT_event"};
   TTreeReaderValue<Int_t> EVENT_run = {fReader, "EVENT_run"};
   TTreeReaderValue<Int_t> EVENT_lumiBlock = {fReader, "EVENT_lumiBlock"};
   TTreeReaderValue<Int_t> PV_N = {fReader, "PV_N"};
   TTreeReaderValue<Bool_t> PV_filter = {fReader, "PV_filter"};
   TTreeReaderArray<float> PV_chi2 = {fReader, "PV_chi2"};
   TTreeReaderArray<float> PV_ndof = {fReader, "PV_ndof"};
   TTreeReaderArray<float> PV_rho = {fReader, "PV_rho"};
   TTreeReaderArray<float> PV_z = {fReader, "PV_z"};


   test(TTree * /*tree*/ =0) { }
   virtual ~test() { }
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

   ClassDef(test,0);

};

#endif

#ifdef test_cxx
void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef test_cxx
