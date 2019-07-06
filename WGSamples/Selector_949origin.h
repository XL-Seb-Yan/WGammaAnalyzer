//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  2 11:07:06 2019 by ROOT version 6.10/09
// from TTree tree/tree
// found on file: flatTuple_MC1600_949.root
//////////////////////////////////////////////////////////

#ifndef Selector_h
#define Selector_h

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



class Selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> genParticle_N = {fReader, "genParticle_N"};
   TTreeReaderArray<float> genParticle_pt = {fReader, "genParticle_pt"};
   TTreeReaderArray<float> genParticle_eta = {fReader, "genParticle_eta"};
   TTreeReaderArray<float> genParticle_phi = {fReader, "genParticle_phi"};
   TTreeReaderArray<float> genParticle_mass = {fReader, "genParticle_mass"};
   TTreeReaderArray<int> genParticle_pdgId = {fReader, "genParticle_pdgId"};
   TTreeReaderArray<int> genParticle_status = {fReader, "genParticle_status"};
   TTreeReaderArray<int> genParticle_isPrompt = {fReader, "genParticle_isPrompt"};
   TTreeReaderArray<int> genParticle_isDirectPromptTauDecayProduct = {fReader, "genParticle_isDirectPromptTauDecayProduct"};
   TTreeReaderArray<int> genParticle_isDirectHardProcessTauDecayProductFinalState = {fReader, "genParticle_isDirectHardProcessTauDecayProductFinalState"};
   TTreeReaderArray<int> genParticle_fromHardProcessFinalState = {fReader, "genParticle_fromHardProcessFinalState"};
   TTreeReaderArray<vector<int>> genParticle_mother = {fReader, "genParticle_mother"};
   TTreeReaderArray<int> genParticle_nMoth = {fReader, "genParticle_nMoth"};
   TTreeReaderArray<int> genParticle_nDau = {fReader, "genParticle_nDau"};
   TTreeReaderArray<vector<int>> genParticle_dau = {fReader, "genParticle_dau"};
   TTreeReaderArray<float> genParticle_tauvispt = {fReader, "genParticle_tauvispt"};
   TTreeReaderArray<float> genParticle_tauviseta = {fReader, "genParticle_tauviseta"};
   TTreeReaderArray<float> genParticle_tauvisphi = {fReader, "genParticle_tauvisphi"};
   TTreeReaderArray<float> genParticle_tauvismass = {fReader, "genParticle_tauvismass"};
   TTreeReaderArray<int> genParticle_taudecay = {fReader, "genParticle_taudecay"};
   TTreeReaderValue<Float_t> lheV_pt = {fReader, "lheV_pt"};
   TTreeReaderValue<Float_t> lheHT = {fReader, "lheHT"};
   TTreeReaderValue<Int_t> lheNj = {fReader, "lheNj"};
   TTreeReaderValue<Int_t> lheNb = {fReader, "lheNb"};
   TTreeReaderValue<Int_t> lheNl = {fReader, "lheNl"};
   TTreeReaderValue<Float_t> lheV_mass = {fReader, "lheV_mass"};
   TTreeReaderValue<Float_t> genWeight = {fReader, "genWeight"};
   TTreeReaderValue<Float_t> genFacWeightUp = {fReader, "genFacWeightUp"};
   TTreeReaderValue<Float_t> genFacWeightDown = {fReader, "genFacWeightDown"};
   TTreeReaderValue<Float_t> genRenWeightUp = {fReader, "genRenWeightUp"};
   TTreeReaderValue<Float_t> genRenWeightDown = {fReader, "genRenWeightDown"};
   TTreeReaderValue<Float_t> genFacRenWeightUp = {fReader, "genFacRenWeightUp"};
   TTreeReaderValue<Float_t> genFacRenWeightDown = {fReader, "genFacRenWeightDown"};
   TTreeReaderValue<Float_t> qScale = {fReader, "qScale"};
   TTreeReaderValue<Float_t> PDF_rms = {fReader, "PDF_rms"};
   TTreeReaderArray<float> PDF_x = {fReader, "PDF_x"};
   TTreeReaderArray<float> PDF_xPDF = {fReader, "PDF_xPDF"};
   TTreeReaderArray<int> PDF_id = {fReader, "PDF_id"};
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
   TTreeReaderValue<vector<bool>> jetAK4_IDTightWithoutLepVeto = {fReader, "jetAK4_IDTightWithoutLepVeto"};
   TTreeReaderArray<float> jetAK4_PUIDdiscriminat = {fReader, "jetAK4_PUIDdiscriminat"};
   TTreeReaderValue<vector<bool>> jetAK4_PUIDloose = {fReader, "jetAK4_PUIDloose"};
   TTreeReaderValue<vector<bool>> jetAK4_PUIDmedium = {fReader, "jetAK4_PUIDmedium"};
   TTreeReaderValue<vector<bool>> jetAK4_PUIDtight = {fReader, "jetAK4_PUIDtight"};
   TTreeReaderArray<float> jetAK4_muf = {fReader, "jetAK4_muf"};
   TTreeReaderArray<float> jetAK4_phf = {fReader, "jetAK4_phf"};
   TTreeReaderArray<float> jetAK4_emf = {fReader, "jetAK4_emf"};
   TTreeReaderArray<float> jetAK4_nhf = {fReader, "jetAK4_nhf"};
   TTreeReaderArray<float> jetAK4_chf = {fReader, "jetAK4_chf"};
   TTreeReaderArray<float> jetAK4_area = {fReader, "jetAK4_area"};
   TTreeReaderArray<int> jetAK4_cm = {fReader, "jetAK4_cm"};
   TTreeReaderArray<int> jetAK4_nm = {fReader, "jetAK4_nm"};
   TTreeReaderArray<float> jetAK4_che = {fReader, "jetAK4_che"};
   TTreeReaderArray<float> jetAK4_ne = {fReader, "jetAK4_ne"};
   TTreeReaderArray<float> jetAK4_hf_hf = {fReader, "jetAK4_hf_hf"};
   TTreeReaderArray<float> jetAK4_hf_emf = {fReader, "jetAK4_hf_emf"};
   TTreeReaderArray<float> jetAK4_hof = {fReader, "jetAK4_hof"};
   TTreeReaderArray<int> jetAK4_chm = {fReader, "jetAK4_chm"};
   TTreeReaderArray<int> jetAK4_neHadMult = {fReader, "jetAK4_neHadMult"};
   TTreeReaderArray<int> jetAK4_phoMult = {fReader, "jetAK4_phoMult"};
   TTreeReaderArray<float> jetAK4_nemf = {fReader, "jetAK4_nemf"};
   TTreeReaderArray<float> jetAK4_cemf = {fReader, "jetAK4_cemf"};
   TTreeReaderArray<int> jetAK4_charge = {fReader, "jetAK4_charge"};
   TTreeReaderArray<float> jetAK4_csv = {fReader, "jetAK4_csv"};
   TTreeReaderArray<float> jetAK4_deep_csv_b = {fReader, "jetAK4_deep_csv_b"};
   TTreeReaderArray<float> jetAK4_deep_csv_bb = {fReader, "jetAK4_deep_csv_bb"};
   TTreeReaderArray<float> jetAK4_vtxMass = {fReader, "jetAK4_vtxMass"};
   TTreeReaderArray<float> jetAK4_vtxNtracks = {fReader, "jetAK4_vtxNtracks"};
   TTreeReaderArray<float> jetAK4_vtx3DVal = {fReader, "jetAK4_vtx3DVal"};
   TTreeReaderArray<float> jetAK4_vtx3DSig = {fReader, "jetAK4_vtx3DSig"};
   TTreeReaderArray<float> jetAK4_etaAxis = {fReader, "jetAK4_etaAxis"};
   TTreeReaderArray<float> jetAK4_phiAxis = {fReader, "jetAK4_phiAxis"};
   TTreeReaderArray<float> jetAK4_phiT = {fReader, "jetAK4_phiT"};
   TTreeReaderArray<float> jetAK4_qg_axis1 = {fReader, "jetAK4_qg_axis1"};
   TTreeReaderArray<float> jetAK4_qg_axis2 = {fReader, "jetAK4_qg_axis2"};
   TTreeReaderArray<int> jetAK4_qg_charged = {fReader, "jetAK4_qg_charged"};
   TTreeReaderArray<float> jetAK4_qg_ptD = {fReader, "jetAK4_qg_ptD"};
   TTreeReaderArray<float> jetAK4_qg_pt_dr = {fReader, "jetAK4_qg_pt_dr"};
   TTreeReaderArray<int> jetAK4_partonFlavour = {fReader, "jetAK4_partonFlavour"};
   TTreeReaderArray<int> jetAK4_hadronFlavour = {fReader, "jetAK4_hadronFlavour"};
   TTreeReaderArray<int> jetAK4_genParton_pdgID = {fReader, "jetAK4_genParton_pdgID"};
   TTreeReaderArray<int> jetAK4_nbHadrons = {fReader, "jetAK4_nbHadrons"};
   TTreeReaderArray<int> jetAK4_ncHadrons = {fReader, "jetAK4_ncHadrons"};
   TTreeReaderArray<float> jetAK4_jer_sf = {fReader, "jetAK4_jer_sf"};
   TTreeReaderArray<float> jetAK4_jer_sf_up = {fReader, "jetAK4_jer_sf_up"};
   TTreeReaderArray<float> jetAK4_jer_sf_down = {fReader, "jetAK4_jer_sf_down"};
   TTreeReaderArray<float> jetAK4_jer_sigma_pt = {fReader, "jetAK4_jer_sigma_pt"};
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
   TTreeReaderArray<float> jetAK8_muf = {fReader, "jetAK8_muf"};
   TTreeReaderArray<float> jetAK8_phf = {fReader, "jetAK8_phf"};
   TTreeReaderArray<float> jetAK8_emf = {fReader, "jetAK8_emf"};
   TTreeReaderArray<float> jetAK8_nhf = {fReader, "jetAK8_nhf"};
   TTreeReaderArray<float> jetAK8_chf = {fReader, "jetAK8_chf"};
   TTreeReaderArray<float> jetAK8_area = {fReader, "jetAK8_area"};
   TTreeReaderArray<int> jetAK8_cm = {fReader, "jetAK8_cm"};
   TTreeReaderArray<int> jetAK8_nm = {fReader, "jetAK8_nm"};
   TTreeReaderArray<float> jetAK8_che = {fReader, "jetAK8_che"};
   TTreeReaderArray<float> jetAK8_ne = {fReader, "jetAK8_ne"};
   TTreeReaderArray<float> jetAK8_hf_hf = {fReader, "jetAK8_hf_hf"};
   TTreeReaderArray<float> jetAK8_hf_emf = {fReader, "jetAK8_hf_emf"};
   TTreeReaderArray<float> jetAK8_hof = {fReader, "jetAK8_hof"};
   TTreeReaderArray<int> jetAK8_chm = {fReader, "jetAK8_chm"};
   TTreeReaderArray<int> jetAK8_neHadMult = {fReader, "jetAK8_neHadMult"};
   TTreeReaderArray<int> jetAK8_phoMult = {fReader, "jetAK8_phoMult"};
   TTreeReaderArray<float> jetAK8_nemf = {fReader, "jetAK8_nemf"};
   TTreeReaderArray<float> jetAK8_cemf = {fReader, "jetAK8_cemf"};
   TTreeReaderArray<int> jetAK8_charge = {fReader, "jetAK8_charge"};
   TTreeReaderArray<int> jetAK8_partonFlavour = {fReader, "jetAK8_partonFlavour"};
   TTreeReaderArray<int> jetAK8_hadronFlavour = {fReader, "jetAK8_hadronFlavour"};
   TTreeReaderArray<int> jetAK8_genParton_pdgID = {fReader, "jetAK8_genParton_pdgID"};
   TTreeReaderArray<int> jetAK8_nbHadrons = {fReader, "jetAK8_nbHadrons"};
   TTreeReaderArray<int> jetAK8_ncHadrons = {fReader, "jetAK8_ncHadrons"};
   TTreeReaderArray<float> jetAK8_jer_sf = {fReader, "jetAK8_jer_sf"};
   TTreeReaderArray<float> jetAK8_jer_sf_up = {fReader, "jetAK8_jer_sf_up"};
   TTreeReaderArray<float> jetAK8_jer_sf_down = {fReader, "jetAK8_jer_sf_down"};
   TTreeReaderArray<float> jetAK8_jer_sigma_pt = {fReader, "jetAK8_jer_sigma_pt"};
   TTreeReaderArray<float> jetAK8_csv = {fReader, "jetAK8_csv"};
   TTreeReaderArray<float> jetAK8_deep_csv_b = {fReader, "jetAK8_deep_csv_b"};
   TTreeReaderArray<float> jetAK8_deep_csv_bb = {fReader, "jetAK8_deep_csv_bb"};
   TTreeReaderArray<float> jetAK8_tau1 = {fReader, "jetAK8_tau1"};
   TTreeReaderArray<float> jetAK8_tau2 = {fReader, "jetAK8_tau2"};
   TTreeReaderArray<float> jetAK8_tau3 = {fReader, "jetAK8_tau3"};
   TTreeReaderArray<float> jetAK8_pull1 = {fReader, "jetAK8_pull1"};
   TTreeReaderArray<float> jetAK8_pull2 = {fReader, "jetAK8_pull2"};
   TTreeReaderArray<float> jetAK8_chs_pruned_mass = {fReader, "jetAK8_chs_pruned_mass"};
   TTreeReaderArray<float> jetAK8_chs_softdrop_mass = {fReader, "jetAK8_chs_softdrop_mass"};
   TTreeReaderArray<float> jetAK8_chs_pruned_massCorr = {fReader, "jetAK8_chs_pruned_massCorr"};
   TTreeReaderArray<float> jetAK8_chs_pruned_jec = {fReader, "jetAK8_chs_pruned_jec"};
   TTreeReaderArray<float> jetAK8_chs_pruned_jecUp = {fReader, "jetAK8_chs_pruned_jecUp"};
   TTreeReaderArray<float> jetAK8_chs_pruned_jecDown = {fReader, "jetAK8_chs_pruned_jecDown"};
   TTreeReaderArray<float> jetAK8_chs_tau1 = {fReader, "jetAK8_chs_tau1"};
   TTreeReaderArray<float> jetAK8_chs_tau2 = {fReader, "jetAK8_chs_tau2"};
   TTreeReaderArray<float> jetAK8_chs_tau3 = {fReader, "jetAK8_chs_tau3"};
   TTreeReaderArray<float> jetAK8_softdrop_mass = {fReader, "jetAK8_softdrop_mass"};
   TTreeReaderArray<float> jetAK8_softdrop_massCorr = {fReader, "jetAK8_softdrop_massCorr"};
   TTreeReaderArray<float> jetAK8_softdrop_jec = {fReader, "jetAK8_softdrop_jec"};
   TTreeReaderArray<float> jetAK8_softdrop_jecUp = {fReader, "jetAK8_softdrop_jecUp"};
   TTreeReaderArray<float> jetAK8_softdrop_jecDown = {fReader, "jetAK8_softdrop_jecDown"};
   TTreeReaderArray<int> jetAK8_subjet_puppi_softdrop_N = {fReader, "jetAK8_subjet_puppi_softdrop_N"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_pt = {fReader, "jetAK8_subjet_puppi_softdrop_pt"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_eta = {fReader, "jetAK8_subjet_puppi_softdrop_eta"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_mass = {fReader, "jetAK8_subjet_puppi_softdrop_mass"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_phi = {fReader, "jetAK8_subjet_puppi_softdrop_phi"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_e = {fReader, "jetAK8_subjet_puppi_softdrop_e"};
   TTreeReaderArray<vector<int>> jetAK8_subjet_puppi_softdrop_charge = {fReader, "jetAK8_subjet_puppi_softdrop_charge"};
   TTreeReaderArray<vector<int>> jetAK8_subjet_puppi_softdrop_genParton_pdgID = {fReader, "jetAK8_subjet_puppi_softdrop_genParton_pdgID"};
   TTreeReaderArray<vector<int>> jetAK8_subjet_puppi_softdrop_nbHadrons = {fReader, "jetAK8_subjet_puppi_softdrop_nbHadrons"};
   TTreeReaderArray<vector<int>> jetAK8_subjet_puppi_softdrop_ncHadrons = {fReader, "jetAK8_subjet_puppi_softdrop_ncHadrons"};
   TTreeReaderArray<vector<int>> jetAK8_subjet_puppi_softdrop_partonFlavour = {fReader, "jetAK8_subjet_puppi_softdrop_partonFlavour"};
   TTreeReaderArray<vector<int>> jetAK8_subjet_puppi_softdrop_hadronFlavour = {fReader, "jetAK8_subjet_puppi_softdrop_hadronFlavour"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_csv = {fReader, "jetAK8_subjet_puppi_softdrop_csv"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_deep_csv_b = {fReader, "jetAK8_subjet_puppi_softdrop_deep_csv_b"};
   TTreeReaderArray<vector<float>> jetAK8_subjet_puppi_softdrop_deep_csv_bb = {fReader, "jetAK8_subjet_puppi_softdrop_deep_csv_bb"};


   Selector(TTree * /*tree*/ =0) { }
   virtual ~Selector() { }
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

   ClassDef(Selector,0);

};

#endif

#ifdef Selector_cxx
void Selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Selector_cxx
