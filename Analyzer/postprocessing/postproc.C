//MC selection
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include "TLorentzVector.h"         // 4-vector class
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TEfficiency.h>
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TRandom.h"
#include <algorithm>
#include <map>
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

void postproc(TString dataset, int runondata)
{

  gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 5;
  int iPos = 0;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.025);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);
  // Plots
   
  // Local variables to store to outfile
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21;
  // System
  float sys_costhetastar, sys_ptoverm, m;
  float xsec_weight, xsec_kfactor, xsec_puweight;
  int sys_pvn;
  
  // Open input file
  Float_t p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight, evt_pdfunc;
  int s_pv; 
  
  TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2018/ntuple_data/"+dataset+"_nominal_pileup_WGamma_full_full_Mar17.root");
  TTree* theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("sys_invmass", &s_mass);
  theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("PV_N", &s_pv);
  
  TH1D* pileup_MC = (TH1D*)input->Get("pileup");
  pileup_MC->Scale(1/(double)pileup_MC->Integral());

  // Create output file
  TFile *outFile = TFile::Open(dataset+"_postproc_WGamma_full_full_Mar17.root", "RECREATE");
  TTree *outTree = new TTree("Events","Events"); 
  outTree->Branch("sys_pvn",       &sys_pvn,      "sys_pvn/I");
  outTree->Branch("photon_pt",       &photon_pt,      "photon_pt/F");
  outTree->Branch("photon_eta",      &photon_eta,      "photon_eta/F");
  outTree->Branch("photon_phi",      &photon_phi,      "photon_phi/F");
  outTree->Branch("photon_e",        &photon_e,      "photon_e/F");
  outTree->Branch("ak8puppijet_pt",       &ak8puppijet_pt,      "ak8puppijet_pt/F");
  outTree->Branch("ak8puppijet_eta",      &ak8puppijet_eta,      "ak8puppijet_eta/F");
  outTree->Branch("ak8puppijet_phi",      &ak8puppijet_phi,      "ak8puppijet_phi/F");
  outTree->Branch("ak8puppijet_e",        &ak8puppijet_e,      "ak8puppijet_e/F");
  outTree->Branch("ak8puppijet_masssoftdropcorr",   &ak8puppijet_masssoftdropcorr,  "ak8puppijet_masssoftdropcorr/F");
  outTree->Branch("ak8puppijet_tau21",              &ak8puppijet_tau21,             "ak8puppijet_tau21/F");
  outTree->Branch("sys_costhetastar",        &sys_costhetastar,      "sys_costhetastar/F");
  outTree->Branch("sys_ptoverm",             &sys_ptoverm,           "sys_ptoverm/F");
  outTree->Branch("m",                       &m,                     "m/F");
  outTree->Branch("xsec_weight",             &xsec_weight,           "xsec_weight/F");
  outTree->Branch("xsec_kfactor",            &xsec_kfactor,           "xsec_kfactor/F");
  outTree->Branch("xsec_puweight",          &xsec_puweight,         "xsec_puweight/F");
  
  TFile* pileup_central = TFile::Open("/afs/cern.ch/user/x/xuyan/WGProj/PROD17/DATA/pileup/Pileup_18.root", "READ");
  TH1F* pileup_Data_central = (TH1F*)pileup_central->Get("pileup");
  pileup_Data_central->Scale(1/(double)pileup_Data_central->Integral());
  pileup_Data_central->Divide(pileup_MC);
 
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    
    photon_pt = p_pt;
    photon_eta = p_eta;
    photon_phi = p_phi;
    photon_e = p_e;
    ak8puppijet_pt = j_pt;
    ak8puppijet_eta = j_eta;
    ak8puppijet_phi = j_phi;
    ak8puppijet_e = j_e;
    ak8puppijet_masssoftdropcorr = j_mass;
    ak8puppijet_tau21 = j_tau21;
    sys_costhetastar = s_cos;
    sys_ptoverm = s_ptm;
    m = s_mass;
    xsec_weight = x_weight;
	sys_pvn = s_pv;
    xsec_kfactor = 1; 
	
	//add pileup-weight
	int binnum = pileup_Data_central->GetXaxis()->FindBin(s_pv);
	double pileupweight_central = pileup_Data_central->GetBinContent(binnum);
	if(runondata)
		xsec_puweight  = 1;
	else
		xsec_puweight = pileupweight_central ;

	outTree->Fill();
  }
  cout<<"Entries: "<<outTree->GetEntries()<<endl;
  outFile->Write();
  outFile->Close();
}