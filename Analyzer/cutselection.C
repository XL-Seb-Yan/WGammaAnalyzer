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
#include "../Utils/interface/ConfParse.hh"             // input conf file parser
#include "../Utils/interface/CSample.hh"      // helper class to handle samples
#include <algorithm>
#include <map>
#endif

void cutselection(int sigm = 1600)
{

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  // Plots
  
  // Local variables to store to outfile
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21;
  // System
  float sys_costhetastar, sys_ptoverm, m;
  float xsec_weight;

  // Create output file
  TString signal = std::to_string(sigm);
  //TFile *outFile = TFile::Open("SinglePhoton2017_WGamma_Wsidebandband_full_finalcut.root", "RECREATE");
  // TFile *outFile = TFile::Open("BackgroundCombinedMC_WGamma_Wband_full_finalcut.root", "RECREATE");
  TFile *outFile = TFile::Open("Signal"+signal+"N_Wwindow_sigrange_finalcut.root", "RECREATE");
  TTree *outTree = new TTree("Events","Events"); 
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

  // Open input file
  Float_t p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight; 
  //TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017_WGamma_full_full.root");
  //TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/BackgroundCombinedMC_WGamma_full_full_weightedTo41p54_fitData.root");
  TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC"+signal+"N_WGamma_full_full.root");
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
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);

    //narrow signals
    //if(s_mass < 0.75*sigm || s_mass > 1.25*sigm) continue;
    
    if(j_mass < 65 || j_mass > 105) continue;
    if(abs(p_eta) > 1.44) continue;
    if(abs(j_eta) > 2) continue;
    if(j_tau21 > 0.3) continue;
    if(s_ptm < 0.37) continue;
    if(s_cos > 0.6) continue;
    
    
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
    //xsec_weight = x_weight;//*1.41; //Additional kfactor on top of xsec weight: QCD:0.57, GJets:1.41
    xsec_weight = 1;
    
    outTree->Fill();
  }
  cout<<"Entries: "<<outTree->GetEntries()<<endl;
  outFile->Write();
  outFile->Close();

  /*
  hist1->SetLineColor(2);
  hist1->SetLineWidth(3);
  hist1->Scale(1/double(hist1->GetEntries()));

  TLegend *legend = new TLegend(0.7,0.75,0.85,0.85);
  TCanvas *c01 = new TCanvas("c01","",1200,900);
  TAxis *xaxis = hist1->GetXaxis();
  TAxis *yaxis = hist1->GetYaxis();
  xaxis->SetTitle("BDT");
  yaxis->SetTitle("Entries / 0.1");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.0001,1);
  c01->SetLogy();
  c01->cd();
  hist1->Draw("HIST");
  c01->Print("BDT.png");
  */
  
  /*
  hist1->SetLineColor(2);
  hist2->SetLineColor(4);
  hist1->SetLineWidth(3);
  hist2->SetLineWidth(3);
  hist1->Scale(1/double(hist1->GetEntries()));
  hist2->Scale(1/double(hist2->GetEntries()));

  TLegend *legend = new TLegend(0.7,0.75,0.85,0.85);
  TCanvas *c01 = new TCanvas("c01","",1200,900);
  TAxis *xaxis = hist2->GetXaxis();
  TAxis *yaxis = hist2->GetYaxis();
  xaxis->SetTitle("pt_{#gamma}/M");
  yaxis->SetTitle("Entries / 0.1");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.2);
  //yaxis->SetRangeUser(0.0001,1);
  //c01->SetLogy();
  c01->cd();
  hist2->Draw("HIST");
  hist1->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hist1,"65-70 GeV" ,"f");
  legend->AddEntry(hist2,"70-100 GeV","f");
  legend->Draw();
  c01->Print("BDT.png");
  */
}
