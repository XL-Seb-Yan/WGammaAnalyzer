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

void cutselection(int sigm = 1600)
{

  gROOT->SetBatch(1);
  lumi_13TeV = "59.74 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 6;
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
  TH1F* invmass1 = new TH1F("mass1","mass",85,600,4000);
  TH1F* invmass2 = new TH1F("mass2","mass",85,600,4000);
  TH1F* invmass3 = new TH1F("mass3","mass",85,600,4000);
  TH1F* invmass4 = new TH1F("mass4","mass",85,600,4000);
  
  // Local variables to store to outfile
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21;
  // System
  float sys_costhetastar, sys_ptoverm, m;
  float xsec_weight;

  TString signal = std::to_string(sigm);
  
  // Open input file
  Float_t p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight;
  int PV_N; 
  //TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/AnalysisNtuples_Jan12/SinglePhoton2017_WGamma_full_full_Jan12.root");
  //TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/AnalysisNtuples_Jan12/GJets_WGamma_full_full_Jan12.root");
  //TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/AnalysisNtuples_Jan12/QCD_WGamma_full_full_Jan12.root");
  TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2018/pre_sel_signal/SignalMC"+signal+"W_nominal_pileup_WGamma_full_full_Mar17.root");
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

  //for pile-up study
  theTree->SetBranchAddress("PV_N", &PV_N);
  TH1D* pileup_MC = (TH1D*)input->Get("pileup");
  pileup_MC->Scale(1/(double)pileup_MC->Integral());

  // Create output file
  //TFile *outFile = TFile::Open("SinglePhoton2017_WGamma_full_full_Jan12_weightedforplot.root", "RECREATE");
  //TFile *outFile = TFile::Open("GJets_WGamma_full_full_Jan12_tau21.root", "RECREATE");
  //TFile *outFile = TFile::Open("QCD_WGamma_full_full_Jan12_tau21.root", "RECREATE");
  TFile *outFile = TFile::Open("M"+signal+"N_WGamma_sigrange_Wband_Feb26.root", "RECREATE");
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
  
   //for pile-up study
  TFile* pileup_central = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/pileup/Pileup_18.root", "READ");
  TH1D* pileup_Data_central = (TH1D*)pileup_central->Get("pileup");
  pileup_Data_central->Scale(1/(double)pileup_Data_central->Integral());
  pileup_Data_central->Divide(pileup_MC);
  
  TFile* pileup_up = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/pileup/PileupUp_18.root", "READ");
  TH1D* pileup_Data_up = (TH1D*)pileup_up->Get("pileup");
  pileup_Data_up->Scale(1/(double)pileup_Data_up->Integral());
  pileup_Data_up->Divide(pileup_MC);
  
  TFile* pileup_down = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/pileup/PileupDown_18.root", "READ");
  TH1D* pileup_Data_down = (TH1D*)pileup_down->Get("pileup");
  pileup_Data_down->Scale(1/(double)pileup_Data_down->Integral());
  pileup_Data_down->Divide(pileup_MC);
  
  double SumWight_central = 0;
  double SumWight_up = 0;
  double SumWight_down = 0;
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    
    if(s_mass < 0.75*sigm || s_mass > 1.25*sigm) continue;
    if(j_mass < 65 || j_mass > 105) continue;
    if(abs(p_eta) > 1.44) continue;
    if(abs(j_eta) > 2) continue;
    if(j_tau21 > 0.35) continue;
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
    //xsec_weight = x_weight*1.39;//*1.41; //Additional kfactor on top of xsec weight: QCD:0.50, GJets:1.39
    xsec_weight = x_weight*1;
    outTree->Fill();
	
	int binnum = pileup_Data_central->GetXaxis()->FindBin(PV_N);
	double pileupweight_central = pileup_Data_central->GetBinContent(binnum);
	SumWight_central+=pileupweight_central;
	double pileupweight_up = pileup_Data_up->GetBinContent(binnum);
	SumWight_up+=pileupweight_up;
	double pileupweight_down = pileup_Data_down->GetBinContent(binnum);
	SumWight_down+=pileupweight_down;
	
	invmass1->Fill(m);
	invmass2->Fill(m,pileupweight_central);
	invmass3->Fill(m,pileupweight_up);
	invmass4->Fill(m,pileupweight_down);
	
  }
  cout<<"Entries: "<<outTree->GetEntries()<<endl;
  cout<<"SumWeight_central: "<<SumWight_central<<" SumWeight_up: "<<SumWight_up<<" SumWeight_down: "<<SumWight_down<<endl;
  outFile->Write();
  outFile->Close();
  
  TCanvas *c1 = new TCanvas("c1","",2400,1400);
  c1->cd();
  c1->SetLogy();
  c1->SetBottomMargin(0.12);
  invmass2->GetYaxis()->SetRangeUser(0.1,10000);
  invmass2->GetYaxis()->SetTitle("Events / 40 GeV");
  invmass2->GetYaxis()->SetTitleOffset(1);
  invmass2->GetXaxis()->SetTitle("M_{j#gamma}");
  invmass2->GetXaxis()->SetLimits(sigm*0.7,sigm*1.3);
  invmass3->GetXaxis()->SetLimits(sigm*0.7,sigm*1.3);
  invmass4->GetXaxis()->SetLimits(sigm*0.7,sigm*1.3);
  invmass2->SetLineColor(4);
  invmass3->SetLineColor(2);
  invmass4->SetLineColor(8);
  invmass2->Draw();
  invmass3->Draw("SAME");
  invmass4->Draw("SAME");
  TLegend *l = new TLegend(0.60,0.8,0.9,0.9);
  l->AddEntry(invmass2, "Pileup reweighting - Central");
  l->AddEntry(invmass3, "Pileup reweighting - Up");
  l->AddEntry(invmass4, "Pileup reweighting - Down");
  l->Draw();
  CMS_lumi(c1,iPeriod,iPos);
  c1->Print("m"+signal+".png");
  c1->Print("m"+signal+".pdf");
  c1->Print("m"+signal+".svg");
}
