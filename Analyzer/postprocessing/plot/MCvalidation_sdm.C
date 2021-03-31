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
//#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

void MCvalidation_sdm()
{
  gROOT->SetBatch(1);
  lumi_13TeV = "137.19 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "";
  lumiTextSize = 0.45;
  cmsTextSize = 0.6;
  int iPeriod = 12;
  int iPos = 0;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.025);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetEndErrorSize(10);
  //gStyle->SetErrorX(0.0000001);

  TH1 *hist16 = new TH1F("16","masssoftdrop_{j}",26,68,94);
  TH1 *hist26 = new TH1F("26","masssoftdrop_{j}",26,68,94);
  TH1 *hist36 = new TH1F("36","masssoftdrop_{j}",26,68,94);
  TH1 *hist46 = new TH1F("46","masssoftdrop_{j}",26,68,94);
  TH1 *hist56 = new TH1F("56","masssoftdrop_{j}",26,68,94);
  TH1 *hist66 = new TH1F("66","masssoftdrop_{j}",26,68,94);
  TH1 *hist76 = new TH1F("76","masssoftdrop_{j}",26,68,94);
  TH1 *hist86 = new TH1F("86","masssoftdrop_{j}",26,68,94);
  TH1 *hist96 = new TH1F("96","masssoftdrop_{j}",26,68,94);
  
  // Local variables to store to outfile
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21;
  // System
  float sys_costhetastar, sys_ptoverm, m;
  
  // Open input file
  Float_t p_pt, p_eta, p_phi, p_e, p_mva, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight, x_kfactor, x_puweight, x_sf;
  
  TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/GJets_postproc_WGamma17_full_full_jmcorr_May22.root");
  TTree* theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_kfactor", &x_kfactor);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
	if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if((j_mass > 68 && j_mass < 94)){
    //f((j_mass > 68 && j_mass < 94)){
    //if((j_mass > 68 && j_mass < 94) || (j_mass > 68 && j_mass < 94)){
	  x_kfactor = 1.4;
    hist16->Fill(j_mass, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
	}
  }

  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/QCD_postproc_WGamma17_full_full_jmcorr_May22.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_kfactor", &x_kfactor);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
	if(p_pt < 225) continue;
	if(j_pt < 225) continue;
	x_kfactor = 1.2;
    if((j_mass > 68 && j_mass < 94)){
    //if((j_mass > 68 && j_mass < 94)){
    //if((j_mass > 68 && j_mass < 94) || (j_mass > 68 && j_mass < 94)){
    hist26->Fill(j_mass, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);  
	}
  }

  int SB = 0;
  int WB = 0;
  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/Run2Data_postproc_WGammaRun2_full_full_jmcorr_May22.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  //theTree->SetBranchAddress("xsec_weight", &x_weight);
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
	if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if(j_mass > 68 && j_mass < 94) {
		SB++;
		
    
    hist36->Fill(j_mass, 1);

    }
  }
  
  int count11 = 0;
  int count12 = 0;
  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/spin-0_scalar_21Mar23/presel/1000_N_postproc_WGamma17_full_full_jmcorr_Mar23.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  //theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if(s_mass < 600) continue;
    if(j_mass > 40 && j_mass < 65)
      count11++;
    if(j_mass > 68 && j_mass < 94)
      count12++;
    
    hist46->Fill(j_mass, 400*0.0343 * x_sf * x_puweight);
  }
  cout<<count11<<" "<<count12<<endl;

  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/spin-0_scalar_21Mar23/presel/1000_W_postproc_WGamma17_full_full_jmcorr_Mar23.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  //theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if(s_mass < 600) continue;
    
    hist56->Fill(j_mass, 400*0.0343 * x_sf * x_puweight);  
  }

  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/SignalMC3500N_postproc_WGamma17_full_full_jmcorr_May22.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  //theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if(s_mass < 600) continue;
    
    hist66->Fill(j_mass, 400*0.0343 * x_sf * x_puweight);
  }
  
    input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/SignalMC3500W_postproc_WGamma17_full_full_jmcorr_May22.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  //theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if(s_mass < 600) continue;
    
    hist76->Fill(j_mass, 400*0.0343 * x_sf * x_puweight);
  }
  
  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/SignalMC5000N_postproc_WGamma17_full_full_jmcorr_May22.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  //theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if(s_mass < 600) continue;
    
    hist86->Fill(j_mass, 400*0.0343 * x_sf * x_puweight * 1.6);

  }
  
  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/SignalMC5000W_postproc_WGamma17_full_full_jmcorr_May22.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  //theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if(s_mass < 600) continue;

    hist96->Fill(j_mass, 400*0.0343 * x_sf * x_puweight * 1.0256);
  
  }


  //=================================================================================
  int color1 = 3;
  int color2 = 7;
  int color3 = kOrange;
  
 //bkg
  
  hist16->SetFillColor(color1);
  
  hist26->SetFillColor(color2);

  
  
  hist26->SetFillStyle(3144);

  // Create stack for bkg MC
  THStack *stack5 = new THStack("stack5","Jet mass");

  
  stack5->Add(hist26); stack5->Add(hist16);

  
  // data
  hist36->SetLineColor(kRed);

  
  
  hist36->SetLineWidth(1);

  
  
  hist36->SetMarkerStyle(20);

  
  
  hist36->SetMarkerSize(5);

  
  
  hist36->SetMarkerColor(kRed);

  
  // signal
  
  hist46->SetLineColor(kMagenta+1);
  
  hist56->SetLineColor(kMagenta+1);
  
  hist66->SetLineColor(kOrange+6);
  
  hist76->SetLineColor(kOrange+6);
  
  hist86->SetLineColor(kOrange+10);
  
  hist96->SetLineColor(kOrange+10);


  
  hist46->SetLineStyle(1);
  
  hist56->SetLineStyle(3);
  
  hist66->SetLineStyle(1);
  
  hist76->SetLineStyle(3);
  
  hist86->SetLineStyle(1);
  
  hist96->SetLineStyle(3);
  
  hist46->SetLineWidth(2);
  
  hist56->SetLineWidth(3);
  
  hist66->SetLineWidth(2);
  
  hist76->SetLineWidth(3);
  
  hist86->SetLineWidth(2);
  
  hist96->SetLineWidth(3);

  TLegend *legend = new TLegend(0.19,0.72,0.87,0.87);
  legend->SetBorderSize(0);
  legend->SetNColumns(2);
  legend->SetFillStyle(0);
  TAxis *xaxis1 = NULL;
  TAxis *yaxis1 = NULL;
  // Residual plot
  TH1 *data = NULL;
  TH1 *bkg = NULL;
  TH1 *pull = NULL; 
  TH1 *bkgerr = NULL; 
  TH1 *bkgerrline = NULL;
  TAxis *xaxis2 = NULL;
  TAxis *yaxis2 = NULL;
  TAxis *xaxis3 = NULL;
  TAxis *yaxis3 = NULL;
  
  double lp1 = 0, hp1 = 0.21, lp2 = 0.15, hp2 = 0.34, lp3 = 0.32;
  double main_pad_lower_margin = 0.05;
  double lower_pad1_upper_margin = 0.05;
  double lower_pad1_lower_margin = 0.4;
  double lower_pad2_upper_margin = 0.06;
  double lower_pad2_lower_margin = 0.45;
  
  //===========================================================
  TCanvas *c05 = new TCanvas("c05","",4500,5000);
  c05->cd();
  TPad *p05a = new TPad("p05a","p05a",0.02,0.24,1,1.0);
  TPad *p05b = new TPad("p05b","p05b",0.02,lp1,1,0.25);
  p05a->Draw();
  p05b->Draw();
  p05a->cd();
  p05a->SetBottomMargin(main_pad_lower_margin);
  p05a->SetLeftMargin(0.16);
  p05a->SetRightMargin(0.07);
  xaxis1 = hist36->GetXaxis();
  yaxis1 = hist36->GetYaxis();
  xaxis1->SetTitle("m^{SD}_{J} [GeV]");
  yaxis1->SetTitle("Events / 1 GeV");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.8);
  yaxis1->SetRangeUser(0,12000);
  hist36->Draw("PE1");
  stack5->Draw("SAMEHIST");
  hist46->Draw("SAMEHIST");
  hist56->Draw("SAMEHIST");
  hist66->Draw("SAMEHIST");
  hist76->Draw("SAMEHIST");
  // hist86->Draw("SAMEHIST");
  // hist96->Draw("SAMEHIST");
  hist36->Draw("E1SAME");
  hist36->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist36,"Data SR","e1p");
  legend->AddEntry(hist36,"Data SR","e1p");
  legend->AddEntry(hist16,"#gamma + jet","f");
  legend->AddEntry(hist26,"QCD","f"); 
  legend->AddEntry(hist46,"Sig. 1.0 TeV, #Gamma_{X} / m_{X} = 0.01%");
  legend->AddEntry(hist56,"Sig. 1.0 TeV, #Gamma_{X} / m_{X} = 5%");
  legend->AddEntry(hist66,"Sig. 3.5 TeV, #Gamma_{X} / m_{X} = 0.01%");
  legend->AddEntry(hist76,"Sig. 3.5 TeV, #Gamma_{X} / m_{X} = 5%");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");
  legend->Draw();
  CMS_lumi(p05a,iPeriod,iPos);
  
  p05b->cd();
  p05b->SetTopMargin(lower_pad1_upper_margin);
  p05b->SetBottomMargin(lower_pad1_lower_margin);
  p05b->SetLeftMargin(0.159);
  p05b->SetRightMargin(0.07);
  pull = (TH1*)hist36->Clone();
  bkg = (TH1*)hist16->Clone();
  bkg->Add((TH1*)hist26->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        cout<<bkg->GetBinError(i+1)/bkg->GetBinContent(i+1)<<" "<<bkg->GetBinError(i+1)<<" "<<bkg->GetBinContent(i+1)<<endl;
        bkgerr->SetBinError(i+1,bkg->GetBinError(i+1)/bkg->GetBinContent(i+1) > 0.97 ? 0.97 : bkg->GetBinError(i+1)/bkg->GetBinContent(i+1));
      }
      else{
        bkgerr->SetBinContent(i+1,1);
        bkgerr->SetBinError(i+1,0);
      }
  }
  bkgerr->SetFillStyle(3144);
  bkgerr->SetFillColor(kGray);
  bkgerr->SetMarkerSize(0);
  bkgerr->SetLineWidth(1);
  bkgerrline = (TH1*)bkgerr->Clone();
  bkgerrline->SetFillStyle(0);
  bkgerrline->SetLineColor(kRed);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("m^{SD}_{J} [GeV]");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0.7,1.3);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(0.03);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.13);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.13);
  yaxis2->SetTitleOffset(0.5);
  p05b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  c05->Print("j_m_SR.pdf");
  c05->Print("j_m_SR.svg");
  //==========================================================
}