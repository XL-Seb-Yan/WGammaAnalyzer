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

void MCvalidation_S01()
{
  gROOT->SetBatch(1);
  lumi_13TeV = "137 fb^{-1}";
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
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetEndErrorSize(10);
  //gStyle->SetErrorX(0.0000001);

  TH1 *hist11 = new TH1F("11","pt_{#gamma}",70,0,2800);
  TH1 *hist12 = new TH1F("12","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist13 = new TH1F("13","pt_{j}",70,0,2800);
  TH1 *hist14 = new TH1F("14","eta_{j}",50,-2.5,2.5);
  TH1 *hist15 = new TH1F("15","E_{j}",100,0,3000);
  TH1 *hist16 = new TH1F("16","masssoftdrop_{j}",80,30,110);
  TH1 *hist17 = new TH1F("17","tau21_{j}",50,0,1);
  TH1 *hist18 = new TH1F("18","cos(#theta*)_{p}",50,0,1);
  TH1 *hist19 = new TH1F("19","pt/M",50,0,2);
  TH1 *hist110 = new TH1F("110","invariant mass",140,400,6000);
  TH1 *hist111 = new TH1F("111","seperation",50,0,8);
  TH1 *hist21 = new TH1F("21","pt_{#gamma}",70,0,2800);
  TH1 *hist22 = new TH1F("22","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist23 = new TH1F("23","pt_{j}",70,0,2800);
  TH1 *hist24 = new TH1F("24","eta_{j}",50,-2.5,2.5);
  TH1 *hist25 = new TH1F("25","E_{j}",100,0,3000);
  TH1 *hist26 = new TH1F("26","masssoftdrop_{j}",80,30,110);
  TH1 *hist27 = new TH1F("27","tau21_{j}",50,0,1);
  TH1 *hist28 = new TH1F("28","cos(#theta*)_{p}",50,0,1);
  TH1 *hist29 = new TH1F("29","pt/M",50,0,2);
  TH1 *hist210 = new TH1F("210","invariant mass",140,400,6000);
  TH1 *hist211 = new TH1F("211","seperation",50,0,8);
  TH1 *hist31 = new TH1F("31","pt_{#gamma}",70,0,2800);
  TH1 *hist32 = new TH1F("32","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist33 = new TH1F("33","pt_{j}",70,0,2800);
  TH1 *hist34 = new TH1F("34","eta_{j}",50,-2.5,2.5);
  TH1 *hist35 = new TH1F("35","E_{j}",100,0,3000);
  TH1 *hist36 = new TH1F("36","masssoftdrop_{j}",80,30,110);
  TH1 *hist37 = new TH1F("37","tau21_{j}",50,0,1);
  TH1 *hist38 = new TH1F("38","cos(#theta*)_{p}",50,0,1);
  TH1 *hist39 = new TH1F("39","pt/M",50,0,2);
  TH1 *hist310 = new TH1F("310","invariant mass",140,400,6000);
  TH1 *hist311 = new TH1F("311","seperation",50,0,8);
  TH1 *hist3a1 = new TH1F("3a1","pt_{#gamma}",70,0,2800);
  TH1 *hist3a2 = new TH1F("3a2","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist3a3 = new TH1F("3a3","pt_{j}",70,0,2800);
  TH1 *hist3a4 = new TH1F("3a4","eta_{j}",50,-2.5,2.5);
  TH1 *hist3a5 = new TH1F("3a5","E_{j}",100,0,3000);
  TH1 *hist3a6 = new TH1F("3a6","masssoftdrop_{j}",80,30,110);
  TH1 *hist3a7 = new TH1F("3a7","tau21_{j}",50,0,1);
  TH1 *hist3a8 = new TH1F("3a8","cos(#theta*)_{p}",50,0,1);
  TH1 *hist3a9 = new TH1F("3a9","pt/M",50,0,2);
  TH1 *hist3a10 = new TH1F("3a10","invariant mass",140,400,6000);
  TH1 *hist3a11 = new TH1F("3a11","seperation",50,0,8);
  TH1 *hist41 = new TH1F("41","pt_{#gamma}",70,0,2800);
  TH1 *hist42 = new TH1F("42","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist43 = new TH1F("43","pt_{j}",70,0,2800);
  TH1 *hist44 = new TH1F("44","eta_{j}",50,-2.5,2.5);
  TH1 *hist45 = new TH1F("45","E_{j}",100,0,3000);
  TH1 *hist46 = new TH1F("46","masssoftdrop_{j}",80,30,110);
  TH1 *hist47 = new TH1F("47","tau21_{j}",50,0,1);
  TH1 *hist48 = new TH1F("48","cos(#theta*)_{p}",50,0,1);
  TH1 *hist49 = new TH1F("49","pt/M",50,0,2);
  TH1 *hist410 = new TH1F("410","invariant mass",140,400,6000);
  TH1 *hist411 = new TH1F("411","seperation",50,0,8);
  TH1 *hist51 = new TH1F("51","pt_{#gamma}",70,0,2800);
  TH1 *hist52 = new TH1F("52","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist53 = new TH1F("53","pt_{j}",70,0,2800);
  TH1 *hist54 = new TH1F("54","eta_{j}",50,-2.5,2.5);
  TH1 *hist55 = new TH1F("55","E_{j}",100,0,3000);
  TH1 *hist56 = new TH1F("56","masssoftdrop_{j}",80,30,110);
  TH1 *hist57 = new TH1F("57","tau21_{j}",50,0,1);
  TH1 *hist58 = new TH1F("58","cos(#theta*)_{p}",50,0,1);
  TH1 *hist59 = new TH1F("59","pt/M",50,0,2);
  TH1 *hist510 = new TH1F("510","invariant mass",140,400,6000);
  TH1 *hist511 = new TH1F("511","seperation",50,0,8);
  TH1 *hist61 = new TH1F("61","pt_{#gamma}",70,0,2800);
  TH1 *hist62 = new TH1F("62","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist63 = new TH1F("63","pt_{j}",70,0,2800);
  TH1 *hist64 = new TH1F("64","eta_{j}",50,-2.5,2.5);
  TH1 *hist65 = new TH1F("65","E_{j}",100,0,3000);
  TH1 *hist66 = new TH1F("66","masssoftdrop_{j}",80,30,110);
  TH1 *hist67 = new TH1F("67","tau21_{j}",50,0,1);
  TH1 *hist68 = new TH1F("68","cos(#theta*)_{p}",50,0,1);
  TH1 *hist69 = new TH1F("69","pt/M",50,0,2);
  TH1 *hist610 = new TH1F("610","invariant mass",140,400,6000);
  TH1 *hist611 = new TH1F("611","seperation",50,0,8);
  TH1 *hist71 = new TH1F("71","pt_{#gamma}",70,0,2800);
  TH1 *hist72 = new TH1F("72","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist73 = new TH1F("73","pt_{j}",70,0,2800);
  TH1 *hist74 = new TH1F("74","eta_{j}",50,-2.5,2.5);
  TH1 *hist75 = new TH1F("75","E_{j}",100,0,3000);
  TH1 *hist76 = new TH1F("76","masssoftdrop_{j}",80,30,110);
  TH1 *hist77 = new TH1F("77","tau21_{j}",50,0,1);
  TH1 *hist78 = new TH1F("78","cos(#theta*)_{p}",50,0,1);
  TH1 *hist79 = new TH1F("79","pt/M",50,0,2);
  TH1 *hist710 = new TH1F("710","invariant mass",140,400,6000);
  TH1 *hist711 = new TH1F("711","seperation",50,0,8);
  TH1 *hist81 = new TH1F("81","pt_{#gamma}",70,0,2800);
  TH1 *hist82 = new TH1F("82","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist83 = new TH1F("83","pt_{j}",70,0,2800);
  TH1 *hist84 = new TH1F("84","eta_{j}",50,-2.5,2.5);
  TH1 *hist85 = new TH1F("85","E_{j}",100,0,3000);
  TH1 *hist86 = new TH1F("86","masssoftdrop_{j}",80,30,110);
  TH1 *hist87 = new TH1F("87","tau21_{j}",50,0,1);
  TH1 *hist88 = new TH1F("88","cos(#theta*)_{p}",50,0,1);
  TH1 *hist89 = new TH1F("89","pt/M",50,0,2);
  TH1 *hist810 = new TH1F("810","invariant mass",140,400,6000);
  TH1 *hist811 = new TH1F("811","seperation",50,0,8);
  TH1 *hist91 = new TH1F("91","pt_{#gamma}",70,0,2800);
  TH1 *hist92 = new TH1F("92","eta_{#gamma}",58,-1.44,1.44);
  TH1 *hist93 = new TH1F("93","pt_{j}",70,0,2800);
  TH1 *hist94 = new TH1F("94","eta_{j}",50,-2.5,2.5);
  TH1 *hist95 = new TH1F("95","E_{j}",100,0,3000);
  TH1 *hist96 = new TH1F("96","masssoftdrop_{j}",80,30,110);
  TH1 *hist97 = new TH1F("97","tau21_{j}",50,0,1);
  TH1 *hist98 = new TH1F("98","cos(#theta*)_{p}",50,0,1);
  TH1 *hist99 = new TH1F("99","pt/M",50,0,2);
  TH1 *hist910 = new TH1F("910","invariant mass",140,400,6000);
  TH1 *hist911 = new TH1F("911","seperation",50,0,8);
  
  TH1 *hist112 = new TH1F("112","pt_{j#gamma}",100,0,3000);
  TH1 *hist212 = new TH1F("212","pt_{j#gamma}",100,0,3000);
  TH1 *hist312 = new TH1F("312","pt_{j#gamma}",100,0,3000);
  TH1 *hist3a12 = new TH1F("3a12","pt_{j#gamma}",100,0,3000);
  TH1 *hist412 = new TH1F("412","pt_{j#gamma}",100,0,3000);
  TH1 *hist512 = new TH1F("512","pt_{j#gamma}",100,0,3000);
  TH1 *hist612 = new TH1F("612","pt_{j#gamma}",100,0,3000);
  TH1 *hist712 = new TH1F("712","pt_{j#gamma}",100,0,3000);
  TH1 *hist812 = new TH1F("812","pt_{j#gamma}",100,0,3000);
  TH1 *hist912 = new TH1F("912","pt_{j#gamma}",100,0,3000);
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
    //if((j_mass > 68 && j_mass < 94) || (j_mass > 40 && j_mass < 65)){
	  x_kfactor = 1.4;
    hist11->Fill(p_pt, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist12->Fill(p_eta, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist13->Fill(j_pt, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist14->Fill(j_eta, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist15->Fill(j_e, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist16->Fill(j_mass, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist17->Fill(j_tau21, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist18->Fill(s_cos, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist19->Fill(s_ptm, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist110->Fill(s_mass, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist111->Fill(p_mva, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);  
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist112->Fill((p+j).Pt(), x_weight * x_kfactor * x_sf * x_puweight * 3.30339); 
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
    //if((j_mass > 68 && j_mass < 94) || (j_mass > 40 && j_mass < 65)){
    hist21->Fill(p_pt, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist22->Fill(p_eta, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist23->Fill(j_pt, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist24->Fill(j_eta, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist25->Fill(j_e, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist26->Fill(j_mass, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist27->Fill(j_tau21, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist28->Fill(s_cos, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist29->Fill(s_ptm, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist210->Fill(s_mass, x_weight * x_kfactor * x_sf * x_puweight * 3.30339);
    hist211->Fill(p_mva, x_weight * x_kfactor * x_sf * x_puweight * 3.30339); 
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist212->Fill((p+j).Pt(), x_weight * x_kfactor * x_sf * x_puweight * 3.30339);     
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
    if(j_mass > 40 && j_mass < 65) {
		SB++;
		hist31->Fill(p_pt, 0.8597);
		hist32->Fill(p_eta, 0.8597);
		hist33->Fill(j_pt, 0.8597);
		hist34->Fill(j_eta, 0.8597);
		hist35->Fill(j_e, 0.8597);
		hist36->Fill(j_mass, 0.8597);
		hist37->Fill(j_tau21, 0.8597);
		hist38->Fill(s_cos, 0.8597);
		hist39->Fill(s_ptm, 0.8597);
		hist310->Fill(s_mass, 0.8597);
		hist311->Fill(p_mva, 0.8597); 
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist312->Fill((p+j).Pt(), 0.8597); 
		
    // hist31->Fill(p_pt, 1);
    // hist32->Fill(p_eta, 1);
    // hist33->Fill(j_pt, 1);
    // hist34->Fill(j_eta, 1);
    // hist35->Fill(j_e, 1);
    // hist36->Fill(j_mass, 1);
    // hist37->Fill(j_tau21, 1);
    // hist38->Fill(s_cos, 1);
    // hist39->Fill(s_ptm, 1);
    // hist310->Fill(s_mass, 1);
    // hist311->Fill(p_mva, 1);  
    }
	// }
	if(j_mass > 68 && j_mass < 94) {
		WB++;
		hist3a1->Fill(p_pt);
		hist3a2->Fill(p_eta);
		hist3a3->Fill(j_pt);
		hist3a4->Fill(j_eta);
		hist3a5->Fill(j_e);
		hist3a6->Fill(j_mass);
		hist3a7->Fill(j_tau21);
		hist3a8->Fill(s_cos);
		hist3a9->Fill(s_ptm);
		hist3a10->Fill(s_mass);
		hist3a11->Fill(p_mva); 
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist3a12->Fill((p+j).Pt()); 
	}
  }
  

  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/spin-0_scalar_21Mar23/presel/2000_N_postproc_WGamma17_full_full_jmcorr_Mar23.root");
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
    hist41->Fill(p_pt, 0.0343 * x_sf * x_puweight); //* 1.538 due to lack of finished MC samples
    hist42->Fill(p_eta, 0.0343 * x_sf * x_puweight);
    hist43->Fill(j_pt, 0.0343 * x_sf * x_puweight);
    hist44->Fill(j_eta, 0.0343 * x_sf * x_puweight);
    hist45->Fill(j_e, 0.0343 * x_sf * x_puweight);
    hist46->Fill(j_mass, 400*0.0343 * x_sf * x_puweight);
    hist47->Fill(j_tau21, 400*0.0343 * x_sf * x_puweight);
    hist48->Fill(s_cos, 0.0343 * x_sf * x_puweight);
    hist49->Fill(s_ptm, 0.0343 * x_sf * x_puweight);
    hist410->Fill(s_mass, 0.0343 * x_sf * x_puweight);
    hist411->Fill(p_mva, 0.0343 * x_sf * x_puweight);  
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist412->Fill((p+j).Pt(), 0.0343 * x_sf * x_puweight); 
  }
  cout<<hist41->GetSumOfWeights()<<endl;
  cout<<hist46->GetSumOfWeights()<<endl;
  cout<<hist410->GetSumOfWeights()<<endl;

  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/2000_N_S1_postproc_WGamma17_full_full_jmcorr_May22.root");
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
    hist51->Fill(p_pt, 0.0343 * x_sf * x_puweight * 1.333);
    hist52->Fill(p_eta, 0.0343 * x_sf * x_puweight * 1.333);
    hist53->Fill(j_pt, 0.0343 * x_sf * x_puweight * 1.333);
    hist54->Fill(j_eta, 0.0343 * x_sf * x_puweight * 1.333);
    hist55->Fill(j_e, 0.0343 * x_sf * x_puweight * 1.333);
    hist56->Fill(j_mass, 400*0.0343 * x_sf * x_puweight * 1.333);
    hist57->Fill(j_tau21, 400*0.0343 * x_sf * x_puweight * 1.333);
    hist58->Fill(s_cos, 0.0343 * x_sf * x_puweight * 1.333);
    hist59->Fill(s_ptm, 0.0343 * x_sf * x_puweight * 1.333);
    hist510->Fill(s_mass, 0.0343 * x_sf * x_puweight * 1.333);
    hist511->Fill(p_mva, 0.0343 * x_sf * x_puweight * 1.333); 
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist512->Fill((p+j).Pt(), 0.0343 * x_sf * x_puweight * 1.333);     
  }
  cout<<hist56->GetSumOfWeights()<<endl;

  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/spin-0_scalar_21Mar23/presel/4000_N_postproc_WGamma17_full_full_jmcorr_Mar23.root");
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
    hist61->Fill(p_pt, 0.0343 * x_sf * x_puweight);
    hist62->Fill(p_eta, 0.0343 * x_sf * x_puweight);
    hist63->Fill(j_pt, 0.0343 * x_sf * x_puweight);
    hist64->Fill(j_eta, 0.0343 * x_sf * x_puweight);
    hist65->Fill(j_e, 0.0343 * x_sf * x_puweight);
    hist66->Fill(j_mass, 400*0.0343 * x_sf * x_puweight);
    hist67->Fill(j_tau21, 400*0.0343 * x_sf * x_puweight);
    hist68->Fill(s_cos, 0.0343 * x_sf * x_puweight);
    hist69->Fill(s_ptm, 0.0343 * x_sf * x_puweight);
    hist610->Fill(s_mass, 0.0343 * x_sf * x_puweight);
    hist611->Fill(p_mva, 0.0343 * x_sf * x_puweight);  
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist612->Fill((p+j).Pt(), 0.0343 * x_sf * x_puweight); 
  }
  
    input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/4000_N_S1_postproc_WGamma17_full_full_jmcorr_May22.root");
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
    hist71->Fill(p_pt, 0.0343 * x_sf * x_puweight);
    hist72->Fill(p_eta, 0.0343 * x_sf * x_puweight);
    hist73->Fill(j_pt, 0.0343 * x_sf * x_puweight);
    hist74->Fill(j_eta, 0.0343 * x_sf * x_puweight);
    hist75->Fill(j_e, 0.0343 * x_sf * x_puweight);
    hist76->Fill(j_mass, 400*0.0343 * x_sf * x_puweight);
    hist77->Fill(j_tau21, 400*0.0343 * x_sf * x_puweight);
    hist78->Fill(s_cos, 0.0343 * x_sf * x_puweight);
    hist79->Fill(s_ptm, 0.0343 * x_sf * x_puweight);
    hist710->Fill(s_mass, 0.0343 * x_sf * x_puweight);
    hist711->Fill(p_mva, 0.0343 * x_sf * x_puweight);  
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist712->Fill((p+j).Pt(), 0.0343 * x_sf * x_puweight); 
  }
  
  /*
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
    hist81->Fill(p_pt, 0.0343 * x_sf * x_puweight * 1.6); //extra 1.6 factor due to lack of finished FULLSIM sample
    hist82->Fill(p_eta, 0.0343 * x_sf * x_puweight * 1.6);
    hist83->Fill(j_pt, 0.0343 * x_sf * x_puweight * 1.6);
    hist84->Fill(j_eta, 0.0343 * x_sf * x_puweight * 1.6);
    hist85->Fill(j_e, 0.0343 * x_sf * x_puweight * 1.6);
    hist86->Fill(j_mass, 400*0.0343 * x_sf * x_puweight * 1.6);
    hist87->Fill(j_tau21, 400*0.0343 * x_sf * x_puweight * 1.6);
    hist88->Fill(s_cos, 0.0343 * x_sf * x_puweight * 1.6);
    hist89->Fill(s_ptm, 0.0343 * x_sf * x_puweight * 1.6);
    hist810->Fill(s_mass, 0.0343 * x_sf * x_puweight * 1.6);
    hist811->Fill(p_mva, 0.0343 * x_sf * x_puweight * 1.6);
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist812->Fill((p+j).Pt(), 0.0343 * x_sf * x_puweight * 1.6);     
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
    hist91->Fill(p_pt, 0.0343 * x_sf * x_puweight * 1.0256);
    hist92->Fill(p_eta, 0.0343 * x_sf * x_puweight * 1.0256);
    hist93->Fill(j_pt, 0.0343 * x_sf * x_puweight * 1.0256);
    hist94->Fill(j_eta, 0.0343 * x_sf * x_puweight * 1.0256);
    hist95->Fill(j_e, 0.0343 * x_sf * x_puweight * 1.0256);
    hist96->Fill(j_mass, 400*0.0343 * x_sf * x_puweight * 1.0256);
    hist97->Fill(j_tau21, 400*0.0343 * x_sf * x_puweight * 1.0256);
    hist98->Fill(s_cos, 0.0343 * x_sf * x_puweight * 1.0256);
    hist99->Fill(s_ptm, 0.0343 * x_sf * x_puweight * 1.0256);
    hist910->Fill(s_mass, 0.0343 * x_sf * x_puweight * 1.0256);
    hist911->Fill(p_mva, 0.0343 * x_sf * x_puweight * 1.0256);  
    TLorentzVector p, j;
    p.SetPtEtaPhiE(p_pt, p_eta, p_phi, p_e);
    j.SetPtEtaPhiE(j_pt, j_eta, j_phi, j_e);
    hist912->Fill((p+j).Pt(), 0.0343 * x_sf * x_puweight * 1.0256);   
  }
  */
  
  
  // Create stack for bkg MC
  THStack *stack1 = new THStack("stack1","pt_{p}");
  THStack *stack2 = new THStack("stack2","eta_{p}");
  THStack *stack3 = new THStack("stack3","pt_{j}");
  THStack *stack4 = new THStack("stack4","eta_{j}");
  THStack *stack5 = new THStack("stack5","Jet softdrop mass");
  THStack *stack6 = new THStack("stack6","Jet tau21");
  THStack *stack7 = new THStack("stack7","cos(#theta*)");
  THStack *stack8 = new THStack("stack8","pt_{j} / M");
  THStack *stack9 = new THStack("stack9","Invariant mass");
  THStack *stack10 = new THStack("stack10","pt_{j#gamma}");
  stack1->Add(hist21); stack1->Add(hist11);
  stack2->Add(hist22); stack2->Add(hist12);
  stack3->Add(hist23); stack3->Add(hist13);
  stack4->Add(hist24); stack4->Add(hist14);
  stack5->Add(hist26); stack5->Add(hist16);
  stack6->Add(hist27); stack6->Add(hist17);
  stack7->Add(hist28); stack7->Add(hist18);
  stack8->Add(hist29); stack8->Add(hist19);
  stack9->Add(hist210); stack9->Add(hist110);
  stack10->Add(hist212); stack10->Add(hist112);


  //=================================================================================
  int color1 = 3;
  int color2 = 7;
  int color3 = kOrange;
  
  //bkg
  hist11->SetFillColor(color1);
  hist12->SetFillColor(color1);
  hist13->SetFillColor(color1);
  hist14->SetFillColor(color1);
  hist15->SetFillColor(color1);
  hist16->SetFillColor(color1);
  hist17->SetFillColor(color1);
  hist18->SetFillColor(color1);
  hist19->SetFillColor(color1);
  hist110->SetFillColor(color1);
  hist111->SetFillColor(color1);
  hist112->SetFillColor(color1);
  hist21->SetFillColor(color2);
  hist22->SetFillColor(color2);
  hist23->SetFillColor(color2);
  hist24->SetFillColor(color2);
  hist25->SetFillColor(color2);
  hist26->SetFillColor(color2);
  hist27->SetFillColor(color2);
  hist28->SetFillColor(color2);
  hist29->SetFillColor(color2);
  hist210->SetFillColor(color2);
  hist211->SetFillColor(color2);
  hist212->SetFillColor(color2);
  
  hist21->SetFillStyle(3144);
  hist22->SetFillStyle(3144);
  hist23->SetFillStyle(3144);
  hist24->SetFillStyle(3144);
  hist25->SetFillStyle(3144);
  hist26->SetFillStyle(3144);
  hist27->SetFillStyle(3144);
  hist28->SetFillStyle(3144);
  hist29->SetFillStyle(3144);
  hist210->SetFillStyle(3144);
  hist211->SetFillStyle(3144);
  hist212->SetFillStyle(3144);
  
  // data
  hist31->SetLineColor(kBlack);
  hist32->SetLineColor(kBlack);
  hist33->SetLineColor(kBlack);
  hist34->SetLineColor(kBlack);
  hist35->SetLineColor(kBlack);
  hist36->SetLineColor(kBlack);
  hist37->SetLineColor(kBlack);
  hist38->SetLineColor(kBlack);
  hist39->SetLineColor(kBlack);
  hist310->SetLineColor(kBlack);
  hist311->SetLineColor(kBlack);
  hist312->SetLineColor(kBlack);
  
  hist31->SetLineWidth(1);
  hist32->SetLineWidth(1);
  hist33->SetLineWidth(1);
  hist34->SetLineWidth(1);
  hist35->SetLineWidth(1);
  hist36->SetLineWidth(1);
  hist37->SetLineWidth(1);
  hist38->SetLineWidth(1);
  hist39->SetLineWidth(1);
  hist310->SetLineWidth(1);
  hist311->SetLineWidth(1);
  hist312->SetLineWidth(1);
  
  hist31->SetMarkerStyle(21);
  hist32->SetMarkerStyle(21);
  hist33->SetMarkerStyle(21);
  hist34->SetMarkerStyle(21);
  hist35->SetMarkerStyle(21);
  hist36->SetMarkerStyle(21);
  hist37->SetMarkerStyle(21);
  hist38->SetMarkerStyle(21);
  hist39->SetMarkerStyle(21);
  hist310->SetMarkerStyle(21);
  hist311->SetMarkerStyle(21);
  hist312->SetMarkerStyle(21);
  
  hist31->SetMarkerSize(5);
  hist32->SetMarkerSize(5);
  hist33->SetMarkerSize(5);
  hist34->SetMarkerSize(5);
  hist35->SetMarkerSize(5);
  hist36->SetMarkerSize(5);
  hist37->SetMarkerSize(5);
  hist38->SetMarkerSize(5);
  hist39->SetMarkerSize(5);
  hist310->SetMarkerSize(5);
  hist311->SetMarkerSize(5);
  hist312->SetMarkerSize(5);
  
  hist31->SetMarkerColor(kBlack);
  hist32->SetMarkerColor(kBlack);
  hist33->SetMarkerColor(kBlack);
  hist34->SetMarkerColor(kBlack);
  hist35->SetMarkerColor(kBlack);
  hist36->SetMarkerColor(kBlack);
  hist37->SetMarkerColor(kBlack);
  hist38->SetMarkerColor(kBlack);
  hist39->SetMarkerColor(kBlack);
  hist310->SetMarkerColor(kBlack);
  hist311->SetMarkerColor(kBlack);
  hist312->SetMarkerColor(kBlack);
  
  hist3a1->SetLineColor(kRed);
  hist3a2->SetLineColor(kRed);
  hist3a3->SetLineColor(kRed);
  hist3a4->SetLineColor(kRed);
  hist3a5->SetLineColor(kRed);
  hist3a6->SetLineColor(kRed);
  hist3a7->SetLineColor(kRed);
  hist3a8->SetLineColor(kRed);
  hist3a9->SetLineColor(kRed);
  hist3a10->SetLineColor(kRed);
  hist3a11->SetLineColor(kRed);
  hist3a12->SetLineColor(kRed);
  
  hist3a1->SetLineWidth(1);
  hist3a2->SetLineWidth(1);
  hist3a3->SetLineWidth(1);
  hist3a4->SetLineWidth(1);
  hist3a5->SetLineWidth(1);
  hist3a6->SetLineWidth(1);
  hist3a7->SetLineWidth(1);
  hist3a8->SetLineWidth(1);
  hist3a9->SetLineWidth(1);
  hist3a10->SetLineWidth(1);
  hist3a11->SetLineWidth(1);
  hist3a12->SetLineWidth(1);
  
  hist3a1->SetMarkerStyle(20);
  hist3a2->SetMarkerStyle(20);
  hist3a3->SetMarkerStyle(20);
  hist3a4->SetMarkerStyle(20);
  hist3a5->SetMarkerStyle(20);
  hist3a6->SetMarkerStyle(20);
  hist3a7->SetMarkerStyle(20);
  hist3a8->SetMarkerStyle(20);
  hist3a9->SetMarkerStyle(20);
  hist3a10->SetMarkerStyle(20);
  hist3a11->SetMarkerStyle(20);
  hist3a12->SetMarkerStyle(20);
  
  hist3a1->SetMarkerSize(5);
  hist3a2->SetMarkerSize(5);
  hist3a3->SetMarkerSize(5);
  hist3a4->SetMarkerSize(5);
  hist3a5->SetMarkerSize(5);
  hist3a6->SetMarkerSize(5);
  hist3a7->SetMarkerSize(5);
  hist3a8->SetMarkerSize(5);
  hist3a9->SetMarkerSize(5);
  hist3a10->SetMarkerSize(5);
  hist3a11->SetMarkerSize(5);
  hist3a12->SetMarkerSize(5);
  
  hist3a1->SetMarkerColor(kRed);
  hist3a2->SetMarkerColor(kRed);
  hist3a3->SetMarkerColor(kRed);
  hist3a4->SetMarkerColor(kRed);
  hist3a5->SetMarkerColor(kRed);
  hist3a6->SetMarkerColor(kRed);
  hist3a7->SetMarkerColor(kRed);
  hist3a8->SetMarkerColor(kRed);
  hist3a9->SetMarkerColor(kRed);
  hist3a10->SetMarkerColor(kRed);
  hist3a11->SetMarkerColor(kRed);
  hist3a12->SetMarkerColor(kRed);
  
  // signal
  int sig_col_1 = kOrange+6;
  int sig_col_2 = kMagenta+1;
  hist41->SetLineColor(sig_col_1);
  hist42->SetLineColor(sig_col_1);
  hist43->SetLineColor(sig_col_1);
  hist44->SetLineColor(sig_col_1);
  hist45->SetLineColor(sig_col_1);
  hist46->SetLineColor(sig_col_1);
  hist47->SetLineColor(sig_col_1);
  hist48->SetLineColor(sig_col_1);
  hist49->SetLineColor(sig_col_1);
  hist410->SetLineColor(sig_col_1);
  hist411->SetLineColor(sig_col_1);
  hist412->SetLineColor(sig_col_1);
  hist51->SetLineColor(sig_col_2);
  hist52->SetLineColor(sig_col_2);
  hist53->SetLineColor(sig_col_2);
  hist54->SetLineColor(sig_col_2);
  hist55->SetLineColor(sig_col_2);
  hist56->SetLineColor(sig_col_2);
  hist57->SetLineColor(sig_col_2);
  hist58->SetLineColor(sig_col_2);
  hist59->SetLineColor(sig_col_2);
  hist510->SetLineColor(sig_col_2);
  hist511->SetLineColor(sig_col_2);
  hist512->SetLineColor(sig_col_2);
  hist61->SetLineColor(sig_col_1);
  hist62->SetLineColor(sig_col_1);
  hist63->SetLineColor(sig_col_1);
  hist64->SetLineColor(sig_col_1);
  hist65->SetLineColor(sig_col_1);
  hist66->SetLineColor(sig_col_1);
  hist67->SetLineColor(sig_col_1);
  hist68->SetLineColor(sig_col_1);
  hist69->SetLineColor(sig_col_1);
  hist610->SetLineColor(sig_col_1);
  hist611->SetLineColor(sig_col_1);
  hist612->SetLineColor(sig_col_1);
  hist71->SetLineColor(sig_col_2);
  hist72->SetLineColor(sig_col_2);
  hist73->SetLineColor(sig_col_2);
  hist74->SetLineColor(sig_col_2);
  hist75->SetLineColor(sig_col_2);
  hist76->SetLineColor(sig_col_2);
  hist77->SetLineColor(sig_col_2);
  hist78->SetLineColor(sig_col_2);
  hist79->SetLineColor(sig_col_2);
  hist710->SetLineColor(sig_col_2);
  hist711->SetLineColor(sig_col_2);
  hist712->SetLineColor(sig_col_2);
  hist81->SetLineColor(kOrange+10);
  hist82->SetLineColor(kOrange+10);
  hist83->SetLineColor(kOrange+10);
  hist84->SetLineColor(kOrange+10);
  hist85->SetLineColor(kOrange+10);
  hist86->SetLineColor(kOrange+10);
  hist87->SetLineColor(kOrange+10);
  hist88->SetLineColor(kOrange+10);
  hist89->SetLineColor(kOrange+10);
  hist810->SetLineColor(kOrange+10);
  hist811->SetLineColor(kOrange+10);
  hist812->SetLineColor(kOrange+10);
  hist91->SetLineColor(kOrange+10);
  hist92->SetLineColor(kOrange+10);
  hist93->SetLineColor(kOrange+10);
  hist94->SetLineColor(kOrange+10);
  hist95->SetLineColor(kOrange+10);
  hist96->SetLineColor(kOrange+10);
  hist97->SetLineColor(kOrange+10);
  hist98->SetLineColor(kOrange+10);
  hist99->SetLineColor(kOrange+10);
  hist910->SetLineColor(kOrange+10);
  hist911->SetLineColor(kOrange+10);
  hist912->SetLineColor(kOrange+10);

  hist41->SetLineStyle(1);
  hist42->SetLineStyle(1);
  hist43->SetLineStyle(1);
  hist44->SetLineStyle(1);
  hist45->SetLineStyle(1);
  hist46->SetLineStyle(1);
  hist47->SetLineStyle(1);
  hist48->SetLineStyle(1);
  hist49->SetLineStyle(1);
  hist410->SetLineStyle(1);
  hist411->SetLineStyle(1);
  hist412->SetLineStyle(1);
  hist51->SetLineStyle(3);
  hist52->SetLineStyle(3);
  hist53->SetLineStyle(3);
  hist54->SetLineStyle(3);
  hist55->SetLineStyle(3);
  hist56->SetLineStyle(3);
  hist57->SetLineStyle(3);
  hist58->SetLineStyle(3);
  hist59->SetLineStyle(3);
  hist510->SetLineStyle(3);
  hist511->SetLineStyle(3);
  hist512->SetLineStyle(3);
  hist61->SetLineStyle(2);
  hist62->SetLineStyle(2);
  hist63->SetLineStyle(2);
  hist64->SetLineStyle(2);
  hist65->SetLineStyle(2);
  hist66->SetLineStyle(2);
  hist67->SetLineStyle(2);
  hist68->SetLineStyle(2);
  hist69->SetLineStyle(2);
  hist610->SetLineStyle(2);
  hist611->SetLineStyle(2);
  hist612->SetLineStyle(2);
  hist71->SetLineStyle(2);
  hist72->SetLineStyle(2);
  hist73->SetLineStyle(2);
  hist74->SetLineStyle(2);
  hist75->SetLineStyle(2);
  hist76->SetLineStyle(2);
  hist77->SetLineStyle(2);
  hist78->SetLineStyle(2);
  hist79->SetLineStyle(2);
  hist710->SetLineStyle(2);
  hist711->SetLineStyle(2);
  hist712->SetLineStyle(2);
  hist81->SetLineStyle(1);
  hist82->SetLineStyle(1);
  hist83->SetLineStyle(1);
  hist84->SetLineStyle(1);
  hist85->SetLineStyle(1);
  hist86->SetLineStyle(1);
  hist87->SetLineStyle(1);
  hist88->SetLineStyle(1);
  hist89->SetLineStyle(1);
  hist810->SetLineStyle(1);
  hist811->SetLineStyle(1);
  hist812->SetLineStyle(1);
  hist91->SetLineStyle(3);
  hist92->SetLineStyle(3);
  hist93->SetLineStyle(3);
  hist94->SetLineStyle(3);
  hist95->SetLineStyle(3);
  hist96->SetLineStyle(3);
  hist97->SetLineStyle(3);
  hist98->SetLineStyle(3);
  hist99->SetLineStyle(3);
  hist910->SetLineStyle(3);
  hist911->SetLineStyle(3);
  hist912->SetLineStyle(3);
  
  hist41->SetLineWidth(2);
  hist42->SetLineWidth(2);
  hist43->SetLineWidth(2);
  hist44->SetLineWidth(2);
  hist45->SetLineWidth(2);
  hist46->SetLineWidth(2);
  hist47->SetLineWidth(2);
  hist48->SetLineWidth(2);
  hist49->SetLineWidth(2);
  hist410->SetLineWidth(2);
  hist411->SetLineWidth(2);
  hist412->SetLineWidth(2);
  hist51->SetLineWidth(3);
  hist52->SetLineWidth(3);
  hist53->SetLineWidth(3);
  hist54->SetLineWidth(3);
  hist55->SetLineWidth(3);
  hist56->SetLineWidth(3);
  hist57->SetLineWidth(3);
  hist58->SetLineWidth(3);
  hist59->SetLineWidth(3);
  hist510->SetLineWidth(3);
  hist511->SetLineWidth(3);
  hist512->SetLineWidth(3);
  hist61->SetLineWidth(2);
  hist62->SetLineWidth(2);
  hist63->SetLineWidth(2);
  hist64->SetLineWidth(2);
  hist65->SetLineWidth(2);
  hist66->SetLineWidth(2);
  hist67->SetLineWidth(2);
  hist68->SetLineWidth(2);
  hist69->SetLineWidth(2);
  hist610->SetLineWidth(2);
  hist611->SetLineWidth(2);
  hist612->SetLineWidth(2);
  hist71->SetLineWidth(3);
  hist72->SetLineWidth(3);
  hist73->SetLineWidth(3);
  hist74->SetLineWidth(3);
  hist75->SetLineWidth(3);
  hist76->SetLineWidth(3);
  hist77->SetLineWidth(3);
  hist78->SetLineWidth(3);
  hist79->SetLineWidth(3);
  hist710->SetLineWidth(3);
  hist711->SetLineWidth(3);
  hist712->SetLineWidth(3);
  hist81->SetLineWidth(2);
  hist82->SetLineWidth(2);
  hist83->SetLineWidth(2);
  hist84->SetLineWidth(2);
  hist85->SetLineWidth(2);
  hist86->SetLineWidth(2);
  hist87->SetLineWidth(2);
  hist88->SetLineWidth(2);
  hist89->SetLineWidth(2);
  hist810->SetLineWidth(2);
  hist811->SetLineWidth(2);
  hist812->SetLineWidth(2);
  hist91->SetLineWidth(3);
  hist92->SetLineWidth(3);
  hist93->SetLineWidth(3);
  hist94->SetLineWidth(3);
  hist95->SetLineWidth(3);
  hist96->SetLineWidth(3);
  hist97->SetLineWidth(3);
  hist98->SetLineWidth(3);
  hist99->SetLineWidth(3);
  hist910->SetLineWidth(3);
  hist911->SetLineWidth(3);
  hist912->SetLineWidth(3);


  TLegend *legend = new TLegend(0.3,0.72,0.92,0.87);
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
  TCanvas *c01 = new TCanvas("c01","",5000,5000);
  c01->cd();
  TPad *p01a = new TPad("p01a","p01a",0,lp3,1,1.0);
  TPad *p01b = new TPad("p01b","p01b",0,lp2,1,hp2);
  TPad *p01c = new TPad("p01c","p01c",0,lp1,1,hp1);
  p01a->Draw();
  p01b->Draw();
  p01c->Draw();
  p01a->cd();
  p01a->SetBottomMargin(main_pad_lower_margin);
  p01a->SetLeftMargin(0.145);
  p01a->SetRightMargin(0.07);
  p01a->SetLogy();
  xaxis1 = hist31->GetXaxis();
  yaxis1 = hist31->GetYaxis();
  xaxis1->SetTitle("p^{T}_{#gamma} [GeV]");
  yaxis1->SetTitle("Events / 40 GeV");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0.2,1000000);
  hist31->Draw("HIST");
  hist3a1->Draw("SAMEHIST");
  stack1->Draw("SAMEHIST");
  hist41->Draw("SAMEHIST");
  hist51->Draw("SAMEHIST");
  hist61->Draw("SAMEHIST");
  hist71->Draw("SAMEHIST");
  // hist81->Draw("SAMEHIST");
  // hist91->Draw("SAMEHIST");
  hist31->Draw("E1SAME");
  hist3a1->Draw("E1SAME");
  hist31->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");


  legend->Draw();
  CMS_lumi(p01a,iPeriod,iPos);
  
  p01b->cd();
  p01b->SetTopMargin(lower_pad1_upper_margin);
  p01b->SetBottomMargin(lower_pad1_lower_margin);
  p01b->SetLeftMargin(0.144);
  p01b->SetRightMargin(0.07);
  pull = (TH1*)hist31->Clone();
  bkg = (TH1*)hist11->Clone();
  bkg->Add((TH1*)hist21->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("p^{T}_{#gamma} [GeV]");
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p01b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p01c->cd();
  p01c->SetTopMargin(lower_pad2_upper_margin);
  p01c->SetBottomMargin(lower_pad2_lower_margin);
  p01c->SetLeftMargin(0.144);
  p01c->SetRightMargin(0.07);
  pull = (TH1*)hist3a1->Clone();
  bkg = (TH1*)hist11->Clone();
  bkg->Add((TH1*)hist21->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("p^{T}_{#gamma} [GeV]");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p01c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c01->Print("p_pt.pdf");
  c01->Print("p_pt.svg");
  //=========================================================
  
  //===========================================================
  TCanvas *c02 = new TCanvas("c02","",5000,5000);
  c02->cd();
  TPad *p02a = new TPad("p02a","p02a",0,lp3,1,1.0);
  TPad *p02b = new TPad("p02b","p02b",0,lp2,1,hp2);
  TPad *p02c = new TPad("p02c","p02c",0,lp1,1,hp1);
  p02a->Draw();
  p02b->Draw();
  p02c->Draw();
  p02a->cd();
  p02a->SetBottomMargin(main_pad_lower_margin);
  p02a->SetLeftMargin(0.145);
  p02a->SetRightMargin(0.07);
  p02a->SetLogy();
  xaxis1 = hist32->GetXaxis();
  yaxis1 = hist32->GetYaxis();
  xaxis1->SetTitle("#eta_{#gamma}");
  yaxis1->SetTitle("Events / 0.08");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0.2,1000000);
  hist32->Draw("HIST");
  hist3a2->Draw("SAMEHIST");
  stack2->Draw("SAMEHIST");
  hist42->Draw("SAMEHIST");
  hist52->Draw("SAMEHIST");
  hist62->Draw("SAMEHIST");
  hist72->Draw("SAMEHIST");
  // hist82->Draw("SAMEHIST");
  // hist92->Draw("SAMEHIST");
  hist32->Draw("E1SAME");
  hist3a2->Draw("E1SAME");
  hist32->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");

  legend->Draw();
  CMS_lumi(p02a,iPeriod,iPos);
  
  p02b->cd();
  p02b->SetTopMargin(lower_pad1_upper_margin);
  p02b->SetBottomMargin(lower_pad1_lower_margin);
  p02b->SetLeftMargin(0.144);
  p02b->SetRightMargin(0.07);
  pull = (TH1*)hist32->Clone();
  bkg = (TH1*)hist12->Clone();
  bkg->Add((TH1*)hist22->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("#eta_{#gamma}");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p02b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p02c->cd();
  p02c->SetTopMargin(lower_pad2_upper_margin);
  p02c->SetBottomMargin(lower_pad2_lower_margin);
  p02c->SetLeftMargin(0.144);
  p02c->SetRightMargin(0.07);
  pull = (TH1*)hist3a2->Clone();
  bkg = (TH1*)hist12->Clone();
  bkg->Add((TH1*)hist22->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("#eta_{#gamma}");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p02c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c02->Print("p_eta.pdf");
  c02->Print("p_eta.svg");
  //=========================================================
  
    //===========================================================
  TCanvas *c03 = new TCanvas("c03","",5000,5000);
  c03->cd();
  TPad *p03a = new TPad("p03a","p03a",0,lp3,1,1.0);
  TPad *p03b = new TPad("p03b","p03b",0,lp2,1,hp2);
  TPad *p03c = new TPad("p03c","p03c",0,lp1,1,hp1);
  p03a->Draw();
  p03b->Draw();
  p03c->Draw();
  p03a->cd();
  p03a->SetBottomMargin(main_pad_lower_margin);
  p03a->SetLeftMargin(0.145);
  p03a->SetRightMargin(0.07);
  p03a->SetLogy();
  xaxis1 = hist33->GetXaxis();
  yaxis1 = hist33->GetYaxis();
  xaxis1->SetTitle("p_{T J} [GeV]");
  yaxis1->SetTitle("Events / 40 Gev");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0.2,1000000);
  hist33->Draw("HIST");
  hist3a3->Draw("SAMEHIST");
  stack3->Draw("SAMEHIST");
  hist43->Draw("SAMEHIST");
  hist53->Draw("SAMEHIST");
  hist63->Draw("SAMEHIST");
  hist73->Draw("SAMEHIST");
  // hist83->Draw("SAMEHIST");
  // hist93->Draw("SAMEHIST");
  hist33->Draw("E1SAME");
  hist3a3->Draw("E1SAME");
  hist33->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");

  legend->Draw();
  CMS_lumi(p03a,iPeriod,iPos);
  
  p03b->cd();
  p03b->SetTopMargin(lower_pad1_upper_margin);
  p03b->SetBottomMargin(lower_pad1_lower_margin);
  p03b->SetLeftMargin(0.144);
  p03b->SetRightMargin(0.07);
  pull = (TH1*)hist33->Clone();
  bkg = (TH1*)hist13->Clone();
  bkg->Add((TH1*)hist23->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("p_{T J} [GeV]");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p03b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p03c->cd();
  p03c->SetTopMargin(lower_pad2_upper_margin);
  p03c->SetBottomMargin(lower_pad2_lower_margin);
  p03c->SetLeftMargin(0.144);
  p03c->SetRightMargin(0.07);
  pull = (TH1*)hist3a3->Clone();
  bkg = (TH1*)hist13->Clone();
  bkg->Add((TH1*)hist23->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("p_{T J} [GeV]");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p03c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c03->Print("j_pt.pdf");
  c03->Print("j_pt.svg");
  //=========================================================
  
  //===========================================================
  TCanvas *c04 = new TCanvas("c04","",5000,5000);
  c04->cd();
  TPad *p04a = new TPad("p04a","p04a",0,lp3,1,1.0);
  TPad *p04b = new TPad("p04b","p04b",0,lp2,1,hp2);
  TPad *p04c = new TPad("p04c","p04c",0,lp1,1,hp1);
  p04a->Draw();
  p04b->Draw();
  p04c->Draw();
  p04a->cd();
  p04a->SetBottomMargin(main_pad_lower_margin);
  p04a->SetLeftMargin(0.145);
  p04a->SetRightMargin(0.07);
  p04a->SetLogy();
  xaxis1 = hist34->GetXaxis();
  yaxis1 = hist34->GetYaxis();
  xaxis1->SetTitle("#eta_{J}");
  yaxis1->SetTitle("Events / 0.1");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0.2,1000000);
  hist34->Draw("PE1");
  hist3a4->Draw("SAMEPE1");
  stack4->Draw("SAMEHIST");
  hist44->Draw("SAMEHIST");
  hist54->Draw("SAMEHIST");
  hist64->Draw("SAMEHIST");
  hist74->Draw("SAMEHIST");
  // hist84->Draw("SAMEHIST");
  // hist94->Draw("SAMEHIST");
  hist34->Draw("E1SAME");
  hist3a4->Draw("E1SAME");
  hist34->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");
  legend->Draw();
  CMS_lumi(p04a,iPeriod,iPos);
  
  p04b->cd();
  p04b->SetTopMargin(lower_pad1_upper_margin);
  p04b->SetBottomMargin(lower_pad1_lower_margin);
  p04b->SetLeftMargin(0.144);
  p04b->SetRightMargin(0.07);
  pull = (TH1*)hist34->Clone();
  bkg = (TH1*)hist14->Clone();
  bkg->Add((TH1*)hist24->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("#eta_{J}");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p04b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p04c->cd();
  p04c->SetTopMargin(lower_pad2_upper_margin);
  p04c->SetBottomMargin(lower_pad2_lower_margin);
  p04c->SetLeftMargin(0.144);
  p04c->SetRightMargin(0.07);
  pull = (TH1*)hist3a4->Clone();
  bkg = (TH1*)hist14->Clone();
  bkg->Add((TH1*)hist24->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("#eta_{J}");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p04c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c04->Print("j_eta.pdf");
  c04->Print("j_eta.svg");
  //==========================================================
  
  
  //===========================================================
  TCanvas *c05 = new TCanvas("c05","",5000,5000);
  c05->cd();
  TPad *p05a = new TPad("p05a","p05a",0,lp3,1,1.0);
  TPad *p05b = new TPad("p05b","p05b",0,lp2,1,hp2);
  TPad *p05c = new TPad("p05c","p05c",0,lp1,1,hp1);
  p05a->Draw();
  p05b->Draw();
  p05c->Draw();
  p05a->cd();
  p05a->SetBottomMargin(main_pad_lower_margin);
  p05a->SetLeftMargin(0.145);
  p05a->SetRightMargin(0.07);
  xaxis1 = hist36->GetXaxis();
  yaxis1 = hist36->GetYaxis();
  xaxis1->SetTitle("m^{SD}_{J} [GeV]");
  yaxis1->SetTitle("Events / 1 GeV");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0,12000);
  hist36->Draw("PE1");
  hist3a6->Draw("SAMEPE1");
  stack5->Draw("SAMEHIST");
  hist46->Draw("SAMEHIST");
  hist56->Draw("SAMEHIST");
  hist66->Draw("SAMEHIST");
  hist76->Draw("SAMEHIST");
  // hist86->Draw("SAMEHIST");
  // hist96->Draw("SAMEHIST");
  hist36->Draw("E1SAME");
  hist3a6->Draw("E1SAME");
  hist36->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");
  legend->Draw();
  CMS_lumi(p05a,iPeriod,iPos);
  
  p05b->cd();
  p05b->SetTopMargin(lower_pad1_upper_margin);
  p05b->SetBottomMargin(lower_pad1_lower_margin);
  p05b->SetLeftMargin(0.144);
  p05b->SetRightMargin(0.07);
  pull = (TH1*)hist36->Clone();
  bkg = (TH1*)hist16->Clone();
  bkg->Add((TH1*)hist26->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("m^{SD}_{J} [GeV]");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0.6,1.4);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.13);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p05b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p05c->cd();
  p05c->SetTopMargin(lower_pad2_upper_margin);
  p05c->SetBottomMargin(lower_pad2_lower_margin);
  p05c->SetLeftMargin(0.144);
  p05c->SetRightMargin(0.07);
  pull = (TH1*)hist3a6->Clone();
  bkg = (TH1*)hist16->Clone();
  bkg->Add((TH1*)hist26->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("m^{SD}_{J} [GeV]");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0.6,1.4);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.11);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.5);
  p05c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c05->Print("j_m.pdf");
  c05->Print("j_m.svg");
  //==========================================================
  
  //===========================================================
  TCanvas *c06 = new TCanvas("c06","",5000,5000);
  c06->cd();
  TPad *p06a = new TPad("p06a","p06a",0,lp3,1,1.0);
  TPad *p06b = new TPad("p06b","p06b",0,lp2,1,hp2);
  TPad *p06c = new TPad("p06c","p06c",0,lp1,1,hp1);
  p06a->Draw();
  p06b->Draw();
  p06c->Draw();
  p06a->cd();
  p06a->SetBottomMargin(main_pad_lower_margin);
  p06a->SetLeftMargin(0.145);
  p06a->SetRightMargin(0.07);
  //p06a->SetLogy();
  xaxis1 = hist37->GetXaxis();
  yaxis1 = hist37->GetYaxis();
  xaxis1->SetTitle("#tau_{21}");
  yaxis1->SetTitle("Events / 0.02");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0,12000);
  hist37->Draw("PE1");
  hist3a7->Draw("SAMEPE1");
  stack6->Draw("SAMEHIST");
  hist47->Draw("SAMEHIST");
  hist57->Draw("SAMEHIST");
  hist67->Draw("SAMEHIST");
  hist77->Draw("SAMEHIST");
  // hist87->Draw("SAMEHIST");
  // hist97->Draw("SAMEHIST");
  hist37->Draw("E1SAME");
  hist3a7->Draw("E1SAME");
  hist37->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");
  legend->Draw();
  CMS_lumi(p06a,iPeriod,iPos);
  
  p06b->cd();
  p06b->SetTopMargin(lower_pad1_upper_margin);
  p06b->SetBottomMargin(lower_pad1_lower_margin);
  p06b->SetLeftMargin(0.144);
  p06b->SetRightMargin(0.07);
  pull = (TH1*)hist37->Clone();
  bkg = (TH1*)hist17->Clone();
  bkg->Add((TH1*)hist27->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("#tau_{21}");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p06b->SetGrid();
 bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p06c->cd();
  p06c->SetTopMargin(lower_pad2_upper_margin);
  p06c->SetBottomMargin(lower_pad2_lower_margin);
  p06c->SetLeftMargin(0.144);
  p06c->SetRightMargin(0.07);
  pull = (TH1*)hist3a7->Clone();
  bkg = (TH1*)hist17->Clone();
  bkg->Add((TH1*)hist27->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("#tau_{21}");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p06c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c06->Print("j_tau21.pdf");
  c06->Print("j_tau21.svg");
  //==========================================================
  
  //===========================================================
  TCanvas *c07 = new TCanvas("c07","",5000,5000);
  c07->cd();
  TPad *p07a = new TPad("p07a","p07a",0,lp3,1,1.0);
  TPad *p07b = new TPad("p07b","p07b",0,lp2,1,hp2);
  TPad *p07c = new TPad("p07c","p07c",0,lp1,1,hp1);
  p07a->Draw();
  p07b->Draw();
  p07c->Draw();
  p07a->cd();
  p07a->SetBottomMargin(main_pad_lower_margin);
  p07a->SetLeftMargin(0.145);
  p07a->SetRightMargin(0.07);
  p07a->SetLogy();
  xaxis1 = hist38->GetXaxis();
  yaxis1 = hist38->GetYaxis();
  xaxis1->SetTitle("cos(#theta*_{#gamma})");
  yaxis1->SetTitle("Events / 0.02");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0.2,1000000);
  hist38->Draw("PE1");
  hist3a8->Draw("SAMEPE1");
  stack7->Draw("SAMEHIST");
  hist48->Draw("SAMEHIST");
  hist58->Draw("SAMEHIST");
  hist68->Draw("SAMEHIST");
  hist78->Draw("SAMEHIST");
  // hist88->Draw("SAMEHIST");
  // hist98->Draw("SAMEHIST");
  hist38->Draw("E1SAME");
  hist3a8->Draw("E1SAME");
  hist38->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");
  legend->Draw();
  CMS_lumi(p07a,iPeriod,iPos);
  
  p07b->cd();
  p07b->SetTopMargin(lower_pad1_upper_margin);
  p07b->SetBottomMargin(lower_pad1_lower_margin);
  p07b->SetLeftMargin(0.144);
  p07b->SetRightMargin(0.07);
  pull = (TH1*)hist38->Clone();
  bkg = (TH1*)hist18->Clone();
  bkg->Add((TH1*)hist28->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("cos(#theta*_{#gamma})");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p07b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p07c->cd();
  p07c->SetTopMargin(lower_pad2_upper_margin);
  p07c->SetBottomMargin(lower_pad2_lower_margin);
  p07c->SetLeftMargin(0.144);
  p07c->SetRightMargin(0.07);
  pull = (TH1*)hist3a8->Clone();
  bkg = (TH1*)hist18->Clone();
  bkg->Add((TH1*)hist28->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("cos(#theta*_{#gamma})");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p07c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c07->Print("s_cos.pdf");
  c07->Print("s_cos.svg");
  //==========================================================
  
  //===========================================================
  TCanvas *c08 = new TCanvas("c08","",5000,5000);
  c08->cd();
  TPad *p08a = new TPad("p08a","p08a",0,lp3,1,1.0);
  TPad *p08b = new TPad("p08b","p08b",0,lp2,1,hp2);
  TPad *p08c = new TPad("p08c","p08c",0,lp1,1,hp1);
  p08a->Draw();
  p08b->Draw();
  p08c->Draw();
  p08a->cd();
  p08a->SetBottomMargin(main_pad_lower_margin);
  p08a->SetLeftMargin(0.145);
  p08a->SetRightMargin(0.07);
  p08a->SetLogy();
  xaxis1 = hist39->GetXaxis();
  yaxis1 = hist39->GetYaxis();
  xaxis1->SetTitle("pT_{#gamma} / m_{J#gamma}");
  yaxis1->SetTitle("Events / 0.04");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0.2,1000000);
  hist39->Draw("PE1");
  hist3a9->Draw("SAMEPE1");
  stack8->Draw("SAMEHIST");
  hist49->Draw("SAMEHIST");
  hist59->Draw("SAMEHIST");
  hist69->Draw("SAMEHIST");
  hist79->Draw("SAMEHIST");
  // hist89->Draw("SAMEHIST");
  // hist99->Draw("SAMEHIST");
  hist39->Draw("E1SAME");
  hist3a9->Draw("E1SAME");
  hist39->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");
  legend->Draw();
  CMS_lumi(p08a,iPeriod,iPos);
  
  p08b->cd();
  p08b->SetTopMargin(lower_pad1_upper_margin);
  p08b->SetBottomMargin(lower_pad1_lower_margin);
  p08b->SetLeftMargin(0.144);
  p08b->SetRightMargin(0.07);
  pull = (TH1*)hist39->Clone();
  bkg = (TH1*)hist19->Clone();
  bkg->Add((TH1*)hist29->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("p_{T}^{#gamma} / m_{J#gamma}");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p08b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p08c->cd();
  p08c->SetTopMargin(lower_pad2_upper_margin);
  p08c->SetBottomMargin(lower_pad2_lower_margin);
  p08c->SetLeftMargin(0.144);
  p08c->SetRightMargin(0.07);
  pull = (TH1*)hist3a9->Clone();
  bkg = (TH1*)hist19->Clone();
  bkg->Add((TH1*)hist29->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("p_{T}^{#gamma} / m_{J#gamma}");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p08c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c08->Print("s_ptm.pdf");
  c08->Print("s_ptm.svg");
  //==========================================================
  
  //===========================================================
  TCanvas *c09 = new TCanvas("c09","",5000,5000);
  c09->cd();
  TPad *p09a = new TPad("p09a","p09a",0,lp3,1,1.0);
  TPad *p09b = new TPad("p09b","p09b",0,lp2,1,hp2);
  TPad *p09c = new TPad("p09c","p09c",0,lp1,1,hp1);
  p09a->Draw();
  p09b->Draw();
  p09c->Draw();
  p09a->cd();
  p09a->SetBottomMargin(main_pad_lower_margin);
  p09a->SetLeftMargin(0.145);
  p09a->SetRightMargin(0.07);
  p09a->SetLogy();
  xaxis1 = hist310->GetXaxis();
  yaxis1 = hist310->GetYaxis();
  xaxis1->SetTitle("m_{J#gamma}");
  yaxis1->SetTitle("Events / 40 GeV");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0.2,100000);
  hist310->Draw("PE1");
  hist3a10->Draw("SAMEPE1");
  stack9->Draw("SAMEHIST");
  hist410->Draw("SAMEHIST");
  hist510->Draw("SAMEHIST");
  hist610->Draw("SAMEHIST");
  hist710->Draw("SAMEHIST");
  // hist810->Draw("SAMEHIST");
  // hist910->Draw("SAMEHIST");
  hist310->Draw("E1SAME");
  hist3a10->Draw("E1SAME");
  hist310->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist31,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a1,"Data SR","e1p");
  legend->AddEntry(hist11,"#gamma + jet","f");
  legend->AddEntry(hist21,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");
  legend->Draw();
  CMS_lumi(p09a,iPeriod,iPos);
  
  p09b->cd();
  p09b->SetTopMargin(lower_pad1_upper_margin);
  p09b->SetBottomMargin(lower_pad1_lower_margin);
  p09b->SetLeftMargin(0.144);
  p09b->SetRightMargin(0.07);
  pull = (TH1*)hist310->Clone();
  bkg = (TH1*)hist110->Clone();
  bkg->Add((TH1*)hist210->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("m_{J#gamma}");
  xaxis2->SetTitleOffset(1.15);
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p09b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p09c->cd();
  p09c->SetTopMargin(lower_pad2_upper_margin);
  p09c->SetBottomMargin(lower_pad2_lower_margin);
  p09c->SetLeftMargin(0.144);
  p09c->SetRightMargin(0.07);
  pull = (TH1*)hist3a10->Clone();
  bkg = (TH1*)hist110->Clone();
  bkg->Add((TH1*)hist210->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("m_{J#gamma}");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p09c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c09->Print("s_M.pdf");
  c09->Print("s_M.svg");
  //==========================================================
  
  //===========================================================
  TCanvas *c10 = new TCanvas("c10","",5000,5000);
  c10->cd();
  TPad *p10a = new TPad("p10a","p10a",0,lp3,1,1.0);
  TPad *p10b = new TPad("p10b","p10b",0,lp2,1,hp2);
  TPad *p10c = new TPad("p10c","p10c",0,lp1,1,hp1);
  p10a->Draw();
  p10b->Draw();
  p10c->Draw();
  p10a->cd();
  p10a->SetBottomMargin(main_pad_lower_margin);
  p10a->SetLeftMargin(0.145);
  p10a->SetRightMargin(0.07);
  p10a->SetLogy();
  xaxis1 = hist312->GetXaxis();
  yaxis1 = hist312->GetYaxis();
  xaxis1->SetTitle("p_{T j#gamma} [GeV]");
  yaxis1->SetTitle("Events / 30 GeV");
  xaxis1->SetLabelOffset(2);
  yaxis1->SetTitleOffset(1.4);
  yaxis1->SetRangeUser(0.2,1000000);
  hist312->Draw("HIST");
  hist3a12->Draw("SAMEHIST");
  stack10->Draw("SAMEHIST");
  hist412->Draw("SAMEHIST");
  hist512->Draw("SAMEHIST");
  hist612->Draw("SAMEHIST");
  hist712->Draw("SAMEHIST");
  // hist81->Draw("SAMEHIST");
  // hist91->Draw("SAMEHIST");
  hist312->Draw("E1SAME");
  hist3a12->Draw("E1SAME");
  hist312->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist312,"Data SB (Norm)","e1p");
  legend->AddEntry(hist3a12,"Data SR","e1p");
  legend->AddEntry(hist112,"#gamma + jet","f");
  legend->AddEntry(hist212,"QCD","f"); 
  legend->AddEntry(hist41,"Sig. 2.0 TeV, spin 0");
  legend->AddEntry(hist51,"Sig. 2.0 TeV, spin 1");
  legend->AddEntry(hist61,"Sig. 4.0 TeV, spin 0");
  legend->AddEntry(hist71,"Sig. 4.0 TeV, spin 1");
  // legend->AddEntry(hist81,"M-5000 N");
  // legend->AddEntry(hist91,"M-5000 W");


  legend->Draw();
  CMS_lumi(p10a,iPeriod,iPos);
  
  p10b->cd();
  p10b->SetTopMargin(lower_pad1_upper_margin);
  p10b->SetBottomMargin(lower_pad1_lower_margin);
  p10b->SetLeftMargin(0.144);
  p10b->SetRightMargin(0.07);
  pull = (TH1*)hist312->Clone();
  bkg = (TH1*)hist112->Clone();
  bkg->Add((TH1*)hist212->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  bkgerrline->SetLineColor(kBlack);
  bkgerrline->SetLineWidth(1);
  pull->Divide(bkg);
  xaxis2 = bkgerr->GetXaxis();
  yaxis2 = bkgerr->GetYaxis();
  xaxis2->SetTitle("p_{T j#gamma} [GeV]");
  yaxis2->SetTitle("Data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.15);
  xaxis2->SetLabelOffset(1);
  xaxis2->SetTitleSize(0.15);
  yaxis2->SetLabelSize(0.15);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.15);
  yaxis2->SetTitleOffset(0.4);
  p10b->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  
  p10c->cd();
  p10c->SetTopMargin(lower_pad2_upper_margin);
  p10c->SetBottomMargin(lower_pad2_lower_margin);
  p10c->SetLeftMargin(0.144);
  p10c->SetRightMargin(0.07);
  pull = (TH1*)hist3a12->Clone();
  bkg = (TH1*)hist112->Clone();
  bkg->Add((TH1*)hist212->Clone());
  bkgerr = (TH1*)bkg->Clone();
  bkgerr->Divide((TH1*)bkg->Clone());
  for(int i=0; i<bkgerr->GetNbinsX(); i++){
      if(bkg->GetBinContent(i+1) > 0){
        bkgerr->SetBinError(i+1,bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1) > 0.97 ? 0.97 : bkgerr->GetBinError(i+1)/bkgerr->GetBinContent(i+1));
        
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
  xaxis3 = bkgerr->GetXaxis();
  yaxis3 = bkgerr->GetYaxis();
  xaxis3->SetTitle("p_{T j#gamma} [GeV]");
  xaxis3->SetTitleOffset(1.1);
  yaxis3->SetTitle("Data/MC");
  yaxis3->SetRangeUser(0,2);
  xaxis3->SetLabelSize(0.15);
  xaxis3->SetLabelOffset(0.05);
  xaxis3->SetTitleSize(0.18);
  yaxis3->SetLabelSize(0.13);
  yaxis3->SetNdivisions(5);
  yaxis3->SetTitleSize(0.13);
  yaxis3->SetTitleOffset(0.46);
  p10c->SetGrid();
  bkgerr->Draw("E2");
  bkgerrline->Draw("hist SAME");
  pull->Draw("PE1 SAME");
  c10->Print("pt_jgamma.pdf");
  c10->Print("pt_jgamma.svg");
  //=========================================================
  
  cout<<SB<<endl;
  cout<<WB<<endl;
}