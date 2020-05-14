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

void postproc(std::string dataset, int runondata, int runyear)
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
  
  //TString dataset = "SignalMC"+std::to_string(mass)+"N";
  
  TString year_str = std::to_string(runyear);
  
 TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/20"+year_str+"/ntuple/ntuple_data/"+dataset+"_nominal_pileup_WGamma_full_full_Mar17.root");
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
  //TFile *outFile = TFile::Open(dataset+"_postproc_WGamma"+year_str+"_full_full_presel_jmcorr_kfactor_Mar17.root", "RECREATE");
  TFile *outFile = TFile::Open(dataset+"_postproc_WGamma"+year_str+"_full_full_presel_jmcorr_Mar17.root", "RECREATE");
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
  
  TFile* pileup_central = TFile::Open("/afs/cern.ch/user/x/xuyan/WGProj/PROD17/DATA/pileup/Pileup_"+year_str+".root", "READ");
  TH1F* pileup_Data_central = (TH1F*)pileup_central->Get("pileup");
  pileup_Data_central->Scale(1/(double)pileup_Data_central->Integral());
  pileup_Data_central->Divide(pileup_MC);
  
  //jet mass correction
  //16
  TF1* f1 = new TF1("f1","[0]/pow((x+[5])/[6],[1])+[2]/pow((x+[5])/[6]*100,[3])+[4]",400,3000);
  f1->SetParameters(3.610,3.840,-458.30,0.4894,126.0119,1783.41,3403.76);
  TF1* f1a = new TF1("f1a","pol1",200,400);
  f1a->SetParameters(77.9253,0.019194);
  
  //17
  TF1* f2 = new TF1("f2","[0]/pow((x+[5])/[6],[1])+[2]/pow((x+[5])/[6]*100,[3])+[4]",400,3000);
  f2->SetParameters(2.2210,14.0254,-438.18,0.5065,119.99,4428.19,5359.17);
  TF1* f2a = new TF1("f2a","pol1",200,400);
  f2a->SetParameters(76.332,0.020022);
  
   
  float sumW = 0;
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
	
	//Apply mass correction
	double masscorr = 1;
	
	if(runyear == 16){
	
	if(j_pt < 400)
			masscorr = 80.379 / f1a->Eval(j_pt);
	else
			masscorr = 80.379 / f1->Eval(j_pt);
	}
	else{
		
	if(j_pt < 400)
			masscorr = 80.379 / f2a->Eval(j_pt);
	else
			masscorr = 80.379 / f2->Eval(j_pt);
	}
	if(masscorr > 5)
		cout<<"ERROR IN MASS CORR CALCULATION"<<endl;
	
	//if(s_mass < 0.75*mass || s_mass > 1.25*mass) continue;
    //if(j_mass * masscorr < 68 || j_mass * masscorr > 94) continue;
	//if(j_mass * masscorr < 38 || j_mass * masscorr > 64) continue;
	//if(j_mass * masscorr < 88 || j_mass * masscorr > 114) continue;
	//if(j_mass * masscorr < 40 || j_mass * masscorr > 65) continue;
	if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    // if(abs(p_eta) > 1.44) continue;
    // if(abs(j_eta) > 2) continue;
    // if(j_tau21 > 0.35) continue;
    // if(s_ptm < 0.37) continue;
    // if(s_cos > 0.6) continue;
    
    photon_pt = p_pt;
    photon_eta = p_eta;
    photon_phi = p_phi;
    photon_e = p_e;
    ak8puppijet_pt = j_pt;
    ak8puppijet_eta = j_eta;
    ak8puppijet_phi = j_phi;
    ak8puppijet_e = j_e;
	ak8puppijet_masssoftdropcorr = j_mass * masscorr;
    ak8puppijet_tau21 = j_tau21;
    sys_costhetastar = s_cos;
    sys_ptoverm = s_ptm;
    m = s_mass;
    xsec_weight = x_weight;
	sys_pvn = s_pv;
	std::string str = dataset;
	if(dataset.find("GJets") != std::string::npos)
		xsec_kfactor = 1.21;
	else if(dataset.find("QCD") != std::string::npos)
		xsec_kfactor = 1.03;
	else
		xsec_kfactor = 1;
	
	//add pileup-weight
	int binnum = pileup_Data_central->GetXaxis()->FindBin(s_pv);
	double pileupweight_central = pileup_Data_central->GetBinContent(binnum);
	if(runondata)
		xsec_puweight  = 1;
	else
		xsec_puweight = pileupweight_central ;

	outTree->Fill();
	sumW += xsec_weight * xsec_kfactor * xsec_puweight;
  }
  cout<<"Entries: "<<outTree->GetEntries()<<endl;
  // cout<<"++++,"<<mass<<","<<sumW<<endl;
  outFile->Write();
  outFile->Close();
}