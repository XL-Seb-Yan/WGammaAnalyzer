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

void analyzer(const TString conf="samples.conf"){

  gBenchmark->Start("analyzerWG");
  gSystem->Load("lib/libMylib.so");
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //=============================================================================================================
  TString type = "CONTROL"; //DATA, SIGNAL, BKG, CONTROL
  TString formate = "DIVIDED"; //DIVIDED RANGE
  TString histotitle;
  if(type == "SIGNAL"){
    histotitle = "MC(2017) W#gamma";
  }
  else if(type == "BKG"){
    histotitle = "MC(2017) QCD Background";
  }
  else if(type == "DATA"){
    histotitle = "DATA(2017 Single Photon)";
  }
  else if(type == "CONTROL"){
    histotitle = "DATA (2017 Single Muon)";
  }
  else
    cout<<"Wrong data type"<<endl;


  TH1F* hist01 = new TH1F("WGamma01",histotitle+", S/sqrt(B)",50,0,1); //Pass: tight ID and EleVeto
  //Photons
  /*
  TH1F* hist01pa = new TH1F("WGamma01pa",histotitle+", pt_{#gamma}",60,0,1200); //Pass: tight ID and EleVeto
  TH1F* hist02pa = new TH1F("WGamma02pa",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pa = new TH1F("WGamma03pa",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pa = new TH1F("WGamma04pa",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pa = new TH1F("WGamma05pa",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pa = new TH1F("WGamma06pa",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pa = new TH1F("WGamma07pa",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pa = new TH1F("WGamma08pa",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pa = new TH1F("WGamma09pa",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pa = new TH1F("WGamma10pa",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pa = new TH1F("WGamma11pa",histotitle+", Photon mvaID Category",10,-1,3);

  TH1F* hist01pb = new TH1F("WGamma01pb",histotitle+", pt_{#gamma}",60,0,1200); //Fired trigger 300
  TH1F* hist02pb = new TH1F("WGamma02pb",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pb = new TH1F("WGamma03pb",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pb = new TH1F("WGamma04pb",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pb = new TH1F("WGamma05pb",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pb = new TH1F("WGamma06pb",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pb = new TH1F("WGamma07pb",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pb = new TH1F("WGamma08pb",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pb = new TH1F("WGamma09pb",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pb = new TH1F("WGamma10pb",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pb = new TH1F("WGamma11pb",histotitle+", Photon mvaID Category",10,-1,3);

  TH1F* hist01pc = new TH1F("WGamma01pc",histotitle+", pt_{#gamma}",60,0,1200); //Fired trigger 200
  TH1F* hist02pc = new TH1F("WGamma02pc",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pc = new TH1F("WGamma03pc",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pc = new TH1F("WGamma04pc",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pc = new TH1F("WGamma05pc",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pc = new TH1F("WGamma06pc",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pc = new TH1F("WGamma07pc",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pc = new TH1F("WGamma08pc",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pc = new TH1F("WGamma09pc",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pc = new TH1F("WGamma10pc",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pc = new TH1F("WGamma11pc",histotitle+", Photon mvaID Category",10,-1,3);

  TH1F* hist01pd = new TH1F("WGamma01pd",histotitle+", pt_{#gamma}",60,0,1200); //Fired trigger 150
  TH1F* hist02pd = new TH1F("WGamma02pd",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pd = new TH1F("WGamma03pd",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pd = new TH1F("WGamma04pd",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pd = new TH1F("WGamma05pd",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pd = new TH1F("WGamma06pd",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pd = new TH1F("WGamma07pd",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pd = new TH1F("WGamma08pd",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pd = new TH1F("WGamma09pd",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pd = new TH1F("WGamma10pd",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pd = new TH1F("WGamma11pd",histotitle+", Photon mvaID Category",10,-1,3);

  TH1F* hist01pe = new TH1F("WGamma01pe",histotitle+", pt_{#gamma}",60,0,1200); //Fired trigger 165
  TH1F* hist02pe = new TH1F("WGamma02pe",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pe = new TH1F("WGamma03pe",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pe = new TH1F("WGamma04pe",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pe = new TH1F("WGamma05pe",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pe = new TH1F("WGamma06pe",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pe = new TH1F("WGamma07pe",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pe = new TH1F("WGamma08pe",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pe = new TH1F("WGamma09pe",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pe = new TH1F("WGamma10pe",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pe = new TH1F("WGamma11pe",histotitle+", Photon mvaID Category",10,-1,3);
  */
  
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;
  gStyle->SetOptStat(0);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  // parse .conf file
  confParse(conf, snamev, samplev);
  
  // Data structures to store info from produced flat ntuples
  float photon_pt;
  float photon_eta;
  float photon_phi;
  float photon_e;
  float photon_mvaval;
  float photon_mvacat;
  float ak8puppijet_pt;
  float ak8puppijet_eta;
  float ak8puppijet_phi;
  float ak8puppijet_e;
  float ak8puppijet_masssoftdropcorr;
  float ak8puppijet_tau21;
  float sys_costhetastar;
  float sys_ptoverm;
  float sys_invmass;


  // S/sqrt(B) data
  std::vector<int> sig_N, bkg_N;
  std::vector<float> recordarr;
  
  // loop over samples
  TTree* eventTree = 0;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    CSample* samp = samplev[isam];
    cout<<"begin loop over files"<<endl;
    TStopwatch stopwatch;

    // loop through files
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {  
      
      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile]<<endl; cout.flush();
      TFile *infile = TFile::Open(samp->fnamev[ifile]);
      assert(infile);

      // Access Event Tree
      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);
      eventTree->SetBranchAddress("photon_pt", &photon_pt);                             TBranch *photonPtBr = eventTree->GetBranch("photon_pt");
      eventTree->SetBranchAddress("photon_eta", &photon_eta);                           TBranch *photonEtaBr = eventTree->GetBranch("photon_eta");
      eventTree->SetBranchAddress("photon_phi", &photon_phi);                           TBranch *photonPhiBr = eventTree->GetBranch("photon_phi");
      eventTree->SetBranchAddress("photon_e", &photon_e);                               TBranch *photonEBr = eventTree->GetBranch("photon_e");
      eventTree->SetBranchAddress("photon_mvaval", &photon_mvaval);                     TBranch *photonMvaValBr = eventTree->GetBranch("photon_mvaval");
      eventTree->SetBranchAddress("photon_mvacat", &photon_mvacat);                     TBranch *photonMvaCatBr = eventTree->GetBranch("photon_mvacat");
      eventTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                             TBranch *jetPtBr = eventTree->GetBranch("ak8puppijet_pt");
      eventTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                           TBranch *jetEtaBr = eventTree->GetBranch("ak8puppijet_eta");
      eventTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                           TBranch *jetPhiBr = eventTree->GetBranch("ak8puppijet_phi");
      eventTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               TBranch *jetEBr = eventTree->GetBranch("ak8puppijet_e");
      eventTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr); TBranch *jetMassBr = eventTree->GetBranch("ak8puppijet_masssoftdropcorr");
      eventTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                       TBranch *jetTau21Br = eventTree->GetBranch("ak8puppijet_tau21");
      eventTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                         TBranch *sysCosthetastarBr = eventTree->GetBranch("sys_costhetastar");
      eventTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                   TBranch *sysPtovermEtaBr = eventTree->GetBranch("sys_ptoverm");
      eventTree->SetBranchAddress("sys_invmass", &sys_invmass);                                   TBranch *sysInvmassBr = eventTree->GetBranch("sys_invmass");


      if(isam == 0){
	//TH1F* histS01 = new TH1F("histS01",histotitle+", invmass",60,0,1200);
	for(float cut2 = 0; cut2 < 1; cut2+=0.05){
	  int N = 0;
	  float record = cut2;
	  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	    // Get Events
	    photonPtBr->GetEntry(ientry);
	    photonEtaBr->GetEntry(ientry);
	    photonPhiBr->GetEntry(ientry);
	    photonEBr->GetEntry(ientry);
	    photonMvaValBr->GetEntry(ientry);
	    photonMvaCatBr->GetEntry(ientry);
	    jetPtBr->GetEntry(ientry);
	    jetEtaBr->GetEntry(ientry);
	    jetPhiBr->GetEntry(ientry);
	    jetEBr->GetEntry(ientry);
	    jetMassBr->GetEntry(ientry);
	    jetTau21Br->GetEntry(ientry);
	    sysCosthetastarBr->GetEntry(ientry);
	    sysPtovermEtaBr->GetEntry(ientry);
	    sysInvmassBr->GetEntry(ientry);

	    if(sys_ptoverm > cut2 && (sys_invmass > 1360 && sys_invmass < 1840)){
	      N++;
	    }
	  }
	  sig_N.push_back(N);
	  recordarr.push_back(record);
	}
	
	// Post processing - fitting
	/*
	TF1 *fS01 = new TF1 ("fS01", "gaus", 65, 105);
	fS01->SetParNames("Constant 1", "Mean 1", "Sigma 1"); 
	gStyle->SetOptFit(1); 
	Double_t parS01[3]= {70,78,7};
	fS01->SetParameters(parS01);
	histS01->Fit(fs01);
	*/
      }

      if(isam == 1){
	//TH1F* histS01 = new TH1F("histS01",histotitle+", invmass",60,0,1200);
	for(float cut2 = 0; cut2 < 1; cut2+=0.05){
	  int N = 0;
	  float record = cut2;
	  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	    // Get Events
	    photonPtBr->GetEntry(ientry);
	    photonEtaBr->GetEntry(ientry);
	    photonPhiBr->GetEntry(ientry);
	    photonEBr->GetEntry(ientry);
	    photonMvaValBr->GetEntry(ientry);
	    photonMvaCatBr->GetEntry(ientry);
	    jetPtBr->GetEntry(ientry);
	    jetEtaBr->GetEntry(ientry);
	    jetPhiBr->GetEntry(ientry);
	    jetEBr->GetEntry(ientry);
	    jetMassBr->GetEntry(ientry);
	    jetTau21Br->GetEntry(ientry);
	    sysCosthetastarBr->GetEntry(ientry);
	    sysPtovermEtaBr->GetEntry(ientry);
	    sysInvmassBr->GetEntry(ientry);

	    if(sys_ptoverm > cut2 && (sys_invmass > 1360 && sys_invmass < 1840)){
	      N++;
	    }
	  }
	  bkg_N.push_back(N);
	  recordarr.push_back(record);
	}
	
	// Post processing - fitting
	/*
	TF1 *fS01 = new TF1 ("fS01", "gaus", 65, 105);
	fS01->SetParNames("Constant 1", "Mean 1", "Sigma 1"); 
	gStyle->SetOptFit(1); 
	Double_t parS01[3]= {70,78,7};
	fS01->SetParameters(parS01);
	histS01->Fit(fs01);
	*/
      }
      cout<<double(ifile)/double(nfiles)*100<<" % done with this dataset"<<endl;
      Double_t elapsed_t_file = stopwatch.RealTime() / (ifile+1);
      stopwatch.Start(kFALSE);
      Int_t time_remain = int(elapsed_t_file*(nfiles-ifile-1));
      Int_t hours = 0, minutes = 0, seconds = 0;
      if(time_remain < 60)
	seconds = time_remain;
      else if(time_remain < 3600){
	minutes = time_remain / 60;
	seconds = time_remain % 60;
      }
      else{
	hours = time_remain / 3600;
	seconds = time_remain - hours*3600;
	if(seconds > 60){
	  minutes = seconds / 60;
	  seconds = seconds % 60;
	}
      }
      cout<<"Time remaining: "<<hours<<":"<<minutes<<":"<<seconds<<endl;
      cout<<endl;
      eventTree = 0;
      delete infile;
    }
  }

  TString Graphname ="pt/M cut on photon";

  //Non stacked plots
  TLegend *legend1 = new TLegend(0.6,0.75,0.85,0.85);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  for(int i=0; i<sig_N.size(); i++){
    cout<<sig_N[i]<<" "<<bkg_N[i]<<endl;
    if(bkg_N[i] == 0)
      SBratio.push_back(0);
    else
      SBratio.push_back(float(sig_N[i])/sqrt(float(bkg_N[i])));
  }
  int dim = SBratio.size();
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(2);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);


  TCanvas *c01 = new TCanvas("c01",Graphname,1200,900);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} (to a factor)");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0.5,10000000);
  //c01->SetLogy();
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  legend1->Clear();
  //legend1->Draw();
  c01->Print("p_ptm_SB.png");

  gBenchmark->Show("analyzerWG");

}
