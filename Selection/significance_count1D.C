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

void significance_count1D(const TString conf="samples.conf"){

  gBenchmark->Start("analyzerWG");
  gSystem->Load("lib/libMylib.so");
  gROOT->SetBatch(1);
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


  TH1F* hist01 = new TH1F("WGamma01","BDT response, signal MC (M-800)",50,-0.5,0.5);
  TH1F* hist02 = new TH1F("WGamma02","BDT response, data sideband",50,-0.5,0.5);
  
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;

  Int_t Nbkg = 0, Nsig = 0;
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
  float sys_BDT;


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
      eventTree->SetBranchAddress("photon_pt", &photon_pt);                       
      eventTree->SetBranchAddress("photon_eta", &photon_eta);                       
      eventTree->SetBranchAddress("photon_phi", &photon_phi);                         
      eventTree->SetBranchAddress("photon_e", &photon_e);                                             
      eventTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
      eventTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
      eventTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
      eventTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
      eventTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
      eventTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
      eventTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
      eventTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
      eventTree->SetBranchAddress("sys_invmass", &sys_invmass); 

      if(isam == 0){
	Nsig = eventTree->GetEntries();
	for(float cut2 = 0; cut2 < 1; cut2+=0.01){
	  int N = 0;
	  float record = cut2;
	  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	  //for(UInt_t ientry=0; ientry<1400; ientry++) {
	    // Get Events
	   eventTree->GetEntry(ientry);

	    if(abs(photon_eta) > 1.44) continue;
	    if(abs(ak8puppijet_eta) > 2.0) continue;
	    //if(sys_costhetastar > 0.65) continue;
	    //if(ak8puppijet_tau21 > 0.4) continue;
	    //if(sys_ptoverm < 0.37) continue;
	    //if(abs(ak8puppijet_masssoftdropcorr - 80.38) > 15) continue;

	    //if(abs(photon_eta) < cut2){
	    //if(abs(ak8puppijet_eta) < cut2){
	    //if(sys_costhetastar < cut2){
	    //if(ak8puppijet_tau21 < cut2){
	    if(sys_ptoverm > cut2){
	    //if(abs(ak8puppijet_masssoftdropcorr - 80.38) < cut2){
	      N++;
	    }
	  }
	  sig_N.push_back(N);
	  recordarr.push_back(record);
	}
      }

      if(isam == 1){
	Nbkg = eventTree->GetEntries();
	for(float cut2 = 0; cut2 < 1; cut2+=0.01){
	  int N = 0;
	  float record = cut2;
	  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	  //for(UInt_t ientry=0; ientry<1000; ientry++) {
	    // Get Events
	    eventTree->GetEntry(ientry);
	    
	    if(abs(photon_eta) > 1.44) continue;
	    if(abs(ak8puppijet_eta) > 2.0) continue;
	    //if(sys_costhetastar > 0.65) continue;
	    //if(ak8puppijet_tau21 > 0.4) continue;
	    //if(sys_ptoverm < 0.37) continue;
	    //if(abs(ak8puppijet_masssoftdropcorr - 80.38) > 15) continue;

	    //if(abs(photon_eta) < cut2){
	    //if(abs(ak8puppijet_eta) < cut2){
	    //if(sys_costhetastar < cut2){
	    //if(ak8puppijet_tau21 < cut2){
	    if(sys_ptoverm > cut2){
	    //if(abs(ak8puppijet_masssoftdropcorr - 80.38) < cut2){
	      N++;
	    }
	  }
	  bkg_N.push_back(N);
	  recordarr.push_back(record);
	}
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

  //TString Graphname ="|#eta| cut on jet";
  //TString Graphname ="|#eta| cut on photon";
  //TString Graphname ="cos(#theta*) cut on photon";
  //TString Graphname ="tau21 cut on jet";
  TString Graphname ="pt/M cut on photon";
  //TString Graphname ="mass window cut on jet";

  //Non stacked plots
  TLegend *legend1 = new TLegend(0.15,0.75,0.3,0.85);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
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

  TGraph* gr1 = new TGraph(dim,&recordarr[0],&Sigeff[0]);
  gr1->SetTitle(Graphname+" efficiency");
  gr1->SetLineColor(6);
  gr1->SetLineWidth(4);
  //gr1->SetMarkerColor(4);
  //gr1->SetMarkerSize(1.5);

  TGraph* gr2 = new TGraph(dim,&recordarr[0],&Bkgeff[0]);
  //gr1->SetTitle(Graphname+" efficiency");
  gr2->SetLineColor(4);
  gr2->SetLineWidth(4);
  //gr2->SetMarkerColor(4);
  //gr2->SetMarkerSize(1.5);


  TCanvas *c01 = new TCanvas("c01",Graphname,1200,900);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //xaxis->SetRangeUser(-0.5,0.5);
  xaxis->SetRangeUser(0,1);
  yaxis->SetRangeUser(0,200);
  //yaxis->SetRangeUser(0.5,10000000);
  //c01->SetLogy();
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  //legend1->Clear();
  //legend1->Draw();  
  c01->Print("700_ptm.png");

  TCanvas *c02 = new TCanvas("c02",Graphname+" Efficiency",1200,900);
  xaxis = gr1->GetXaxis();
  yaxis = gr1->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("Eff");
  yaxis->SetTitleOffset(1.1);
  xaxis->SetTitleOffset(1.3);
  //xaxis->SetRangeUser(-0.5,0.5);
  xaxis->SetRangeUser(0,1);
  yaxis->SetRangeUser(0,1);
  //yaxis->SetRangeUser(0.5,10000000);
  //c01->SetLogy();
  c02->cd();
  c02->SetGrid();
  cout<<"OK"<<endl;
  gr1->Draw("ACY+");
  gr2->Draw("SAME");
  legend1->Clear();
  legend1->AddEntry(gr1,"signal efficiency");
  legend1->AddEntry(gr2,"background efficiency");
  //legend1->Draw();
  c02->Print("700_ptm.png");

  gBenchmark->Show("analyzerWG");

}
