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
#include "TH1D.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
//#include "../Utils/MyTools.hh"      // various helper functions
// C++ tool
#include <algorithm>
#include <map>
#endif

void select_trigger(const TString conf="samples.conf", // input file
	      const TString outputDir=".",  // output directory
	      const Bool_t  doScaleCorr=0   // apply energy scale corrections?
	       ) {
  gBenchmark->Start("selectWG");
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

  //Photons
  TH1* hist01 = new TH1F("WGamma00",histotitle+", pt_{#gamma}",50,0,2500);
  TH1* hist02 = new TH1F("WGamma01",histotitle+", pt_{#gamma}",50,0,2500);
  TH1* hist03 = new TH1F("Trigger150",histotitle+", pt_{#gamma}",50,0,2500);
  TH1* hist04 = new TH1F("Trigger165",histotitle+", pt_{#gamma}",50,0,2500);
  TH1* hist05 = new TH1F("Trigger175",histotitle+", pt_{#gamma}",50,0,2500);
  TH1* hist06 = new TH1F("Trigger200",histotitle+", pt_{#gamma}",50,0,2500);
  
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
  Int_t runnum = -999;
  Int_t evtnum = -999;
  //--Photons---
  Int_t ph_N = -99;//for int we need to use this to initialize the container
  std::vector<float> *ph_pt = new std::vector<float>();
  std::vector<int>   *ph_passLooseId = new std::vector<int>();
  //Trigger
  std::map<std::string,bool> *HLT_isFired = new std::map<std::string,bool>();
  
  // loop over samples
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    CSample* samp = samplev[isam];

    // loop through files
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      // Read input file and get the TTrees
      if(ifile%10 == 0)
	cout << "Processing "<<ifile<<" file out of "<<nfiles<<" and this is "<<isam<<" samples"; cout.flush();
      TFile *infile = TFile::Open(samp->fnamev[ifile]);
      assert(infile);
      if(ifile%10 == 0)
	cout<<"infile::"<<infile<<endl;

      // counter reset
      count1 = 0;
      count2 = 0;

      // Access Event Tree
      TDirectory *folder;
      folder = (TDirectory*)infile->Get("ntuplizer");
      TTree *eventTree = (TTree*)folder->Get("tree");
      assert(eventTree);
      //--Photons--
      eventTree->SetBranchAddress("ph_N", &ph_N);                               TBranch *photonNBr = eventTree->GetBranch("ph_N");
      eventTree->SetBranchAddress("ph_pt", &ph_pt);                             TBranch *photonPtBr = eventTree->GetBranch("ph_pt");
       eventTree->SetBranchAddress("ph_passLooseId", &ph_passLooseId);           TBranch *photonPassLooseIdBr = eventTree->GetBranch("ph_passLooseId");
      //--Triggers
      eventTree->SetBranchAddress("HLT_isFired", &HLT_isFired);                 TBranch *HLTisFiredBr = eventTree->GetBranch("HLT_isFired");

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	
	// Get Events
	photonNBr->GetEntry(ientry);
	ph_pt->clear();                photonPtBr->GetEntry(ientry);
	ph_passLooseId->clear();       photonPassLooseIdBr->GetEntry(ientry);
	HLT_isFired->clear();          HLTisFiredBr->GetEntry(ientry);


	std::vector<float> p_pt;;
	std::vector<int> p_looseID;

        for(vector<float>::iterator it = ph_pt->begin(); it != ph_pt->end(); it++)
	  p_pt.push_back(*it);
	for(vector<int>::iterator it = ph_passLooseId->begin(); it != ph_passLooseId->end(); it++)
	  p_looseID.push_back(*it);

	count5++;

	//Only study events contain photons  -- trigger study
	if(ph_N < 1)
	  continue;
	hist01->Fill(p_pt[0]);
	
	if(p_looseID[0] < 1) //validating EM obejects to be photons
	  continue;
	count1++;
	count3++;
		
	hist02->Fill(p_pt[0]);

	bool passTrig = false;
	for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it) {
	  if (it->first.find("HLT_Photon150") != std::string::npos && it->second == 1){
	    hist03->Fill(p_pt[0]);
	    passTrig = true;
	  }
	  if (it->first.find("HLT_Photon165") != std::string::npos && it->second == 1){
	    hist04->Fill(p_pt[0]);
	    passTrig = true;
	  }
	  if (it->first.find("HLT_Photon175") != std::string::npos && it->second == 1){
	    hist05->Fill(p_pt[0]);
	    passTrig = true;
	  }
	  if (it->first.find("HLT_Photon200") != std::string::npos && it->second == 1){
	    hist06->Fill(p_pt[0]);
	    passTrig = true;
	  }
	}
	if (!passTrig) continue;
	count2++;
	count4++;
	
      }//end of event loop
      //cout<<"Number of events in this file: "<<eventTree->GetEntries()<<endl;
      //cout<<"Events with more than 1 photons: "<<count1<<endl;
      //cout<<"Events passed trigger: "<<count2<<endl;
      if(ifile%100 == 0)
	cout<<"Events passed trigger accmulated: "<<count4<<endl;
      infile->Close();
    }//end of file loop
    
  }//end of sample loop
  cout<<"Total number of events contain photon: "<<count3<<endl;
  cout<<"Total number of events passed trigger: "<<count4<<endl;
  cout<<"Total number of events: "<<count5<<endl;
  //Tirgger Efficiency
  TEfficiency *Eff1 = new TEfficiency(*hist03, *hist02);
  TEfficiency *Eff2 = new TEfficiency(*hist04, *hist02);
  TEfficiency *Eff3 = new TEfficiency(*hist05, *hist02);
  TEfficiency *Eff4 = new TEfficiency(*hist06, *hist02);

  //Non stacked plots

  TLegend *legend1 = new TLegend(0.55,0.7,0.85,0.85);
  TLegend *legend2 = new TLegend(0.55,0.7,0.85,0.85);

  TCanvas *c01 = new TCanvas("c01","pt_{#gamma}",1200,900);
  TAxis *xaxis = hist01->GetXaxis();
  TAxis *yaxis = hist01->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c01->SetLogy();
  c01->cd();
  hist01->SetLineWidth(2);
  hist01->SetLineColor(1);
  hist01->Draw("HIST");
  hist02->SetLineWidth(2);
  hist02->SetLineColor(4);
  hist02->Draw("SAMEHIST");
  hist03->SetLineWidth(2);
  hist03->SetLineColor(9);
  hist03->Draw("SAMEHIST");
  hist04->SetLineWidth(2);
  hist04->SetLineColor(6);
  hist04->Draw("SAMEHIST");
  hist05->SetLineWidth(2);
  hist05->SetLineColor(7);
  hist05->Draw("SAMEHIST");
  hist06->SetLineWidth(2);
  hist06->SetLineColor(8);
  hist06->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist01,"2017 SingleMuon, with EM objects","f");
  legend1->AddEntry(hist02,"2017 SingleMuon, pass photon loose ID","f");
  legend1->AddEntry(hist03,"2017 SingleMuon, fired HLT_Photon150_","f");
  legend1->AddEntry(hist04,"2017 SingleMuon, fired HLT_Photon165_","f");
  legend1->AddEntry(hist05,"2017 SingleMuon, fired HLT_Photon175_","f");
  legend1->AddEntry(hist06,"2017 SingleMuon, fired HLT_Photon200_","f");
  legend1->Draw();
  c01->Print("p_pt.png");

  TCanvas *c02 = new TCanvas("c02","Trigger Efficiency",1200,900);
  c02->cd();
  Eff1->SetTitle("Trigger Efficiency");
  Eff1->SetLineWidth(2);
  Eff1->SetLineColor(9);
  Eff1->Draw();
  gPad->Update();
  Eff1->GetPaintedGraph()->GetXaxis()->SetTitle("pt_{#gamma} (GeV)");
  Eff1->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");
  Eff1->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1);
  Eff2->SetLineWidth(2);
  Eff2->SetLineColor(6);
  Eff2->Draw("SAME");
  Eff3->SetLineWidth(2);
  Eff3->SetLineColor(7);
  Eff3->Draw("SAME");
  Eff4->SetLineWidth(2);
  Eff4->SetLineColor(8);
  Eff4->Draw("SAME");
  legend2->Clear();
  legend2->AddEntry(Eff1,"HLT_Photon150_","f");
  legend2->AddEntry(Eff2,"HLT_Photon165_","f");
  legend2->AddEntry(Eff3,"HLT_Photon175_","f");
  legend2->AddEntry(Eff4,"HLT_Photon200_","f");
  legend2->Draw();
  c02->Print("eff.png");

  gBenchmark->Show("selectWG");

}
