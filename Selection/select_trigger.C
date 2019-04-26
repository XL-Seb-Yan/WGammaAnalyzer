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
  TH1* hist01 = new TH1F("WGamma00",histotitle+", pt_{#gamma}",40,0,2000);
  TH1* hist02 = new TH1F("WGamma01",histotitle+", pt_{#gamma}",40,0,2000);
  TH1* hist03 = new TH1F("WGamma02",histotitle+", pt_{#gamma}",40,0,2000);
  
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
  TFile *infile=0;
  TTree *eventTree=0;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    CSample* samp = samplev[isam];
    cout<<"begin loop over files"<<endl;

    // loop through files
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      // Read input file and get the TTrees
      //cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "<<" this is "<<ifile<<" file out of "<<nfiles<<" ..."; cout.flush();
      infile = TFile::Open(samp->fnamev[ifile]);
      assert(infile);
      cout<<"infile::"<<infile<<endl;

      // counter reset
      count1 = 0;
      count2 = 0;

      // Access Event Tree
      TDirectory *folder;
      folder = (TDirectory*)infile->Get("ntuplizer");
      eventTree = (TTree*)folder->Get("tree");
      assert(eventTree);
      eventTree->SetBranchAddress("EVENT_run", &runnum);                        TBranch *runNumBr = eventTree->GetBranch("EVENT_run");
      eventTree->SetBranchAddress("EVENT_event", &evtnum);                      TBranch *evtNumBr = eventTree->GetBranch("EVENT_event");
      //--Photons--
      eventTree->SetBranchAddress("ph_N", &ph_N);                               TBranch *photonNBr = eventTree->GetBranch("ph_N");
      eventTree->SetBranchAddress("ph_pt", &ph_pt);                             TBranch *photonPtBr = eventTree->GetBranch("ph_pt");
       eventTree->SetBranchAddress("ph_passLooseId", &ph_passLooseId);           TBranch *photonPassLooseIdBr = eventTree->GetBranch("ph_passLooseId");
      //--Triggers
      eventTree->SetBranchAddress("HLT_isFired", &HLT_isFired);                 TBranch *HLTisFiredBr = eventTree->GetBranch("HLT_isFired");

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        //if(ientry%50000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;

	// Get Events
	runNumBr->GetEntry(ientry);
	evtNumBr->GetEntry(ientry);
	photonNBr->GetEntry(ientry);
	ph_pt->clear();                photonPtBr->GetEntry(ientry);
	ph_passLooseId->clear();       photonPassLooseIdBr->GetEntry(ientry);
	HLT_isFired->clear();          HLTisFiredBr->GetEntry(ientry);


	std::vector<float> p_pt;
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
	
	if(p_looseID[0] != 1) //validating EM obejects to be photons
	  continue;
	count1++;
	count3++;
	
	if(p_pt[0]>15)
	  count6++;
	
	hist02->Fill(p_pt[0]);

	bool passTrig = false;
	for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it) {
	  if (it->first.find("HLT_Photon175_") != std::string::npos && it->second == 1)
	    passTrig = true;
	  if (it->first.find("HLT_Photon165_HE10_") != std::string::npos && it->second == 1)
	    passTrig = true;
	  if (it->first.find("HLT_Photon200") != std::string::npos && it->second == 1)
	    passTrig = true;
	}
	if (!passTrig) continue;
	count2++;
	count4++;
	hist03->Fill(p_pt[0]);
	
      }//end of event loop
      //cout<<"Number of events in this file: "<<eventTree->GetEntries()<<endl;
      //cout<<"Events with more than 1 photons: "<<count1<<endl;
      //cout<<"Events passed trigger: "<<count2<<endl;
      cout<<"Events passed trigger accmulated: "<<count4<<endl;
      
    }//end of file loop
  }//end of sample loop
  cout<<"Total number of events contain photon: "<<count3<<endl;
  cout<<"Total number of events passed trigger: "<<count4<<endl;
  cout<<"Total number of events: "<<count5<<endl;
  cout<<"Total number of events with pt>15: "<<count6<<endl;
  //Tirgger Efficiency
  TEfficiency *Eff = new TEfficiency(*hist03, *hist02);

  //Non stacked plots

  TLegend *legend = new TLegend(0.65,0.8,0.9,0.9);

  TCanvas *c01 = new TCanvas("c01","pt_{#gamma}",1200,900);
  TAxis *xaxis = hist01->GetXaxis();
  TAxis *yaxis = hist01->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c01->SetLogy();
  c01->cd();
  hist01->SetLineWidth(2);
  hist01->Draw("HIST");
  hist02->SetLineWidth(2);
  hist02->SetLineColor(2);
  hist02->Draw("SAMEHIST");
  hist03->SetFillColor(30);
  hist03->SetLineColor(30);
  hist03->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hist01,"2017 SingleMuon B, with EM objects","f");
  legend->AddEntry(hist02,"2017 SingleMuon B, pass photon loose ID","f");
  legend->AddEntry(hist03,"2017 SingleMuon B, fired trigger","f");
  legend->Draw();
  c01->Print("p_pt.png");

  TCanvas *c02 = new TCanvas("c02","Trigger Efficiency",1200,900);
  c02->cd();
  Eff->Draw();
  c02->Print("eff.png");

  gBenchmark->Show("selectWG");

}
