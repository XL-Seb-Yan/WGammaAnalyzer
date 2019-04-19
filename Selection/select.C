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

void select(const TString conf="samples.conf", // input file
	      const TString outputDir=".",  // output directory
	      const Bool_t  doScaleCorr=0   // apply energy scale corrections?
	       ) {
  gBenchmark->Start("selectWG");
  gSystem->Load("/afs/cern.ch/work/x/xuyan/work5/CMSSW_9_4_0/src/VgammaTuplizer/Analyzer/Selection/lib/libMylib.so");
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
  TH1* hist01 = new TH1F("WGamma00",histotitle+", pt_{#gamma}",48,0,2400);
  TH1* hist02 = new TH1F("WGamma01",histotitle+", pt_{#gamma}",48,0,2400);
  
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;
  gStyle->SetOptStat(0);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  // parse .conf file
  confParse(conf, snamev, samplev);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);
  
  // Data structures to store info from produced flat ntuples
  Int_t runnum = -999;
  Int_t evtnum = -999;
  //--Photons---
  Int_t ph_N = -99;//for int we need to use this to initialize the container
  std::vector<float> *ph_pt = new std::vector<float>();
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
      //--Triggers
      eventTree->SetBranchAddress("HLT_isFired", &HLT_isFired);                           TBranch *HLTisFiredBr = eventTree->GetBranch("HLT_isFired");

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        //if(ientry%50000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;

	// Get Events
	runNumBr->GetEntry(ientry);
	evtNumBr->GetEntry(ientry);
	photonNBr->GetEntry(ientry);
	ph_pt->clear();                photonPtBr->GetEntry(ientry);
	HLTisFiredBr->GetEntry(ientry);


	std::vector<float> p_pt;

        for(vector<float>::iterator it = ph_pt->begin(); it != ph_pt->end(); it++)
	  p_pt.push_back(*it);

	//Only study events contain photons  -- trigger study
	if(ph_N < 1) 
	  continue;
	count1++;
	count3++;
	
	hist01->Fill(p_pt[0]);

	/*
	bool print = false;
	for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it){
	  if(it->second == 1){
	    cout<<it->first<<" "<<it->second<<endl;
	    print =true;	  
	  }
	}
	if (print)
	  cout<<endl;
	*/
	  
	bool passTrig = false;
	for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it) {
	  std::string trigName = it->first;
	  if (trigName.find("HLT_Photon175_v") != std::string::npos){
	    passTrig |= (1==it->second);
	  }
	}
	if (!passTrig) continue;
	count2++;
	count4++;

	hist02->Fill(p_pt[0]);
	/*
	if(p_pt[0]<150){
	  for(int i=0; i<p_pt.size();i++)
	    cout<<p_pt[i]<<" ";
	  cout<<endl;
	}
	*/
	
      }//end of event loop
      //cout<<"Number of events in this file: "<<eventTree->GetEntries()<<endl;
      //cout<<"Events with more than 1 photons: "<<count1<<endl;
      //cout<<"Events passed trigger: "<<count2<<endl;
      cout<<"Events passed trigger accmulated: "<<count4<<endl;
      
    }//end of file loop
  }//end of sample loop
  cout<<"Total number of events contain photon: "<<count3<<endl;
  cout<<"Total number of events passed trigger: "<<count4<<endl;

  //Tirgger Efficiency
  TEfficiency *Eff = new TEfficiency(*hist02, *hist01);
  
  //Non stacked plots
  hist01->SetLineWidth(2);

  TLegend *legend = new TLegend(0.7,0.85,0.9,0.9);

  TCanvas *c01 = new TCanvas("c01","pt_{#gamma}",1200,900);
  TAxis *xaxis = hist01->GetXaxis();
  TAxis *yaxis = hist01->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c01->SetLogy();
  c01->cd();
  hist01->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist01,"2017 SingleMuon B, with photons","f");
  legend->Draw();
  c01->Print("p_pt.png");

  TCanvas *c02 = new TCanvas("c02","Trigger Efficiency (HLT_Photon175)",1200,900);
  c02->cd();
  Eff->Draw();
  c02->Print("eff.png");

  TCanvas *c03 = new TCanvas("c03","pt_{#gamma}",1200,900);
  xaxis = hist02->GetXaxis();
  yaxis = hist02->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c03->SetLogy();
  c03->cd();
  hist02->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist02,"2017 SingleMuon B, fired trigger","f");
  legend->Draw();
  c03->Print("p_pt_fired.png");

  gBenchmark->Show("selectWG");

}
