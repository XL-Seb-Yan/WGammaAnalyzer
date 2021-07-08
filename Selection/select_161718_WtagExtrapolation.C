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
#include <utility>
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
#include <TStopwatch.h>
#include "TH1D.h"
#include "TRandom3.h"
#include "../Utils/interface/ConfParse.hh"             // input conf file parser
#include "../Utils/interface/CSample.hh"      // helper class to handle samples
#include "../Utils/interface/RunLumiRangeMap.hh"
//#include "../Utils/MyTools.hh"      // various helper functions
// C++ tool
#include <algorithm>
#include <map>
#endif

void select_161718_WtagExtrapolation(const TString conf="samples.conf", // input file
	          const TString outputDir=".",  // output directory
	          const Float_t weight=1   // MC weight?
	       ) {
  gBenchmark->Start("selectWG");
  gROOT->SetBatch(1);

  //Variable bin size
  // float xbinlow[5] = {100,300,600,1400,2400};
  float xbinlow[21] = {100,300,500,700,900,1100,1300,1500,1700,1900,2100,2300,2500,2700,2900,3100,3300,3500,3700,3900,4100};
  // float xbinlow[20] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190};
  // ------------------------Photons------------------------
  TH1F* hist00pa = new TH1F("WGamma00pa","Matched",20,&xbinlow[0]);
  TH1F* hist01pa = new TH1F("WGamma01pa","Pass",20,&xbinlow[0]);
  
  TH1F* hist = new TH1F("1","tau21",100,0,1);
 
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
  const TString ntupDir = outputDir;
  gSystem->mkdir(ntupDir,kTRUE); 
  
  // Data structures to store info from produced flat ntuples
  Int_t runnum = -999;
  Int_t evtnum = -999;
  Int_t lumiBlock = -999;
  //--Jets(AK8)--
  Int_t jetAK8_puppi_N = -99;
  std::vector<bool>  *jetAK8_puppi_IDTight = new std::vector<bool>();
  //std::vector<bool>  *jetAK8_puppi_IDTightLepVeto = new std::vector<bool>();
  std::vector<float> *jetAK8_puppi_softdrop_pt = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_eta = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_phi = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_mass = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_E = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_tau1 = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_tau2 = new std::vector<float>();
  //--GenParticle
  std::vector<float> *genPart_pt = new std::vector<float>();
  std::vector<float> *genPart_eta = new std::vector<float>();
  std::vector<float> *genPart_phi = new std::vector<float>();
  std::vector<float> *genPart_mass = new std::vector<float>();
  std::vector<int>   *genPart_pdgId = new std::vector<int>();
  std::vector<int>   *genPart_status = new std::vector<int>();
  std::vector<std::vector<int>>  *genPart_mother = new std::vector<std::vector<int>>;
  
  // loop over samples
  TTree* eventTree = 0;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    CSample* samp = samplev[isam];
    TStopwatch stopwatch;
    
    // loop through files
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      
      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile]<<endl; cout.flush();
      TFile *infile = TFile::Open(samp->fnamev[ifile]);
      assert(infile);

      // Access Event Tree
      TDirectory *folder;
      folder = (TDirectory*)infile->Get("ntuplizer");
      eventTree = (TTree*)folder->Get("tree");
      assert(eventTree);
      
      //--Jets (AK8 PUPPI)
      eventTree->SetBranchAddress("jetAK8_N", &jetAK8_puppi_N);                                  
      eventTree->SetBranchAddress("jetAK8_pt", &jetAK8_puppi_softdrop_pt);             
      eventTree->SetBranchAddress("jetAK8_eta", &jetAK8_puppi_softdrop_eta);            
      eventTree->SetBranchAddress("jetAK8_phi", &jetAK8_puppi_softdrop_phi);           
      eventTree->SetBranchAddress("jetAK8_e", &jetAK8_puppi_softdrop_E);               
      eventTree->SetBranchAddress("jetAK8_softdrop_massCorr", &jetAK8_puppi_softdrop_mass);      
      eventTree->SetBranchAddress("jetAK8_tau1", &jetAK8_puppi_tau1);                           
      eventTree->SetBranchAddress("jetAK8_tau2", &jetAK8_puppi_tau2);                           
      eventTree->SetBranchAddress("jetAK8_IDTight", &jetAK8_puppi_IDTight);
	  
	  //--GenParticles
	  eventTree->SetBranchAddress("genParticle_pt", &genPart_pt);
      eventTree->SetBranchAddress("genParticle_eta", &genPart_eta);
      eventTree->SetBranchAddress("genParticle_phi", &genPart_phi);
      eventTree->SetBranchAddress("genParticle_mass", &genPart_mass);
      eventTree->SetBranchAddress("genParticle_pdgId", &genPart_pdgId);
      eventTree->SetBranchAddress("genParticle_status", &genPart_status);
	  eventTree->SetBranchAddress("genParticle_mother", &genPart_mother);

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //
		UInt_t count0=0, count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;

		// Get Events
		jetAK8_puppi_softdrop_pt->clear();               
		jetAK8_puppi_softdrop_eta->clear();             
		jetAK8_puppi_softdrop_phi->clear();             
		jetAK8_puppi_softdrop_mass->clear();            
		jetAK8_puppi_softdrop_E->clear();               
		jetAK8_puppi_tau1->clear();             
		jetAK8_puppi_tau2->clear();            
		jetAK8_puppi_IDTight->clear();
		genPart_pt->clear();
	    genPart_eta->clear();
	    genPart_phi->clear();
	    genPart_mass->clear();
	    genPart_pdgId->clear();
	    genPart_status->clear();
		genPart_mother->clear();
		
		eventTree->GetEntry(ientry);

		// Only study events contain jets -- DATA -- 1st skimming
		if(jetAK8_puppi_N  == 0) continue;
		
		// Base
		// Status code
		int status_code = -99;
		if(conf.Contains("herwig"))
			status_code = 11;
		else
			status_code = 22;
		std::vector<TLorentzVector> *gen_W = new std::vector<TLorentzVector>();
		TLorentzVector gen_W_temp;
		for(int i=0; i<genPart_pt->size(); i++){
			//cout<<genPart_pdgId->at(i)<<","<<genPart_status->at(i)<<"    ";
			if(abs(genPart_pdgId->at(i)) == 24 && genPart_status->at(i) == status_code && abs(genPart_eta->at(i)) < 2){ //pythia
				//cout<<"++++++++OK++++++++++++"<<endl;
				// cout<<genPart_mother->at(i).size()<<",";
				for(int j=0; j<genPart_mother->at(i).size(); j++){
					//cout<<abs(genPart_mother->at(i).at(j))<<endl;
					 if(abs(genPart_mother->at(i).at(j)) == 9000007){
					//if(abs(genPart_mother->at(i).at(j)) == 39){
				       gen_W_temp.SetPtEtaPhiM(genPart_pt->at(i),genPart_eta->at(i),genPart_phi->at(i),genPart_mass->at(i));
					   gen_W->push_back(gen_W_temp);
					}
			    }
			}
		} 
		//cout<<endl;
		//cout<<"This event has: "<<gen_W->size()<<" GEN W(s)"<<endl;
		if(gen_W->size() == 0) continue;
		
		
		TLorentzVector reco_W_temp;
		for(int i=0; i<jetAK8_puppi_softdrop_pt->size(); i++){
			if(jetAK8_puppi_softdrop_pt->at(i) < 225) continue;
			if(abs(jetAK8_puppi_softdrop_eta->at(i)) > 2) continue;
			if(jetAK8_puppi_softdrop_mass->at(i) < 65 || jetAK8_puppi_softdrop_mass->at(i) > 105) continue;
			if(jetAK8_puppi_IDTight->at(i) != true) continue;
			reco_W_temp.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(i),jetAK8_puppi_softdrop_eta->at(i),jetAK8_puppi_softdrop_phi->at(i),jetAK8_puppi_softdrop_E->at(i));
			bool is_matched = false;
			for(int j=0; j<gen_W->size(); j++){
				if(reco_W_temp.DeltaR(gen_W->at(j)) < 0.3) is_matched = true;
			}
			float tau21 = jetAK8_puppi_tau2->at(i) / jetAK8_puppi_tau1->at(i);
			if(is_matched) {
				hist00pa->Fill(jetAK8_puppi_softdrop_pt->at(i));
				hist->Fill(tau21);
			}
			if(tau21 < 0.35 && is_matched) hist01pa->Fill(jetAK8_puppi_softdrop_pt->at(i));
		} 

      }//end of event loop
      //cout<<double(ifile)/double(nfiles)*100<<" % done with this dataset"<<endl;
      Double_t elapsed_t_file = stopwatch.RealTime() / (ifile+1);
      stopwatch.Start(kFALSE);
      Int_t time_remain = int(elapsed_t_file*(nfiles-ifile-1));
      Int_t hours = 0, minutes = 0, seconds = 0;
      if(time_remain < 60) seconds = time_remain;
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
      eventTree = 0;
      delete infile;
    }//end of file loop
  }//end of sample loop
  TString MCsample;
  if(conf.Contains("herwig"))
	  MCsample = "herwig";
  else
	  MCsample = "pythia";
  TFile *file = new TFile("WtaggingEffHighPt_"+MCsample+"_17.root","RECREATE");
  hist00pa->Write();
  hist01pa->Write();
  file->Close();
  gBenchmark->Show("selectWG");
  
  // TCanvas *c = new TCanvas("","",1200,900);
  // c->cd();
  // c->SetLogy();
  // hist->Scale(1/hist->GetEntries());
  // hist->GetYaxis()->SetRangeUser(0.0001,1);
  // hist->Draw();
  // c->Print("tau21.png");

}
