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

void select201817_GENAcc(const TString conf="samples.conf", // input file
	          const TString outputDir=".",  // output directory
	          const Float_t weight=1   // MC weight?
	       ) {
  gBenchmark->Start("selectWG");
  gROOT->SetBatch(1);
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //=============================================================================================================
  bool isPlot = false; 
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
  //--GenParticle
  std::vector<float> *genPart_pt = new std::vector<float>();
  std::vector<float> *genPart_eta = new std::vector<float>();
  std::vector<float> *genPart_phi = new std::vector<float>();
  std::vector<float> *genPart_mass = new std::vector<float>();
  std::vector<int>   *genPart_pdgId = new std::vector<int>();
  std::vector<int>   *genPart_status = new std::vector<int>();
  std::vector<vector<int>>   *genPart_mother = new std::vector<vector<int>>();
  std::vector<vector<float>> *genPart_mother_pt = new std::vector<vector<float>>();
  std::vector<vector<float>> *genPart_mother_eta = new std::vector<vector<float>>();
  std::vector<vector<float>> *genPart_mother_phi = new std::vector<vector<float>>();
  std::vector<vector<float>> *genPart_mother_e = new std::vector<vector<float>>();
  std::vector<vector<int>>   *genPart_daughter = new std::vector<vector<int>>();
  std::vector<vector<float>> *genPart_daughter_pt = new std::vector<vector<float>>();
  std::vector<vector<float>> *genPart_daughter_eta = new std::vector<vector<float>>();
  std::vector<vector<float>> *genPart_daughter_phi = new std::vector<vector<float>>();
  std::vector<vector<float>> *genPart_daughter_e = new std::vector<vector<float>>();

  
  // loop over samples
  TTree* eventTree = 0;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
	UInt_t count0=0, count1=0;
	CSample* samp = samplev[isam];
	string samplename = snamev[isam].Data();
	std::stringstream ss(samplename);
	string substr;
    getline(ss,substr,'_');
	float GENmass = strtof((substr).c_str(),0);
	cout<<GENmass<<endl;

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
		TDirectory *folder;
		folder = (TDirectory*)infile->Get("ntuplizer");
		eventTree = (TTree*)folder->Get("tree");
		assert(eventTree);
		//--GenPartticle
		eventTree->SetBranchAddress("genParticle_pt", &genPart_pt);
		eventTree->SetBranchAddress("genParticle_eta", &genPart_eta);
		eventTree->SetBranchAddress("genParticle_phi", &genPart_phi);
		eventTree->SetBranchAddress("genParticle_mass", &genPart_mass);
		eventTree->SetBranchAddress("genParticle_pdgId", &genPart_pdgId);
		eventTree->SetBranchAddress("genParticle_status", &genPart_status);
		eventTree->SetBranchAddress("genParticle_mother", &genPart_mother);
		eventTree->SetBranchAddress("genParticle_mother_e", &genPart_mother_e);
		eventTree->SetBranchAddress("genParticle_mother_pt", &genPart_mother_pt);
		eventTree->SetBranchAddress("genParticle_mother_eta", &genPart_mother_eta);
		eventTree->SetBranchAddress("genParticle_mother_phi", &genPart_mother_phi);
		eventTree->SetBranchAddress("genParticle_dau", &genPart_daughter);
		eventTree->SetBranchAddress("genParticle_dau_pt", &genPart_daughter_pt);
		eventTree->SetBranchAddress("genParticle_dau_eta", &genPart_daughter_eta);
		eventTree->SetBranchAddress("genParticle_dau_phi", &genPart_daughter_phi);
		eventTree->SetBranchAddress("genParticle_dau_e", &genPart_daughter_e);

		for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
		// for(UInt_t ientry=0; ientry<10; ientry++) {
			genPart_pt->clear();
			genPart_eta->clear();
			genPart_phi->clear();
			genPart_mass->clear();
			genPart_pdgId->clear();
			genPart_status->clear();
			genPart_mother->clear();
			genPart_mother_pt->clear();
			genPart_mother_eta->clear();
			genPart_mother_phi->clear();
			genPart_mother_e->clear();
			genPart_daughter->clear();
			genPart_daughter_pt->clear();
			genPart_daughter_eta->clear();
			genPart_daughter_phi->clear();
			genPart_daughter_e->clear();
			eventTree->GetEntry(ientry);
			count0++;
			//-------locating candidate photon and W boson
				
			TLorentzVector gen_p;
			int index_p_gen = -99;
			for(int i=0; i<genPart_pt->size(); i++){
				if(genPart_pdgId->at(i) == 22 && genPart_status->at(i) == 1 && abs(genPart_eta->at(i)) < 1.44 && genPart_pt->at(i) > 225){
					for(int j=0; j<genPart_mother->at(i).size(); j++){
						if(abs(genPart_mother->at(i).at(j)) == 9000007){
						   gen_p.SetPtEtaPhiM(genPart_pt->at(i),genPart_eta->at(i),genPart_phi->at(i),genPart_mass->at(i));
						   index_p_gen = i;
						   break;
						}
					}
				}
			} 
			
			TLorentzVector gen_W;
			int index_W_gen = -99;
			for(int i=0; i<genPart_pt->size(); i++){
				if(abs(genPart_pdgId->at(i)) == 24 && abs(genPart_eta->at(i)) < 2 && genPart_pt->at(i) > 225){
					for(int j=0; j<genPart_mother->at(i).size(); j++){
						if(abs(genPart_mother->at(i).at(j)) == 9000007){
						   gen_W.SetPtEtaPhiM(genPart_pt->at(i),genPart_eta->at(i),genPart_phi->at(i),genPart_mass->at(i));
						   index_W_gen = i;
						   break;
						}
					}
				}
			} 
			
			if(index_p_gen == -99 || index_W_gen == -99){
				continue;
			}
			
			//-----------------------------------------W+Gamma system--------------------------------------------------
			TLorentzVector v_sys;
			v_sys = gen_W + gen_p;

			// Calculate invariant mass and pt/m
			Double_t invmass = 0;
			Double_t ptoverM = 0;
			invmass = (v_sys).M();
			ptoverM = genPart_pt->at(index_p_gen)/invmass;

			// Calculate cos(theta*)
			Double_t cosThetaStar = -999;
			TLorentzVector v_boosted_j, v_boosted_p;
			v_boosted_p = gen_p;
			v_boosted_p.Boost(-(v_sys.BoostVector()));
			cosThetaStar = abs(v_boosted_p.Pz()/v_boosted_p.P());

			// Calculate seperation
			Double_t seperation = gen_W.DeltaR(gen_p);
			
			//cuts:
			// cout<<invmass<<" "<<ptoverM<<" "<<cosThetaStar<<" "<<seperation<<endl;
			if(invmass < 600) continue;
			if(ptoverM < 0.37) continue;
			if(cosThetaStar > 0.6) continue;
			if(seperation < 1.1) continue;
			if(invmass < GENmass*0.75 || invmass > GENmass*1.25) continue;
			count1++;
	
      }//end of event loop
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
      eventTree = 0;
      delete infile;
    }//end of file loop
	cout<<"Number of events: "<<count0<<endl;
	cout<<"Number of events within acceptance: "<<count1<<endl;
	cout<<"Acceptance: "<<snamev[isam]<<" "<<float(count1)/float(count0)<<endl;
  }//end of sample loop
  gBenchmark->Show("selectWG");

}