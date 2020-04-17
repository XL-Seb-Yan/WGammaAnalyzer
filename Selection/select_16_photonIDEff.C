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

void select_16_photonIDEff(const TString conf="samples.conf", // input file
	          const TString outputDir=".",  // output directory
	          const Float_t weight=1   // MC weight?
	       ) {
  gBenchmark->Start("selectWG");
  gROOT->SetBatch(1);

  //Variable bin size
  float xbinlow[20] = {0,25,50,75,100,150,200,250,300,350,400,500,600,700,800,1000,1200,1600,2000,3000};
  // ------------------------Photons------------------------
  TH1F* hist00pa = new TH1F("WGamma00pa","Base",19,&xbinlow[0]);
  TH1F* hist01pa = new TH1F("WGamma01pa","LID",19,&xbinlow[0]);
  TH1F* hist02pa = new TH1F("WGamma02pa","MID",19,&xbinlow[0]);
  TH1F* hist03pa = new TH1F("WGamma03pa","TID",19,&xbinlow[0]);
  TH1F* hist04pa = new TH1F("WGamma04pa","80ID",19,&xbinlow[0]);
  TH1F* hist05pa = new TH1F("WGamma05pa","90ID",19,&xbinlow[0]);
 
  
  UInt_t count0=0, count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;
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

  // Data Structure for output skimmed files
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat, photon_corr, photon_energyscale, photon_energyscale_up, photon_energyscale_down, photon_resolution;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21, ak8puppijet_massdiff;
  float ak8puppijet_jec, ak8puppijet_jec_up, ak8puppijet_jec_down, ak8puppijet_jer, ak8puppijet_jer_sf, ak8puppijet_jer_sf_up, ak8puppijet_jer_sf_down;
  // System
  float sys_costhetastar, sys_ptoverm, sys_invmass, sys_seperation;
  // MC xsec weight
  float xsec_weight;
  
  // Data structures to store info from produced flat ntuples
  Int_t runnum = -999;
  Int_t evtnum = -999;
  Int_t lumiBlock = -999;
  
  //--Photons---
  Int_t ph_N = -99;//for int we need to use this to initialize the container
  std::vector<float> *ph_pt = new std::vector<float>();
  std::vector<float> *ph_eta = new std::vector<float>();
  std::vector<float> *ph_phi = new std::vector<float>();
  std::vector<float> *ph_E = new std::vector<float>();
  std::vector<float> *ph_Et = new std::vector<float>();
  std::vector<float> *ph_m = new std::vector<float>();
  std::vector<float> *ph_superCluster_eta = new std::vector<float>();
  std::vector<float> *ph_superCluster_phi = new std::vector<float>();
  std::vector<float> *ph_sigmaIetaIeta = new std::vector<float>();
  std::vector<float> *ph_hOverE = new std::vector<float>();
  std::vector<float> *ph_isoGamma = new std::vector<float>();
  std::vector<float> *ph_isoCh = new std::vector<float>();
  std::vector<bool>  *ph_passEleVeto = new std::vector<bool>();
  std::vector<int>   *ph_passLooseId = new std::vector<int>(); //check this 2.15
  std::vector<int>   *ph_passMediumId = new std::vector<int>();
  std::vector<int>   *ph_passTightId = new std::vector<int>();
  std::vector<float> *ph_mvaVal = new std::vector<float>();
  std::vector<float> *ph_mvaCat = new std::vector<float>();
  std::vector<float> *ph_corr = new std::vector<float>();
  std::vector<float> *ph_energyscale = new std::vector<float>();
  std::vector<float> *ph_energyscale_up = new std::vector<float>();
  std::vector<float> *ph_energyscale_down = new std::vector<float>();
  std::vector<float> *ph_resolution = new std::vector<float>();
  //--Jets(AK8)--
  Int_t jetAK8_puppi_N = -99;
  //--Trigger
  std::map<std::string,bool> *HLT_isFired = new std::map<std::string,bool>();
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
      
      //--Photons--
      eventTree->SetBranchAddress("ph_N", &ph_N);                              
      eventTree->SetBranchAddress("ph_pt", &ph_pt);                            
      eventTree->SetBranchAddress("ph_eta", &ph_eta);                           
      eventTree->SetBranchAddress("ph_phi", &ph_phi);                          
      eventTree->SetBranchAddress("ph_e", &ph_E);                                             
      eventTree->SetBranchAddress("ph_passEleVeto", &ph_passEleVeto);          
      eventTree->SetBranchAddress("ph_passLooseId", &ph_passLooseId);           
      eventTree->SetBranchAddress("ph_passMediumId", &ph_passMediumId);         
      eventTree->SetBranchAddress("ph_passTightId", &ph_passTightId);           
      eventTree->SetBranchAddress("ph_mvaVal", &ph_mvaVal);                    
      eventTree->SetBranchAddress("ph_mvaCat", &ph_mvaCat);  
      eventTree->SetBranchAddress("ph_Corr", &ph_corr);                    
      eventTree->SetBranchAddress("ph_energyscale", &ph_energyscale); 	 
	  eventTree->SetBranchAddress("ph_energyscale_up", &ph_energyscale_up); 	
	  eventTree->SetBranchAddress("ph_energyscale_down", &ph_energyscale_down); 		  
	  eventTree->SetBranchAddress("ph_resolution", &ph_resolution); 
      //--Jets (AK8 PUPPI)
      eventTree->SetBranchAddress("jetAK8_N", &jetAK8_puppi_N);                                  
      eventTree->SetBranchAddress("HLT_isFired", &HLT_isFired);
	  //--GenParticles
	  eventTree->SetBranchAddress("genParticle_pt", &genPart_pt);
      eventTree->SetBranchAddress("genParticle_eta", &genPart_eta);
      eventTree->SetBranchAddress("genParticle_phi", &genPart_phi);
      eventTree->SetBranchAddress("genParticle_mass", &genPart_mass);
      eventTree->SetBranchAddress("genParticle_pdgId", &genPart_pdgId);
      eventTree->SetBranchAddress("genParticle_status", &genPart_status);
	  eventTree->SetBranchAddress("genParticle_mother", &genPart_mother);

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<2; ientry++) {
		count1++;

		// Get Events
		ph_pt->clear();               
		ph_eta->clear();            
		ph_phi->clear();              
		ph_E->clear();                          
		ph_passEleVeto->clear();       
		ph_passLooseId->clear();       
		ph_passMediumId->clear();     
		ph_passTightId->clear();      
		ph_mvaVal->clear();           
		ph_mvaCat->clear(); 
		ph_corr->clear();     
		ph_energyscale->clear();     
		ph_energyscale_up->clear();     
		ph_energyscale_down->clear();     
		ph_resolution->clear();     
		HLT_isFired->clear();
		genPart_pt->clear();
	    genPart_eta->clear();
	    genPart_phi->clear();
	    genPart_mass->clear();
	    genPart_pdgId->clear();
	    genPart_status->clear();
		genPart_mother->clear();
		
		// -------------------------------------Trigger----------------------------------------------
		// HLT Trigger decision first, improve speed
		TBranch *HLTBranch = (TBranch*)eventTree->GetBranch("HLT_isFired");
		HLTBranch->GetEntry(ientry);
		
		bool passTrig = false;
		for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it) {
		  if (it->first.find("HLT_Photon165_HE10_v") != std::string::npos && it->second == 1)
			passTrig = true;
	      if (it->first.find("HLT_Photon165_R9Id90_HE10_IsoM_v") != std::string::npos && it->second == 1)
			passTrig = true;
	      if (it->first.find("HLT_Photon175_v") != std::string::npos && it->second == 1)
			passTrig = true;
		}
		// if (!passTrig) continue;
		count0++;
		eventTree->GetEntry(ientry);

		// Only study events contain photons (EM objects here) and jets -- DATA -- 1st skimming
		if(ph_N  == 0 || jetAK8_puppi_N  == 0) continue;
		count2++;
		
		// -----------------------------------Photon-----------------------------------------------------
		// Photon selection Flags
		Bool_t pass_p1 = false;
		Bool_t pass_p2 = false;
		
		// Base
		TLorentzVector gen_ph, reco_ph;
		int index_p_gen = -99;
		for(int i=0; i<genPart_pt->size(); i++){
			if(genPart_pdgId->at(i) == 22 && genPart_status->at(i) == 1 && abs(genPart_eta->at(i)) < 1.44){
				for(int j=0; j<genPart_mother->at(i).size(); j++){
					if(abs(genPart_mother->at(i).at(j)) == 9000007){
				       gen_ph.SetPtEtaPhiM(genPart_pt->at(i),genPart_eta->at(i),genPart_phi->at(i),genPart_mass->at(i));
				       index_p_gen = i;
				       break;
					}
			    }
			}
		} 
		if(index_p_gen == -99) continue;
		count3++;
		int index_p = -99;
		for(int i=0; i<ph_pt->size(); i++){
			reco_ph.SetPtEtaPhiE(ph_pt->at(i),ph_eta->at(i),ph_phi->at(i),ph_E->at(i));
			if(reco_ph.DeltaR(gen_ph) < 0.3) index_p = i;
		} 
		if(index_p == -99) continue;
		count4++;
		hist00pa->Fill(ph_pt->at(index_p)*ph_corr->at(index_p));
		
		if(ph_passLooseId->at(index_p)) hist01pa->Fill(ph_pt->at(index_p)*ph_corr->at(index_p));
		if(ph_passMediumId->at(index_p)) hist02pa->Fill(ph_pt->at(index_p)*ph_corr->at(index_p));
		if(ph_passTightId->at(index_p)) hist03pa->Fill(ph_pt->at(index_p)*ph_corr->at(index_p));
		if(ph_mvaVal->at(index_p) > 0.42) hist04pa->Fill(ph_pt->at(index_p)*ph_corr->at(index_p));
		if(ph_mvaVal->at(index_p) > -0.02) hist05pa->Fill(ph_pt->at(index_p)*ph_corr->at(index_p));

      }//end of event loop
      cout<<double(ifile)/double(nfiles)*100<<" % done with this dataset"<<endl;
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
      cout<<"Number of events: "<<count0<<endl;
      cout<<"Number of events have photon and AK8 jet: "<<count1<<endl;
      cout<<"Number of events fired trigger: "<<count2<<endl;
      cout<<"Number of events pass photon pre-selection: "<<count3<<endl;
      cout<<"Number of events pass jet pre-selection: "<<count4<<endl;
      cout<<"Number of events pass photon mvaID: "<<count5<<endl;
      cout<<"Time remaining: "<<hours<<":"<<minutes<<":"<<seconds<<endl;
      eventTree = 0;
      delete infile;
    }//end of file loop
  }//end of sample loop
  TFile *file = new TFile("PhotonIDEff.root","RECREATE");
  hist00pa->Write();
  hist01pa->Write();
  hist02pa->Write();
  hist03pa->Write();
  hist04pa->Write();
  hist05pa->Write();
  file->Close();
  gBenchmark->Show("selectWG");

}
