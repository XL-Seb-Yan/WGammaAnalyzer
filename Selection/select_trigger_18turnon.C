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
#include "../Utils/interface/RunLumiRangeMap.hh"
//#include "../Utils/MyTools.hh"      // various helper functions
// C++ tool
#include <algorithm>
#include <map>
#endif

void select_trigger_18turnon(const TString conf="samples.conf", // input file
	      const TString outputDir=".",  // output directory
	      const TString era=""
	       ) {
  gBenchmark->Start("selectWG");
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

  //Photons
  TH1F* hist00lpa = new TH1F("WGamma00lpa",histotitle+", pt_{#gamma}",100,0,2000);//Pass basic ID, pt turn on
  TH1F* hist00lpb = new TH1F("WGamma00lpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist00lpc = new TH1F("WGamma00lpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist00hpa = new TH1F("WGamma00hpa",histotitle+", pt_{#gamma}",100,0,2000);//Pass basic ID, mass turn on
  TH1F* hist00hpb = new TH1F("WGamma00hpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist00hpc = new TH1F("WGamma00hpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist01lpa = new TH1F("WGamma01lpa",histotitle+", pt_{#gamma}",100,0,2000);//Fired 1 HLT, pt turn on
  TH1F* hist01lpb = new TH1F("WGamma01lpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist01lpc = new TH1F("WGamma01lpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist01hpa = new TH1F("WGamma01hpa",histotitle+", pt_{#gamma}",100,0,2000);//Fired 1 HLT, mass turn on
  TH1F* hist01hpb = new TH1F("WGamma01hpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist01hpc = new TH1F("WGamma01hpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist02lpa = new TH1F("WGamma02lpa",histotitle+", pt_{#gamma}",100,0,2000);//Fired 2 HLT, pt turn on
  TH1F* hist02lpb = new TH1F("WGamma02lpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist02lpc = new TH1F("WGamma02lpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist02hpa = new TH1F("WGamma02hpa",histotitle+", pt_{#gamma}",100,0,2000);//Fired 2 HLT, mass turn on
  TH1F* hist02hpb = new TH1F("WGamma02hpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist02hpc = new TH1F("WGamma02hpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist03lpa = new TH1F("WGamma03lpa",histotitle+", pt_{#gamma}",100,0,2000);//Fired 3 HLT, pt turn on
  TH1F* hist03lpb = new TH1F("WGamma03lpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist03lpc = new TH1F("WGamma03lpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist03hpa = new TH1F("WGamma03hpa",histotitle+", pt_{#gamma}",100,0,2000);//Fired 3 HLT, mass turn on
  TH1F* hist03hpb = new TH1F("WGamma03hpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist03hpc = new TH1F("WGamma03hpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist04lpa = new TH1F("WGamma04lpa",histotitle+", pt_{#gamma}",100,0,2000);//Fired all HLT, pt turn on
  TH1F* hist04lpb = new TH1F("WGamma04lpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist04lpc = new TH1F("WGamma04lpc",histotitle+", #m_{j#gamma}",80,0,4000);

  TH1F* hist04hpa = new TH1F("WGamma04hpa",histotitle+", pt_{#gamma}",100,0,2000);//Fired all HLT, mass turn on
  TH1F* hist04hpb = new TH1F("WGamma04hpb",histotitle+", #eta_{#gamma}",50,-1.44,1.44);
  TH1F* hist04hpc = new TH1F("WGamma04hpc",histotitle+", #m_{j#gamma}",80,0,4000);

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
  Int_t lumiBlock = -999;

  //--Photons---
  Int_t ph_N = -99;//for int we need to use this to initialize the container
  std::vector<float> *ph_pt = new std::vector<float>();
  std::vector<float> *ph_eta = new std::vector<float>();
  std::vector<float> *ph_phi = new std::vector<float>();
  std::vector<float> *ph_E = new std::vector<float>();
  std::vector<bool>  *ph_passEleVeto = new std::vector<bool>();
  std::vector<float> *ph_mvaVal = new std::vector<float>();
  //--Jets(AK8)--
  Int_t jetAK8_puppi_N = -99;
  std::vector<bool>  *jetAK8_puppi_IDTight = new std::vector<bool>();
  std::vector<float> *jetAK8_puppi_softdrop_pt = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_eta = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_phi = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_mass = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_E = new std::vector<float>();
  //--Trigger---
  std::map<std::string,bool> *HLT_isFired = new std::map<std::string,bool>();
  
  // loop over samples
  TTree* eventTree = 0;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    CSample* samp = samplev[isam];
    cout<<"begin loop over files"<<endl;
    TStopwatch stopwatch;

    /*
    // Lumi-section
    vector<RunLumiRangeMap*> Lumi_Photon175;
    Int_t nprescales = samp->prescaleJSONv.size();
    for(int iprescale=0; iprescale<nprescales; iprescale++){
      if(samp->prescaletriggernamev[iprescale] != "HLT_Photon175") continue;
      Lumi_Photon175.push_back(new RunLumiRangeMap());
      Lumi_Photon175.back()->addJSONFile(samp->prescaleJSONv[iprescale]);
      cout<<samp->prescaleJSONv[iprescale]<<endl;
    }
    */

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
      eventTree->SetBranchAddress("ph_mvaVal", &ph_mvaVal);                                       
      //--Jets (AK8 PUPPI)
      eventTree->SetBranchAddress("jetAK8_N", &jetAK8_puppi_N);                                  
      eventTree->SetBranchAddress("jetAK8_pt", &jetAK8_puppi_softdrop_pt);             
      eventTree->SetBranchAddress("jetAK8_eta", &jetAK8_puppi_softdrop_eta);            
      eventTree->SetBranchAddress("jetAK8_phi", &jetAK8_puppi_softdrop_phi);           
      eventTree->SetBranchAddress("jetAK8_e", &jetAK8_puppi_softdrop_E);               
      eventTree->SetBranchAddress("jetAK8_softdrop_massCorr", &jetAK8_puppi_softdrop_mass);                                
      eventTree->SetBranchAddress("jetAK8_IDTight", &jetAK8_puppi_IDTight);                      
      eventTree->SetBranchAddress("HLT_isFired", &HLT_isFired);

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	// Get Events
	TBranch *ph_NBranch = (TBranch*)eventTree->GetBranch("ph_N");
	ph_NBranch->GetEntry(ientry);
	TBranch *jetAK8_NBranch = (TBranch*)eventTree->GetBranch("jetAK8_N");
	jetAK8_NBranch->GetEntry(ientry);
	//Only study events contain photons (EM objects here) -- Muon Sample
	if(ph_N < 1 || jetAK8_puppi_N < 1) continue;
	count1++;
	
	ph_pt->clear();               
	ph_eta->clear();            
	ph_phi->clear();
	ph_E->clear();        
	ph_passEleVeto->clear();           
	ph_mvaVal->clear();           
	jetAK8_puppi_softdrop_pt->clear();               
	jetAK8_puppi_softdrop_eta->clear();             
	jetAK8_puppi_softdrop_phi->clear();             
	jetAK8_puppi_softdrop_mass->clear();            
	jetAK8_puppi_softdrop_E->clear();                          
	jetAK8_puppi_IDTight->clear();
	HLT_isFired->clear();
	eventTree->GetEntry(ientry);

	// lower pt threshold base==============================================
	Int_t index_pl = -99;
	for(int i=0; i<ph_pt->size(); i++){
	  if(ph_pt->at(i) < 100) continue;
	  if(ph_passEleVeto->at(i) != true) continue;
	  //if(abs(ph_eta->at(i)) > 2.4) continue;
	  //if(abs(ph_eta->at(i)) > 1.4442 && abs(ph_eta->at(i))< 1.566) continue;
	  if(abs(ph_eta->at(i)) > 1.44) continue; //barrel photon
	  // Photon mvaID
	  // See https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2
	  if(ph_mvaVal->at(i) < -0.02) continue;
	  index_pl = i;
	  break;
	}
	bool haslowp = false;
	TLorentzVector v_pl;
	if(index_pl != -99) {
	  haslowp = true;
	  v_pl.SetPtEtaPhiE(ph_pt->at(index_pl), ph_eta->at(index_pl), ph_phi->at(index_pl), ph_E->at(index_pl));
	}

	Int_t index_jl = -99;
	if(haslowp){
	  Double_t mass_diff = 99999;
	  for(int i=0; i<jetAK8_puppi_softdrop_pt->size(); i++){
	    //if(jetAK8_puppi_IDTightLepVeto->at(i) != true) continue;
	    if(jetAK8_puppi_IDTight->at(i) != true) continue;
	    if(jetAK8_puppi_softdrop_pt->at(i) < 100) continue;
	    if(jetAK8_puppi_softdrop_mass->at(i) < 0) continue;
	    TLorentzVector v_temp1;
	    v_temp1.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(i),jetAK8_puppi_softdrop_eta->at(i),jetAK8_puppi_softdrop_phi->at(i),jetAK8_puppi_softdrop_E->at(i));
	    if(v_temp1.DeltaR(v_pl) < 1.1) continue;
	    if(abs(jetAK8_puppi_softdrop_mass->at(i) - 80.379) < mass_diff){
	      mass_diff = abs(jetAK8_puppi_softdrop_mass->at(i) - 80.379);
	      index_jl = i;
	    }
	  }
	}
	bool haslowj = false;
	TLorentzVector v_jl;
	if(index_jl != -99) {
	  haslowj = true;
	  v_jl.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(index_jl), jetAK8_puppi_softdrop_eta->at(index_jl), jetAK8_puppi_softdrop_phi->at(index_jl), jetAK8_puppi_softdrop_E->at(index_jl));
	}
	
	// higer pt threshold base==============================================
	Int_t index_ph = -99;
	for(int i=0; i<ph_pt->size(); i++){
	  if(ph_pt->at(i) < 225) continue;
	  if(ph_passEleVeto->at(i) != true) continue;
	  //if(abs(ph_eta->at(i)) > 2.4) continue;
	  //if(abs(ph_eta->at(i)) > 1.4442 && abs(ph_eta->at(i))< 1.566) continue;
	  if(abs(ph_eta->at(i)) > 1.44) continue; //barrel photon
	  // Photon mvaID
	  // See https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2
	  if(ph_mvaVal->at(i) < -0.02) continue;
	  index_ph = i;
	  break;
	}
	bool hashighp = false;
	TLorentzVector v_ph;
	if(index_ph != -99) {
	  hashighp = true;
	  v_ph.SetPtEtaPhiE(ph_pt->at(index_ph), ph_eta->at(index_ph), ph_phi->at(index_ph), ph_E->at(index_ph));
	}

	Int_t index_jh = -99;
	if(hashighp){
	  Double_t mass_diff = 99999;
	  for(int i=0; i<jetAK8_puppi_softdrop_pt->size(); i++){
	    //if(jetAK8_puppi_IDTightLepVeto->at(i) != true) continue;
	    if(jetAK8_puppi_IDTight->at(i) != true) continue;
	    if(jetAK8_puppi_softdrop_pt->at(i) < 225) continue;
	    if(jetAK8_puppi_softdrop_mass->at(i) < 0) continue;
	    TLorentzVector v_temp1;
	    v_temp1.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(i),jetAK8_puppi_softdrop_eta->at(i),jetAK8_puppi_softdrop_phi->at(i),jetAK8_puppi_softdrop_E->at(i));
	    if(v_temp1.DeltaR(v_ph) < 1.1) continue;
	    if(abs(jetAK8_puppi_softdrop_mass->at(i) - 80.379) < mass_diff){
	      mass_diff = abs(jetAK8_puppi_softdrop_mass->at(i) - 80.379);
	      index_jh = i;
	    }
	  }
	}
	bool hashighj = false;
	TLorentzVector v_jh;
	if(index_jh != -99) {
	  hashighj = true;
	  v_jh.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(index_jh), jetAK8_puppi_softdrop_eta->at(index_jh), jetAK8_puppi_softdrop_phi->at(index_jh), jetAK8_puppi_softdrop_E->at(index_jh));
	}

	if(!(haslowp && haslowj) && !(hashighp && hashighj)) continue;

	double invmassl = -99;
	if(haslowp && haslowj){
	  invmassl = (v_pl+v_jl).M();
	  hist00lpa->Fill(ph_pt->at(index_pl));
	  hist00lpb->Fill(ph_eta->at(index_pl));
	  hist00lpc->Fill(invmassl);
	  count2++;
	}
	double invmassh = -99;
	if(hashighp && hashighj){
	  invmassh = (v_ph+v_jh).M();
	  hist00hpa->Fill(ph_pt->at(index_ph));
	  hist00hpb->Fill(ph_eta->at(index_ph));
	  hist00hpc->Fill(invmassh);
	  count3++;
	}
	
	/*
	//Trigger pre-scale
	Int_t tri_prescale = -99;
	RunLumiRangeMap::RunLumiPairType rl(runnum, lumiBlock);
	for(int iprescale=0; iprescale<nprescales; iprescale++){
	  if(Lumi_Photon175[iprescale]->hasRunLumi(rl)){
	    if(tri_prescale > -99 && tri_prescale != samp->prescalev[iprescale]){
	      cout<<"Error: Multiple pre-scale conflict"<<endl;
	      cout<<tri_prescale<<" and "<<samp->prescalev[iprescale]<<" Run: "<<runnum<<" Lumiblock: "<<lumiBlock<<" Evtnum: "<<evtnum<<endl;
	    }
	    tri_prescale = samp->prescalev[iprescale];
	  }
	}
	if(tri_prescale == -99)
	  cout<<"Error: No prescale found"<<endl;

	Double_t weight = 1/double(tri_prescale);
	hist01pa->Fill(ph_pt->at(index_p),weight);
	hist02pa->Fill(ph_eta->at(index_p),weight);
	hist03pa->Fill(ph_E->at(index_p),weight);
	hist04pa->Fill(ph_Et->at(index_p),weight);
	hist05pa->Fill(ph_hOverE->at(index_p),weight);
	hist06pa->Fill(ph_isoGamma->at(index_p),weight);
	hist07pa->Fill(ph_isoCh->at(index_p),weight);
	hist08pa->Fill(ph_passTightId->at(index_p),weight);
	hist09pa->Fill(ph_passEleVeto->at(index_p),weight);
	hist10pa->Fill(ph_mvaVal->at(index_p),weight);
	hist11pa->Fill(ph_mvaCat->at(index_p),weight);
	*/

	//Trigger decision
	bool passTrig = false;
	for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it) {
	  if (it->first.find("HLT_Photon110EB_TightID_TightIso") != std::string::npos && it->second == 1){
	    if(haslowp && haslowj){
	      hist01lpa->Fill(ph_pt->at(index_pl));
	      hist01lpb->Fill(ph_eta->at(index_pl));
	      hist01lpc->Fill(invmassl);
	    }
	    if(hashighp && hashighj){
	      hist01hpa->Fill(ph_pt->at(index_ph));
	      hist01hpb->Fill(ph_eta->at(index_ph));
	      hist01hpc->Fill(invmassh);
	    }
	    passTrig = true;
	  }
	  if (it->first.find("HLT_Photon120EB_TightID_TightIso") != std::string::npos && it->second == 1){
	    if(haslowp && haslowj){
	      hist02lpa->Fill(ph_pt->at(index_pl));
	      hist02lpb->Fill(ph_eta->at(index_pl));
	      hist02lpc->Fill(invmassl);
	    }
	    if(hashighp && hashighj){
	      hist02hpa->Fill(ph_pt->at(index_ph));
	      hist02hpb->Fill(ph_eta->at(index_ph));
	      hist02hpc->Fill(invmassh);
	    }
	    passTrig = true;
	  }
	  if (it->first.find("HLT_Photon200") != std::string::npos && it->second == 1){
	    if(haslowp && haslowj){
	      hist03lpa->Fill(ph_pt->at(index_pl));
	      hist03lpb->Fill(ph_eta->at(index_pl));
	      hist03lpc->Fill(invmassl);
	    }
	    if(hashighp && hashighj){
	      hist03hpa->Fill(ph_pt->at(index_ph));
	      hist03hpb->Fill(ph_eta->at(index_ph));
	      hist03hpc->Fill(invmassh);
	    }
	    passTrig = true;
	  }
	}
	if (!passTrig) continue;
	if(invmassl < 0 && invmassh < 0) {
	  cout<<"Problem!!!"<<endl;
	  cout<<invmassl<<" "<<invmassh<<endl;
	}
	if(haslowp && haslowj){
	  hist04lpa->Fill(ph_pt->at(index_pl));
	  hist04lpb->Fill(ph_eta->at(index_pl));
	  hist04lpc->Fill(invmassl);
	  count4++;
	}
	if(hashighp && hashighj){
	  hist04hpa->Fill(ph_pt->at(index_ph));
	  hist04hpb->Fill(ph_eta->at(index_ph));
	  hist04hpc->Fill(invmassh);
	  count5++;
	}
      }//end of event loop
      cout<<"Pass low pt / high pt sel: "<<count2<<" "<<count3<<endl;
      cout<<"Pass low pt and fired trigger / high pt and fired trigger: "<<count4<<" "<<count5<<endl;
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
    }//end of file loop
  }//end of sample loop

  TFile *outFile = TFile::Open("Trigger_"+era+"turnon.root", "RECREATE");
  hist00lpa->Write();
  hist00lpb->Write();
  hist00lpc->Write();
  hist00hpa->Write();
  hist00hpb->Write();
  hist00hpc->Write();
  hist01lpa->Write();
  hist01lpb->Write();
  hist01lpc->Write();
  hist01hpa->Write();
  hist01hpb->Write();
  hist01hpc->Write();
  hist02lpa->Write();
  hist02lpb->Write();
  hist02lpc->Write();
  hist02hpa->Write();
  hist02hpb->Write();
  hist02hpc->Write();
  hist03lpa->Write();
  hist03lpb->Write();
  hist03lpc->Write();
  hist03hpa->Write();
  hist03hpb->Write();
  hist03hpc->Write();
  hist04lpa->Write();
  hist04lpb->Write();
  hist04lpc->Write();
  hist04hpa->Write();
  hist04hpb->Write();
  hist04hpc->Write();
  
  outFile->Close();

  /*
  //Tirgger Efficiency
  TEfficiency *Eff11 = new TEfficiency(*hist01pb, *hist01pa);
  TEfficiency *Eff12 = new TEfficiency(*hist01pc, *hist01pa);
  TEfficiency *Eff13 = new TEfficiency(*hist01pd, *hist01pa);
  TEfficiency *Eff21 = new TEfficiency(*hist02pb, *hist02pa);
  TEfficiency *Eff22 = new TEfficiency(*hist02pc, *hist02pa);
  TEfficiency *Eff23 = new TEfficiency(*hist02pd, *hist02pa);
  TEfficiency *Eff31 = new TEfficiency(*hist03pb, *hist03pa);
  TEfficiency *Eff32 = new TEfficiency(*hist03pc, *hist03pa);
  TEfficiency *Eff33 = new TEfficiency(*hist03pd, *hist03pa);

  //Non stacked plots
  TLegend *legend2 = new TLegend(0.62,0.16,0.87,0.26);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  TCanvas *ceff1 = new TCanvas("ceff1","",1200,900);
  ceff1->cd();
  Eff11->SetTitle("Trigger Efficiency vs photon pT ("+era+")");
  Eff11->SetLineWidth(2);
  Eff11->SetLineColor(kSpring-6);
  Eff11->Draw();
  gPad->Update();
  Eff11->SetTitle("Trigger Efficiency on photon pT ("+era+"); pt_{#gamma} (GeV); Efficiency"); 
  Eff11->GetPaintedGraph()->GetYaxis()->SetRangeUser(-0.1,1.1);
  Eff11->Draw("SAME");
  Eff12->SetLineWidth(2);
  Eff12->SetLineColor(2);
  Eff12->Draw("SAME");
  Eff13->SetLineWidth(2);
  Eff13->SetLineColor(4);
  Eff13->Draw("SAME");
  legend2->Clear();
  legend2->AddEntry(Eff11,"HLT_Photon165_HE10","f");
  legend2->AddEntry(Eff12,"HLT_Photon175","f");
  legend2->AddEntry(Eff13,"Combined","f");
  legend2->Draw();
  ceff1->Print("eff_pt.png");
  ceff1->Print("eff_pt.pdf");

  TCanvas *ceff2 = new TCanvas("ceff2","",1200,900);
  ceff2->cd();
  Eff21->SetTitle("Trigger Efficiency vs photon #eta ("+era+")");
  Eff21->SetLineWidth(2);
  Eff21->SetLineColor(kSpring-6);
  Eff21->Draw();
  gPad->Update();
  Eff21->SetTitle("Trigger Efficiency on photon #eta ("+era+"); #eta_{#gamma}; Efficiency"); 
  Eff21->GetPaintedGraph()->GetYaxis()->SetRangeUser(-0.1,1.1);
  Eff21->Draw("SAME");
  Eff22->SetLineWidth(2);
  Eff22->SetLineColor(2);
  Eff22->Draw("SAME");
  Eff23->SetLineWidth(2);
  Eff23->SetLineColor(4);
  Eff23->Draw("SAME");
  legend2->Clear();
  legend2->AddEntry(Eff21,"HLT_Photon165_HE10","f");
  legend2->AddEntry(Eff22,"HLT_Photon175","f");
  legend2->AddEntry(Eff23,"Combined","f");
  legend2->Draw();
  ceff2->Print("eff_eta.png");
  ceff2->Print("eff_eta.pdf");

  TCanvas *ceff3 = new TCanvas("ceff3","",1200,900);
  ceff3->cd();
  Eff31->SetTitle("Trigger Efficiency vs M_{j#gamma} ("+era+")");
  Eff31->SetLineWidth(2);
  Eff31->SetLineColor(kSpring-6);
  Eff31->Draw();
  gPad->Update();
  Eff31->SetTitle("Trigger Efficiency on M_{j#gamma} ("+era+"); M_{j#gamma} (GeV); Efficiency"); 
  Eff31->GetPaintedGraph()->GetYaxis()->SetRangeUser(-0.1,1.1);
  Eff31->Draw("SAME");
  Eff32->SetLineWidth(2);
  Eff32->SetLineColor(2);
  Eff32->Draw("SAME");
  Eff33->SetLineWidth(2);
  Eff33->SetLineColor(4);
  Eff33->Draw("SAME");
  legend2->Clear();
  legend2->AddEntry(Eff31,"HLT_Photon165_HE10","f");
  legend2->AddEntry(Eff32,"HLT_Photon175","f");
  legend2->AddEntry(Eff33,"Combined","f");
  legend2->Draw();
  ceff3->Print("eff_m.png");
  ceff3->Print("eff_m.pdf");
  */

  gBenchmark->Show("selectWG");

}
