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
#include "Utils/interface/ConfParse.hh"             // input conf file parser
#include "Utils/interface/CSample.hh"      // helper class to handle samples
#include "Utils/interface/RunLumiRangeMap.hh"
//#include "../Utils/MyTools.hh"      // various helper functions
// C++ tool
#include <algorithm>
#include <map>
#endif

#define runmode 0 //11: ph_corr_up 12: ph_corr_down 21: jec_corr_up 22:jec_corr_down 31: jer_sf_up 32: jer_sf_down

void select2016_Final(const TString conf="samples.conf", // input file
	          const TString outputDir=".",  // output directory
	          const Float_t weight=1   // MC weight?
	       ) {
  gBenchmark->Start("selectWG");
  gROOT->SetBatch(1);
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //=============================================================================================================
  TString type = "SIGNAL"; //DATA, SIGNAL, BKG, CONTROL
  TString formate = "DIVIDED"; //DIVIDED RANGE
  TString histotitle;
  if(type == "SIGNAL"){
    histotitle = "MC(2017) W#gamma";
  }
  else if(type == "BKG"){
    histotitle = "MC(2017) QCD Background";
  }
  else if(type == "DATA"){
    histotitle = "DATA (2017 Single Photon)";
  }
  else if(type == "CONTROL"){
    histotitle = "DATA (2017 Single Muon)";
  }
  else
    cout<<"Wrong data type"<<endl;
  bool isPlot = false;

  // ------------------------Photons------------------------
  TH1F* hist01pa = new TH1F("WGamma01pa",histotitle+", pt_{#gamma}",48,0,2400); // Pass Trigger + Photon MVAID + Basic jet selection
  TH1F* hist02pa = new TH1F("WGamma02pa",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pa = new TH1F("WGamma03pa",histotitle+", #varphi_{#gamma}",50,-3.14,3.14);
  TH1F* hist04pa = new TH1F("WGamma04pa",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist05pa = new TH1F("WGamma05pa",histotitle+", cos(#theta^{*})_{#gamma}",50,0,1);
  TH1F* hist06pa = new TH1F("WGamma06pa",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist07pa = new TH1F("WGamma07pa",histotitle+", pt/M",50,0,3);
  TH1F *hist08pa = new TH1F("WGamma08pa",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pa = new TH1F("WGamma09pa",histotitle+", Medium Photon ID",10,-1,3);
  TH1F *hist10pa = new TH1F("WGamma10pa",histotitle+", Tight Photon ID",10,-1,3);
  TH1F *hist11pa = new TH1F("WGamma11pa",histotitle+", Photon mvaID Value",50,-1.5,1.5);
  TH1F *hist12pa = new TH1F("WGamma12pa",histotitle+", Photon mvaID Category",10,-1,3);
  TH1F *hist13pa = new TH1F("WGamma13pa",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist14pa = new TH1F("WGamma14pa",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist15pa = new TH1F("WGamma15pa",histotitle+", Eletron Veto",10,-1,3);

  TH1F* hist01pb = new TH1F("WGamma01pb",histotitle+", pt_{#gamma}",48,0,2400); 
  TH1F* hist02pb = new TH1F("WGamma02pb",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pb = new TH1F("WGamma03pb",histotitle+", #varphi_{#gamma}",50,-3.14,3.14);
  TH1F* hist04pb = new TH1F("WGamma04pb",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist05pb = new TH1F("WGamma05pb",histotitle+", cos(#theta^{*})_{#gamma}",50,0,1);
  TH1F* hist06pb = new TH1F("WGamma06pb",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist07pb = new TH1F("WGamma07pb",histotitle+", pt/M",50,0,3);
  TH1F *hist08pb = new TH1F("WGamma08pb",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pb = new TH1F("WGamma09pb",histotitle+", Medium Photon ID",10,-1,3);
  TH1F *hist10pb = new TH1F("WGamma10pb",histotitle+", Tight Photon ID",10,-1,3);
  TH1F *hist11pb = new TH1F("WGamma11pb",histotitle+", Photon mvaID Value",50,-1.5,1.5);
  TH1F *hist12pb = new TH1F("WGamma12pb",histotitle+", Photon mvaID Category",10,-1,3);
  TH1F *hist13pb = new TH1F("WGamma13pb",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist14pb = new TH1F("WGamma14pb",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist15pb = new TH1F("WGamma15pb",histotitle+", Eletron Veto",10,-1,3);
  //Jets
  TH1F* hist01ja = new TH1F("WGamma01ja",histotitle+", pt_{j} (AK8 Puppi Softdrop)",50,0,2500); // Pass Trigger + Photon MVAID + Basic jet selection
  TH1F* hist02ja = new TH1F("WGamma02ja",histotitle+", #eta_{j} (AK8 Puppi Softdrop)",50,-5,5);
  TH1F* hist03ja = new TH1F("WGamma03ja",histotitle+", #varphi_{j} (AK8 Puppi Softdrop)",50,-3.14,3.14);
  TH1F* hist04ja = new TH1F("WGamma04ja",histotitle+", E_{j} (AK8 Puppi Softdrop)",50,0,3000);
  TH1F* hist05ja = new TH1F("WGamma05ja",histotitle+", m_{j} (AK8 Puppi Softdrop Corrected)",60,50,110);
  //TH1F* hist05ja = new TH1F("WGamma05ja",histotitle+", m_{j} (AK8 Puppi)",50,0,600);
  TH1F* hist06ja = new TH1F("WGamma06ja",histotitle+", #tau1 (AK8 Puppi Softdrop)",50,0,1);
  TH1F* hist07ja = new TH1F("WGamma07ja",histotitle+", #tau2 (AK8 Puppi Softdrop)",50,0,1);
  TH1F* hist08ja = new TH1F("WGamma08ja",histotitle+", #tau2/#tau1 (AK8 Puppi Softdrop)",50,0,1);
  //TH1F* hist09ja = new TH1F("WGamma09ja",histotitle+", Jet Tight Lepton Veto (AK8 Puppi)",10,-1,3);
  TH1F* hist10ja = new TH1F("WGamma10ja",histotitle+", Jet Tight ID (AK8 Puppi Softdrop)",10,-1,3);
  TH1F* hist11ja = new TH1F("WGamma11ja",histotitle+", Jet mass difference to W (AK8 Puppi Softdrop)",50,0,15);

  TH1F* hist01jb = new TH1F("WGamma01jb",histotitle+", pt_{j} (AK8 Puppi Softdrop)",50,0,2500);
  TH1F* hist02jb = new TH1F("WGamma02jb",histotitle+", #eta_{j} (AK8 Puppi Softdrop)",50,-5,5);
  TH1F* hist03jb = new TH1F("WGamma03jb",histotitle+", #varphi_{j} (AK8 Puppi Softdrop)",50,-3.14,3.14);
  TH1F* hist04jb = new TH1F("WGamma04jb",histotitle+", E_{j} (AK8 Puppi Softdrop)",50,0,3000);
  TH1F* hist05jb = new TH1F("WGamma05jb",histotitle+", m_{j} (AK8 Puppi Softdrop Corrected)",60,50,110);
  //TH1F* hist05jb = new TH1F("WGamma05jb",histotitle+", m_{j} (AK8 Puppi)",50,0,600);
  TH1F* hist06jb = new TH1F("WGamma06jb",histotitle+", #tau1 (AK8 Puppi Softdrop)",50,0,1);
  TH1F* hist07jb = new TH1F("WGamma07jb",histotitle+", #tau2 (AK8 Puppi Softdrop)",50,0,1);
  TH1F* hist08jb = new TH1F("WGamma08jb",histotitle+", #tau2/#tau1 (AK8 Puppi Softdrop)",50,0,1);
  //TH1F* hist09jb = new TH1F("WGamma09jb",histotitle+", Jet Tight Lepton Veto (AK8 Puppi)",10,-1,3);
  TH1F* hist10jb = new TH1F("WGamma10jb",histotitle+", Jet Tight ID (AK8 Puppi Softdrop)",10,-1,3);
  TH1F* hist11jb = new TH1F("WGamma11jb",histotitle+", Jet mass difference to W (AK8 Puppi Softdrop)",50,0,15);
  //WGamma system
  TH1F* hist01sa = new TH1F("WGamma01sa",histotitle+", invariant mass",50,200,2200);
  TH1F* hist01sb = new TH1F("WGamma01sb",histotitle+", invariant mass",50,200,2200);
  TH1F* hist02sa = new TH1F("WGamma02sa",histotitle+", photon jet seperation",50,0,10);
  TH1F* hist02sb = new TH1F("WGamma02sb",histotitle+", photon jet seperation",50,0,10);
  
  //Pileup histogram
  TH1F* histpileup = NULL;
  
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
  float photon_pt, photon_eta, photon_phi, photon_e, photon_corr, photon_energyscale, photon_energyscale_up, photon_energyscale_down, photon_resolution;
  int photon_looseid;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21, ak8puppijet_massdiff;
  float ak8puppijet_jec, ak8puppijet_jec_up, ak8puppijet_jec_down, ak8puppijet_jer, ak8puppijet_jer_sf, ak8puppijet_jer_sf_up, ak8puppijet_jer_sf_down;
  // System
  float sys_costhetastar, sys_ptoverm, sys_invmass, sys_seperation;
  // MC xsec weight
  float xsec_weight;
  float event_pdfunc;
  // Pileup
  int PV_N;
  // evt info
  int run_num, evt_num, lumi_block;  
  
  // Data structures to store info from produced flat ntuples
  //--Event--
  Int_t runnum = -999;
  Int_t evtnum = -999;
  Int_t lumiBlock = -999;
  Int_t pv_N = -99;	
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
  std::vector<float> *ph_corr = new std::vector<float>();
  std::vector<float> *ph_energyscale = new std::vector<float>();
  std::vector<float> *ph_energyscale_up = new std::vector<float>();
  std::vector<float> *ph_energyscale_down = new std::vector<float>();
  std::vector<float> *ph_resolution = new std::vector<float>();
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
  std::vector<float> *jetAK8_puppi_jec = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_jec_up = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_jec_down = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_jer = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_jer_sf = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_jer_sf_up = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_jer_sf_down = new std::vector<float>();
  //--Trigger
  std::map<std::string,bool> *HLT_isFired = new std::map<std::string,bool>();
  //--GenInfo
  float evt_pdfunc = -99;			 
  //--GenParticle
  // std::vector<float> *genPart_pt = new std::vector<float>();
  // std::vector<float> *genPart_eta = new std::vector<float>();
  // std::vector<float> *genPart_phi = new std::vector<float>();
  // std::vector<float> *genPart_mass = new std::vector<float>();
  // std::vector<int>   *genPart_pdgId = new std::vector<int>();
  // std::vector<int>   *genPart_status = new std::vector<int>();
  // std::vector<vector<int>>   *genPart_mother = new std::vector<vector<int>>();
  // std::vector<vector<float>> *genPart_mother_pt = new std::vector<vector<float>>();
  // std::vector<vector<float>> *genPart_mother_eta = new std::vector<vector<float>>();
  // std::vector<vector<float>> *genPart_mother_phi = new std::vector<vector<float>>();
  // std::vector<vector<float>> *genPart_mother_e = new std::vector<vector<float>>();
  // std::vector<vector<int>>   *genPart_daughter = new std::vector<vector<int>>();
  // std::vector<vector<float>> *genPart_daughter_pt = new std::vector<vector<float>>();
  // std::vector<vector<float>> *genPart_daughter_eta = new std::vector<vector<float>>();
  // std::vector<vector<float>> *genPart_daughter_phi = new std::vector<vector<float>>();
  // std::vector<vector<float>> *genPart_daughter_e = new std::vector<vector<float>>();

  
  // loop over samples
  TTree* eventTree = 0;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    CSample* samp = samplev[isam];
	
	//pileip
	histpileup = new TH1F("pileup","Pile up",50,0,100);

    //Set up output ntuple
    //TString outfilename1 = ntupDir + TString("/") + snamev[isam] + TString("_WGamma_full_full_weightedTo41p54.root");
#if runmode == 11
    TString outfilename1 = ntupDir + TString("/") + snamev[isam] + TString("_photonSFup_WGamma_full_full_May22.root");
#elif runmode == 12
    TString outfilename1 = ntupDir + TString("/") + snamev[isam] + TString("_photonSFdown_WGamma_full_full_May22.root");
#elif runmode == 21
    TString outfilename1 = ntupDir + TString("/") + snamev[isam] + TString("_jetjecup_WGamma_full_full_May22.root");
#elif runmode == 22
    TString outfilename1 = ntupDir + TString("/") + snamev[isam] + TString("_jetjecdown_WGamma_full_full_May22.root");
#elif runmode == 31
    TString outfilename1 = ntupDir + TString("/") + snamev[isam] + TString("_jetjerup_WGamma_full_full_May22.root");
#elif runmode == 32
	TString outfilename1 = ntupDir + TString("/") + snamev[isam] + TString("_jetjerdown_WGamma_full_full_May22.root");
#else
	TString outfilename1 = ntupDir + TString("/") + snamev[isam] + TString("_nominal_pileup_WGamma_full_full_May22.root");
#endif
    TFile *outFile1 = new TFile(outfilename1,"RECREATE"); 
    TTree *outTree1 = new TTree("Events","Events");
	outTree1->Branch("PV_N",            &PV_N,          "PV_N/I");
    outTree1->Branch("photon_pt",       &photon_pt,      "photon_pt/F");
    outTree1->Branch("photon_eta",      &photon_eta,      "photon_eta/F");
    outTree1->Branch("photon_phi",      &photon_phi,      "photon_phi/F");
    outTree1->Branch("photon_e",        &photon_e,      "photon_e/F");
    outTree1->Branch("photon_looseid",   &photon_looseid,  "photon_looseid/I");
	outTree1->Branch("photon_corr",   &photon_corr,  "photon_corr/F");
	outTree1->Branch("photon_energyscale",   &photon_energyscale,  "photon_energyscale/F");
	outTree1->Branch("photon_resolution",   &photon_resolution,  "photon_resolution/F");
    outTree1->Branch("ak8puppijet_pt",       &ak8puppijet_pt,      "ak8puppijet_pt/F");
    outTree1->Branch("ak8puppijet_eta",      &ak8puppijet_eta,      "ak8puppijet_eta/F");
    outTree1->Branch("ak8puppijet_phi",      &ak8puppijet_phi,      "ak8puppijet_phi/F");
    outTree1->Branch("ak8puppijet_e",        &ak8puppijet_e,      "ak8puppijet_e/F");
    outTree1->Branch("ak8puppijet_masssoftdropcorr",   &ak8puppijet_masssoftdropcorr,  "ak8puppijet_masssoftdropcorr/F");
    outTree1->Branch("ak8puppijet_tau21",              &ak8puppijet_tau21,             "ak8puppijet_tau21/F");
    outTree1->Branch("ak8puppijet_massdiff",           &ak8puppijet_massdiff,          "ak8puppijet_massdiff/F");
	outTree1->Branch("ak8puppijet_jec",                &ak8puppijet_jec,               "ak8puppijet_jec/F");
    outTree1->Branch("ak8puppijet_jecUp",              &ak8puppijet_jec_up,             "ak8puppijet_jec_up/F");
    outTree1->Branch("ak8puppijet_jecDown",            &ak8puppijet_jec_down,           "ak8puppijet_jec_down/F");
	outTree1->Branch("ak8puppijet_jer",                &ak8puppijet_jer,               "ak8puppijet_jer/F");
	outTree1->Branch("ak8puppijet_jer_sf",             &ak8puppijet_jer_sf,            "ak8puppijet_jer_sf/F");
    outTree1->Branch("ak8puppijet_jer_sf_Up",          &ak8puppijet_jer_sf_up,         "ak8puppijet_jer_sf_up/F");
    outTree1->Branch("ak8puppijet_jer_sf_Down",        &ak8puppijet_jer_sf_down,       "ak8puppijet_jer_sf_down/F");
    outTree1->Branch("sys_costhetastar",        &sys_costhetastar,      "sys_costhetastar/F");
    outTree1->Branch("sys_ptoverm",             &sys_ptoverm,           "sys_ptoverm/F");
    outTree1->Branch("sys_invmass",             &sys_invmass,           "sys_invmass/F");
    outTree1->Branch("sys_seperation",          &sys_seperation,        "sys_seperation/F");
    outTree1->Branch("xsec_weight",             &xsec_weight,           "xsec_weight/F");
    outTree1->Branch("event_pdfunc",            &event_pdfunc,          "event_pdfunc/F");
	outTree1->Branch("run_num",                 &run_num,               "run_num/I");
	outTree1->Branch("evt_num",                 &evt_num,               "evt_num/I");
	outTree1->Branch("lumi_block",              &lumi_block,            "lumi_block/I");
   

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
	  //--Run Info--
      eventTree->SetBranchAddress("EVENT_run", &runnum);                      
      eventTree->SetBranchAddress("EVENT_event", &evtnum);                    
      eventTree->SetBranchAddress("EVENT_lumiBlock", &lumiBlock);
	  //--PV--
	  eventTree->SetBranchAddress("PV_N", &pv_N);     
      //--Photons--
      eventTree->SetBranchAddress("ph_N", &ph_N);                              
      eventTree->SetBranchAddress("ph_pt", &ph_pt);                            
      eventTree->SetBranchAddress("ph_eta", &ph_eta);                           
      eventTree->SetBranchAddress("ph_phi", &ph_phi);                          
      eventTree->SetBranchAddress("ph_e", &ph_E);                             
      //eventTree->SetBranchAddress("ph_hOverE", &ph_hOverE);                   
      //eventTree->SetBranchAddress("ph_isoGamma", &ph_isoGamma);                
      //eventTree->SetBranchAddress("ph_isoCh", &ph_isoCh);                       
      eventTree->SetBranchAddress("ph_passEleVeto", &ph_passEleVeto);          
      eventTree->SetBranchAddress("ph_passLooseId", &ph_passLooseId);           
      //eventTree->SetBranchAddress("ph_passMediumId", &ph_passMediumId);         
      eventTree->SetBranchAddress("ph_passTightId", &ph_passTightId);           
      // eventTree->SetBranchAddress("ph_mvaVal", &ph_mvaVal);                    
      // eventTree->SetBranchAddress("ph_mvaCat", &ph_mvaCat);  
      eventTree->SetBranchAddress("ph_Corr", &ph_corr);                    
      eventTree->SetBranchAddress("ph_energyscale", &ph_energyscale); 	 
	  eventTree->SetBranchAddress("ph_energyscale_up", &ph_energyscale_up); 	
	  eventTree->SetBranchAddress("ph_energyscale_down", &ph_energyscale_down); 		  
	  eventTree->SetBranchAddress("ph_resolution", &ph_resolution); 
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
      eventTree->SetBranchAddress("jetAK8_softdrop_jec", &jetAK8_puppi_jec);
      eventTree->SetBranchAddress("jetAK8_softdrop_jecUp", &jetAK8_puppi_jec_up);
      eventTree->SetBranchAddress("jetAK8_softdrop_jecDown", &jetAK8_puppi_jec_down);
	  eventTree->SetBranchAddress("jetAK8_jer_sigma_pt", &jetAK8_puppi_jer);
      eventTree->SetBranchAddress("jetAK8_jer_sf", &jetAK8_puppi_jer_sf);
      eventTree->SetBranchAddress("jetAK8_jer_sf_up", &jetAK8_puppi_jer_sf_up);
      eventTree->SetBranchAddress("jetAK8_jer_sf_down", &jetAK8_puppi_jer_sf_down);
      eventTree->SetBranchAddress("HLT_isFired", &HLT_isFired);
	  //--GenInfo
	  eventTree->SetBranchAddress("PDF_rms", &evt_pdfunc);
      //--GenPartticle
      /*
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
      */

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<2; ientry++) {
	count1++;

	// Get Events
	ph_pt->clear();               
	ph_eta->clear();            
	ph_phi->clear();              
	ph_E->clear();                
	//ph_hOverE->clear();            
	//ph_isoGamma->clear();          
	//ph_isoCh->clear();             
	ph_passEleVeto->clear();       
	ph_passLooseId->clear();       
	//ph_passMediumId->clear();     
	ph_passTightId->clear();      
	// ph_mvaVal->clear();           
	// ph_mvaCat->clear(); 
	ph_corr->clear();     
	ph_energyscale->clear();     
	ph_energyscale_up->clear();     
	ph_energyscale_down->clear();     
	ph_resolution->clear();     
	jetAK8_puppi_softdrop_pt->clear();               
	jetAK8_puppi_softdrop_eta->clear();             
	jetAK8_puppi_softdrop_phi->clear();             
	jetAK8_puppi_softdrop_mass->clear();            
	jetAK8_puppi_softdrop_E->clear();               
	jetAK8_puppi_tau1->clear();             
	jetAK8_puppi_tau2->clear();            
	jetAK8_puppi_IDTight->clear();

	jetAK8_puppi_jec->clear();
	jetAK8_puppi_jec_up->clear();
	jetAK8_puppi_jec_down->clear();
	jetAK8_puppi_jer->clear();
	jetAK8_puppi_jer_sf->clear();
	jetAK8_puppi_jer_sf_up->clear();
	jetAK8_puppi_jer_sf_down->clear();
	HLT_isFired->clear();
	/*
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
	*/
	
	//-------------------------------------pileup----------------------------------------------
	// full pile-up before any selection, used to creat histogram for pileup reweighting after full selection
	TBranch *PileupBranch = (TBranch*)eventTree->GetBranch("PV_N");
	PileupBranch->GetEntry(ientry);
	histpileup->Fill(pv_N);
	
	//-------------------------------------Trigger----------------------------------------------
	// HLT Trigger decision first, improve speed
	TBranch *HLTBranch = (TBranch*)eventTree->GetBranch("HLT_isFired");
	HLTBranch->GetEntry(ientry);
	
	bool passTrig = false;
	for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it) {
	  if (it->first.find("HLT_Photon165_HE10_v") != std::string::npos && it->second == 1)
	    passTrig = true;
	  // if (it->first.find("HLT_Photon165_R9Id90_HE10_IsoM_v") != std::string::npos && it->second == 1)
	    // passTrig = true;
	  if (it->first.find("HLT_Photon175_v") != std::string::npos && it->second == 1)
	    passTrig = true;
	}
	if (!passTrig) continue;
	count0++;
	eventTree->GetEntry(ientry);

	//Only study events contain photons (EM objects here) and jets -- DATA -- 1st skimming
	if(ph_N  == 0 || jetAK8_puppi_N  == 0) continue;
	count2++;
	
	//-----------------------------------Photon-----------------------------------------------------
	// Photon selection Flags
	Bool_t pass_p1 = false;
	Bool_t pass_p2 = false;

	// Locate the first photon in EM objects array (highest pt) and kinetic cut
	Int_t index_p = -99;
	for(int i=0; i<ph_pt->size(); i++){
#if runmode == 11
	  if(ph_pt->at(i)*(ph_energyscale_up->at(i)/ph_E->at(i)) < 225) continue;
#elif runmode == 12
      if(ph_pt->at(i)*(ph_energyscale_down->at(i)/ph_E->at(i)) < 225) continue;
#else
	  if(ph_pt->at(i)*ph_corr->at(i) < 225) continue;
#endif
	  if(ph_passEleVeto->at(i) != true) continue;
	  if(abs(ph_eta->at(i)) > 1.44) continue; //barrel photon
	  if(!(ph_passLooseId->at(i))) continue;
	  // Photon mvaID
	  // See https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2
	  // if(ph_mvaVal->at(i) < -0.02) continue;
	  index_p = i;
	  break;
	} 
	if(index_p == -99) continue;
	count3++;
	pass_p1 = true;
	pass_p2 = true;
	    

	// Assigning leading photon
	TLorentzVector v_p;
#if runmode == 11
	  v_p.SetPtEtaPhiE(ph_pt->at(index_p)*(ph_energyscale_up->at(index_p)/ph_E->at(index_p)), ph_eta->at(index_p), ph_phi->at(index_p), ph_E->at(index_p)*(ph_energyscale_up->at(index_p)/ph_E->at(index_p)));
#elif runmode == 12
      v_p.SetPtEtaPhiE(ph_pt->at(index_p)*(ph_energyscale_down->at(index_p)/ph_E->at(index_p)), ph_eta->at(index_p), ph_phi->at(index_p), ph_E->at(index_p)*(ph_energyscale_down->at(index_p)/ph_E->at(index_p)));
#else
	  v_p.SetPtEtaPhiE(ph_pt->at(index_p)*ph_corr->at(index_p), ph_eta->at(index_p), ph_phi->at(index_p), ph_E->at(index_p)*ph_corr->at(index_p));
#endif
	

	//--------------------------------------Jet----------------------------------------------------

	// Jet selection Flags
	Bool_t pass_j1 = false;
	Bool_t pass_j2 = false;
	
	// Selection
	Int_t index_j = -99;
	Double_t mass_diff = 99999;
#if runmode == 21
	for(int i=0; i<jetAK8_puppi_softdrop_pt->size(); i++){
	  if(jetAK8_puppi_IDTight->at(i) != true) continue;
	  if(jetAK8_puppi_softdrop_pt->at(i)*jetAK8_puppi_jec_up->at(i) < 225) continue;
	  if(jetAK8_puppi_softdrop_mass->at(i)*(jetAK8_puppi_jec_up->at(i)/jetAK8_puppi_jec->at(i)) < 0) continue; //mass is already corrected
	  TLorentzVector v_temp1;
	  v_temp1.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(i)*jetAK8_puppi_jec_up->at(i),jetAK8_puppi_softdrop_eta->at(i),jetAK8_puppi_softdrop_phi->at(i),jetAK8_puppi_softdrop_E->at(i)*jetAK8_puppi_jec_up->at(i));
	  if(v_temp1.DeltaR(v_p) < 1.1) continue;
	  if(abs(jetAK8_puppi_softdrop_mass->at(i)*(jetAK8_puppi_jec_up->at(i)/jetAK8_puppi_jec->at(i)) - 80.379) < mass_diff){
	    mass_diff = abs(jetAK8_puppi_softdrop_mass->at(i)*(jetAK8_puppi_jec_up->at(i)/jetAK8_puppi_jec->at(i)) - 80.379);
	    index_j = i;
	  }
	}
#elif runmode == 22
	for(int i=0; i<jetAK8_puppi_softdrop_pt->size(); i++){
	  if(jetAK8_puppi_IDTight->at(i) != true) continue;
	  if(jetAK8_puppi_softdrop_pt->at(i)*jetAK8_puppi_jec_down->at(i) < 225) continue;
	  if(jetAK8_puppi_softdrop_mass->at(i)*(jetAK8_puppi_jec_down->at(i)/jetAK8_puppi_jec->at(i)) < 0) continue; //mass is already corrected
	  TLorentzVector v_temp1;
	  v_temp1.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(i)*jetAK8_puppi_jec_down->at(i),jetAK8_puppi_softdrop_eta->at(i),jetAK8_puppi_softdrop_phi->at(i),jetAK8_puppi_softdrop_E->at(i)*jetAK8_puppi_jec_down->at(i));
	  if(v_temp1.DeltaR(v_p) < 1.1) continue;
	  if(abs(jetAK8_puppi_softdrop_mass->at(i)*(jetAK8_puppi_jec_down->at(i)/jetAK8_puppi_jec->at(i)) - 80.379) < mass_diff){
	    mass_diff = abs(jetAK8_puppi_softdrop_mass->at(i)*(jetAK8_puppi_jec_down->at(i)/jetAK8_puppi_jec->at(i)) - 80.379);
	    index_j = i;
	  }
	}
#else
	for(int i=0; i<jetAK8_puppi_softdrop_pt->size(); i++){
	  if(jetAK8_puppi_IDTight->at(i) != true) continue;
	  if(jetAK8_puppi_softdrop_pt->at(i)*jetAK8_puppi_jec->at(i) < 225) continue;
	  if(jetAK8_puppi_softdrop_mass->at(i) < 0) continue; //mass is already corrected
	  TLorentzVector v_temp1;
	  v_temp1.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(i)*jetAK8_puppi_jec->at(i),jetAK8_puppi_softdrop_eta->at(i),jetAK8_puppi_softdrop_phi->at(i),jetAK8_puppi_softdrop_E->at(i)*jetAK8_puppi_jec->at(i));
	  if(v_temp1.DeltaR(v_p) < 1.1) continue;
	  if(abs(jetAK8_puppi_softdrop_mass->at(i) - 80.379) < mass_diff){
	    mass_diff = abs(jetAK8_puppi_softdrop_mass->at(i) - 80.379);
	    index_j = i;
	  }
	}
#endif
	if(index_j == -99) continue;
	pass_j1 = true;
	count4++;

	// N-subjettiness
	Double_t tau21 = jetAK8_puppi_tau2->at(index_j) / jetAK8_puppi_tau1->at(index_j);	

	// Assigning candidate jet
	TLorentzVector v_j;
	// Apply JER
	TRandom3 R(ientry);
	double random = R.Gaus(0,jetAK8_puppi_jer->at(index_j));
	double p4scale = 1 + random * sqrt(pow(jetAK8_puppi_jer_sf->at(index_j),2) - 1 > 0? pow(jetAK8_puppi_jer_sf->at(index_j),2) - 1: 0);
	double p4scale_up = 1 + random * sqrt(pow(jetAK8_puppi_jer_sf_up->at(index_j),2) - 1 > 0? pow(jetAK8_puppi_jer_sf_up->at(index_j),2) - 1: 0);
	double p4scale_down = 1 + random * sqrt(pow(jetAK8_puppi_jer_sf_down->at(index_j),2) - 1 > 0? pow(jetAK8_puppi_jer_sf_down->at(index_j),2) - 1: 0);
		
#if runmode == 21
	v_j.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(index_j)*p4scale*jetAK8_puppi_jec_up->at(index_j), jetAK8_puppi_softdrop_eta->at(index_j), jetAK8_puppi_softdrop_phi->at(index_j), jetAK8_puppi_softdrop_E->at(index_j)*p4scale*jetAK8_puppi_jec_up->at(index_j));
#elif runmode == 22
    v_j.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(index_j)*p4scale*jetAK8_puppi_jec_down->at(index_j), jetAK8_puppi_softdrop_eta->at(index_j), jetAK8_puppi_softdrop_phi->at(index_j), jetAK8_puppi_softdrop_E->at(index_j)*p4scale*jetAK8_puppi_jec_down->at(index_j));
#elif runmode == 31
    v_j.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale_up, jetAK8_puppi_softdrop_eta->at(index_j), jetAK8_puppi_softdrop_phi->at(index_j), jetAK8_puppi_softdrop_E->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale_up);
#elif runmode == 32
    v_j.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale_down, jetAK8_puppi_softdrop_eta->at(index_j), jetAK8_puppi_softdrop_phi->at(index_j), jetAK8_puppi_softdrop_E->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale_down);
#else
	v_j.SetPtEtaPhiE(jetAK8_puppi_softdrop_pt->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale, jetAK8_puppi_softdrop_eta->at(index_j), jetAK8_puppi_softdrop_phi->at(index_j), jetAK8_puppi_softdrop_E->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale);
#endif
	// Tight Selection
	// Tau21: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetWtagging
	// Bool_t tau21_cut = (tau21 < 0.45);
	//Bool_t mass_range = (jetAK8_puppi_softdrop_mass->at(index_j) > 50 && jetAK8_puppi_softdrop_mass->at(index_j) < 105);
	//Bool_t mass_range = (jetAK8_puppi_softdrop_mass->at(index_j) > 65 && jetAK8_puppi_softdrop_mass->at(index_j) < 95); //for sig samples derived from data for TMVA
	if(v_j.Pt() < 225) continue;
	pass_j2 = true;

	//-----------------------------------------W+Gamma system--------------------------------------------------
	bool pass_s1 = false;
	TLorentzVector v_sys;
	v_sys = v_j + v_p;
	
	// Calculate invariant mass and pt/m
	Double_t invmass = 0;
	Double_t ptoverM = 0;
	invmass = (v_sys).M();
#if runmode == 11
	  ptoverM = ph_pt->at(index_p)*(ph_energyscale_up->at(index_p)/ph_E->at(index_p))/invmass;
#elif runmode == 12
	  ptoverM = ph_pt->at(index_p)*(ph_energyscale_down->at(index_p)/ph_E->at(index_p))/invmass;
#else
	  ptoverM = ph_pt->at(index_p)*ph_corr->at(index_p)/invmass;
#endif
	
	//if(invmass > 675 && invmass < 1125)
	pass_s1 = true;

	// Calculate cos(theta*)
	Double_t cosThetaStar = -999;
	TLorentzVector v_boosted_j, v_boosted_p;
	v_boosted_p = v_p;
    v_boosted_p.Boost(-(v_sys.BoostVector()));
    cosThetaStar = abs(v_boosted_p.Pz()/v_boosted_p.P());

	// Calculate seperation
	Double_t seperation = v_j.DeltaR(v_p);

	/*
	// Passing p1 and j1 selection
	if(pass_p1 && pass_j1){
	  hist01pa->Fill(ph_pt->at(index_p));
	  hist02pa->Fill(ph_eta->at(index_p));
	  hist03pa->Fill(ph_phi->at(index_p));
	  hist04pa->Fill(ph_E->at(index_p));
	  hist05pa->Fill(cosThetaStar);
	  hist06pa->Fill(ph_hOverE->at(index_p));
	  hist07pa->Fill(ptoverM);
	  hist08pa->Fill(ph_passLooseId->at(index_p));
	  hist09pa->Fill(ph_passMediumId->at(index_p));
	  hist10pa->Fill(ph_passTightId->at(index_p));
	  hist11pa->Fill(ph_mvaVal->at(index_p));
	  hist12pa->Fill(ph_mvaCat->at(index_p));
	  hist13pa->Fill(ph_isoGamma->at(index_p));
	  hist14pa->Fill(ph_isoCh->at(index_p));
	  hist15pa->Fill(ph_passEleVeto->at(index_p));
	  hist01ja->Fill(jetAK8_puppi_softdrop_pt->at(index_j));
	  hist02ja->Fill(jetAK8_puppi_softdrop_eta->at(index_j));
	  hist03ja->Fill(jetAK8_puppi_softdrop_phi->at(index_j));
	  hist04ja->Fill(jetAK8_puppi_softdrop_E->at(index_j));
	  hist05ja->Fill(jetAK8_puppi_softdrop_mass->at(index_j));
	  hist06ja->Fill(jetAK8_puppi_tau1->at(index_j));
	  hist07ja->Fill(jetAK8_puppi_tau2->at(index_j));
	  hist08ja->Fill(tau21);
	  //hist09ja->Fill(jetAK8_puppi_IDTightLepVeto->at(index_j));
	  hist10ja->Fill(jetAK8_puppi_IDTight->at(index_j));
	  hist11ja->Fill(mass_diff);
	  hist01sa->Fill(invmass);
	  hist02sa->Fill(seperation);
	}
	if(pass_p2 && pass_j2){
	  hist01pb->Fill(ph_pt->at(index_p));
	  hist02pb->Fill(ph_eta->at(index_p));
	  hist03pb->Fill(ph_phi->at(index_p));
	  hist04pb->Fill(ph_E->at(index_p));
	  hist05pb->Fill(cosThetaStar);
	  hist06pb->Fill(ph_hOverE->at(index_p));
	  hist07pb->Fill(ptoverM);
	  hist08pb->Fill(ph_passLooseId->at(index_p));
	  hist09pb->Fill(ph_passMediumId->at(index_p));
	  hist10pb->Fill(ph_passTightId->at(index_p));
	  hist11pb->Fill(ph_mvaVal->at(index_p));
	  hist12pb->Fill(ph_mvaCat->at(index_p));
	  hist13pb->Fill(ph_isoGamma->at(index_p));
	  hist14pb->Fill(ph_isoCh->at(index_p));
	  hist15pb->Fill(ph_passEleVeto->at(index_p));
	  hist01jb->Fill(jetAK8_puppi_softdrop_pt->at(index_j));
	  hist02jb->Fill(jetAK8_puppi_softdrop_eta->at(index_j));
	  hist03jb->Fill(jetAK8_puppi_softdrop_phi->at(index_j));
	  hist04jb->Fill(jetAK8_puppi_softdrop_E->at(index_j));
	  hist05jb->Fill(jetAK8_puppi_softdrop_mass->at(index_j));
	  hist06jb->Fill(jetAK8_puppi_tau1->at(index_j));
	  hist07jb->Fill(jetAK8_puppi_tau2->at(index_j));
	  hist08jb->Fill(tau21);
	  //hist09jb->Fill(jetAK8_puppi_IDTightLepVeto->at(index_j));
	  hist10jb->Fill(jetAK8_puppi_IDTight->at(index_j));
	  hist11jb->Fill(mass_diff);
	  hist01sb->Fill(invmass);
	  hist02sb->Fill(seperation);
	}
	*/
	if(pass_p2 && pass_j2 && pass_s1){
#if runmode == 11
	  photon_pt = ph_pt->at(index_p)*(ph_energyscale_up->at(index_p)/ph_E->at(index_p));
      photon_e = ph_E->at(index_p)*(ph_energyscale_up->at(index_p)/ph_E->at(index_p));
#elif runmode == 12
	  photon_pt = ph_pt->at(index_p)*(ph_energyscale_down->at(index_p)/ph_E->at(index_p));
      photon_e = ph_E->at(index_p)*(ph_energyscale_down->at(index_p)/ph_E->at(index_p));
#else
	  photon_pt = ph_pt->at(index_p)*ph_corr->at(index_p);
      photon_e = ph_E->at(index_p)*ph_corr->at(index_p);
#endif  
	  photon_eta = ph_eta->at(index_p);
	  photon_phi = ph_phi->at(index_p);
	  // photon_mvaval = ph_mvaVal->at(index_p);
	  // photon_mvacat = ph_mvaCat->at(index_p);
	  photon_looseid = ph_passLooseId->at(index_p);
	  photon_corr = ph_corr->at(index_p);
	  photon_energyscale = ph_energyscale->at(index_p);
	  photon_resolution = ph_resolution->at(index_p);
#if runmode == 21
	  ak8puppijet_pt = jetAK8_puppi_softdrop_pt->at(index_j)*jetAK8_puppi_jec_up->at(index_j)*p4scale;
	  ak8puppijet_e = jetAK8_puppi_softdrop_E->at(index_j)*jetAK8_puppi_jec_up->at(index_j)*p4scale;
	  ak8puppijet_masssoftdropcorr = jetAK8_puppi_softdrop_mass->at(index_j)*jetAK8_puppi_jec_up->at(index_j)/jetAK8_puppi_jec->at(index_j)*p4scale;
#elif runmode == 22
	  ak8puppijet_pt = jetAK8_puppi_softdrop_pt->at(index_j)*jetAK8_puppi_jec_down->at(index_j)*p4scale;
	  ak8puppijet_e = jetAK8_puppi_softdrop_E->at(index_j)*jetAK8_puppi_jec_down->at(index_j)*p4scale;
	  ak8puppijet_masssoftdropcorr = jetAK8_puppi_softdrop_mass->at(index_j)*jetAK8_puppi_jec_down->at(index_j)/jetAK8_puppi_jec->at(index_j)*p4scale;
#elif runmode == 31
	  ak8puppijet_pt = jetAK8_puppi_softdrop_pt->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale_up;
	  ak8puppijet_e = jetAK8_puppi_softdrop_E->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale_up;
	  ak8puppijet_masssoftdropcorr = jetAK8_puppi_softdrop_mass->at(index_j)*p4scale_up;
#elif runmode == 32
	  ak8puppijet_pt = jetAK8_puppi_softdrop_pt->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale_down;
	  ak8puppijet_e = jetAK8_puppi_softdrop_E->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale_down;
	  ak8puppijet_masssoftdropcorr = jetAK8_puppi_softdrop_mass->at(index_j)*p4scale_down;
#else
	  ak8puppijet_pt = jetAK8_puppi_softdrop_pt->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale;
	  ak8puppijet_e = jetAK8_puppi_softdrop_E->at(index_j)*jetAK8_puppi_jec->at(index_j)*p4scale;
	  ak8puppijet_masssoftdropcorr = jetAK8_puppi_softdrop_mass->at(index_j)*p4scale;
#endif
	  ak8puppijet_eta = jetAK8_puppi_softdrop_eta->at(index_j);
	  ak8puppijet_phi = jetAK8_puppi_softdrop_phi->at(index_j);
	  ak8puppijet_tau21 = tau21;
	  ak8puppijet_massdiff = mass_diff;
	  ak8puppijet_jec = jetAK8_puppi_jec->at(index_j);
	  ak8puppijet_jec_up = jetAK8_puppi_jec_up->at(index_j);
	  ak8puppijet_jec_down = jetAK8_puppi_jec_down->at(index_j);
	  ak8puppijet_jer = jetAK8_puppi_jer->at(index_j);
	  ak8puppijet_jer_sf = jetAK8_puppi_jer_sf->at(index_j);
	  ak8puppijet_jer_sf_up = jetAK8_puppi_jer_sf_up->at(index_j);
	  ak8puppijet_jer_sf_down = jetAK8_puppi_jer_sf_down->at(index_j);
	  sys_costhetastar = cosThetaStar;
	  sys_ptoverm = ptoverM;
	  sys_invmass = invmass;
	  sys_seperation = seperation;
	  xsec_weight = weight;
	  event_pdfunc = evt_pdfunc;
	  PV_N = pv_N;
	  run_num = runnum;
	  evt_num = evtnum;
	  lumi_block = lumiBlock;
	  outTree1->Fill();
	}
	if(pass_p2) count5++;

	/*
	if(invmass > 1000) continue;

	cout<<"information of selected RECO objects:"<<endl;
	cout<<"AK8jet pt is "<<ak8puppijet_pt<<" AK8jet eta is "<<ak8puppijet_eta<<" AK8jet phi is "<<ak8puppijet_phi<<" AK8jet e is "<<ak8puppijet_e<<endl;
	cout<<"photon pt is "<<photon_pt<<" photon eta is "<<photon_eta<<" photon phi is "<<photon_phi<<" photon e is "<<photon_e<<endl;
	TLorentzVector test_p, test_j;
	test_p.SetPtEtaPhiE(photon_pt,photon_eta,photon_phi,photon_e);
	test_j.SetPtEtaPhiE(ak8puppijet_pt,ak8puppijet_eta,ak8puppijet_phi,ak8puppijet_e);
	cout<<"RECO invariant mass "<<(test_p+test_j).M()<<endl;
	//-----------------------------------GEN level-----------------------------------------------
	cout<<"begin of GEN level event dump:"<<endl;
	for(int i=0; i<genPart_pt->size(); i++){
	  if(abs(genPart_pdgId->at(i)) == 24){
	    cout<<"particle is: "<<genPart_pdgId->at(i)<<" pt is: "<<genPart_pt->at(i)<<" eta is: "<<genPart_eta->at(i)<<" phi is: "<<genPart_phi->at(i)<<" status is: "<<genPart_status->at(i)<<endl;
	    test_j.SetPtEtaPhiM(genPart_pt->at(i),genPart_eta->at(i),genPart_phi->at(i),genPart_mass->at(i));
	    int Ndau = genPart_daughter->at(i).size();
	    for(int j=0; j<Ndau; j++){
	      cout<<"  daughter "<<j<<" is "<<genPart_daughter->at(i).at(j)<<" pt is "<<genPart_daughter_pt->at(i).at(j)<<" eta is "<<genPart_daughter_eta->at(i).at(j)<<" phi is "<<genPart_daughter_phi->at(i).at(j)<<endl;
	    }
	    int Nmom = genPart_mother->at(i).size();
	    for(int j=0; j<Nmom; j++){
	      cout<<"  mother "<<j<<" is "<<genPart_mother->at(i).at(j)<<" pt is "<<genPart_mother_pt->at(i).at(j)<<" eta is "<<genPart_mother_eta->at(i).at(j)<<" phi is "<<genPart_mother_phi->at(i).at(j)<<endl;
	      TLorentzVector H_p;
	      H_p.SetPtEtaPhiE(genPart_mother_pt->at(i).at(j),genPart_mother_eta->at(i).at(j),genPart_mother_phi->at(i).at(j),genPart_mother_e->at(i).at(j));
	      cout<<"  mother mass is "<<H_p.M()<<endl;
	    }
	  }
	}
	cout<<"====================================="<<endl;
	for(int i=0; i<genPart_pt->size(); i++){
	  if(abs(genPart_pdgId->at(i)) == 22){
	    cout<<"particle is: "<<genPart_pdgId->at(i)<<" pt is: "<<genPart_pt->at(i)<<" eta is: "<<genPart_eta->at(i)<<" phi is: "<<genPart_phi->at(i)<<" status is: "<<genPart_status->at(i)<<endl;
	    test_p.SetPtEtaPhiM(genPart_pt->at(i),genPart_eta->at(i),genPart_phi->at(i),genPart_mass->at(i));
	    int Nmom = genPart_mother->at(i).size();
	    for(int j=0; j<Nmom; j++){
	      cout<<"  mother "<<j<<" is "<<genPart_mother->at(i).at(j)<<" pt is "<<genPart_mother_pt->at(i).at(j)<<" eta is "<<genPart_mother_eta->at(i).at(j)<<" phi is "<<genPart_mother_phi->at(i).at(j)<<endl;
	      TLorentzVector H_p;
	      H_p.SetPtEtaPhiE(genPart_mother_pt->at(i).at(j),genPart_mother_eta->at(i).at(j),genPart_mother_phi->at(i).at(j),genPart_mother_e->at(i).at(j));
	      cout<<"  mother mass is "<<H_p.M()<<endl;
	    }
	    break;
	  }
	}
	cout<<"GEN invariant mass: "<<(test_p+test_j).M()<<endl;
	cout<<endl;
	*/
	
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
	outFile1->cd();
	histpileup->Write();
    outFile1->Write();
    outFile1->Close();
	delete histpileup;
  }//end of sample loop
  cout<<"Number of events: "<<count0<<endl;
  cout<<"Number of events have photon and AK8 jet: "<<count1<<endl;
  cout<<"Number of events fired trigger: "<<count2<<endl;
  cout<<"Number of events pass photon pre-selection: "<<count3<<endl;
  cout<<"Number of events pass jet pre-selection: "<<count4<<endl;
  cout<<"Number of events pass photon mvaID: "<<count5<<endl;

  /*
  if(isPlot){
    //Non stacked plots
    hist01pa->SetLineWidth(2); hist01pa->SetLineColor(4);
    hist02pa->SetLineWidth(2); hist02pa->SetLineColor(4);
    hist03pa->SetLineWidth(2); hist03pa->SetLineColor(4);
    hist04pa->SetLineWidth(2); hist04pa->SetLineColor(4);
    hist05pa->SetLineWidth(2); hist05pa->SetLineColor(4);
    hist06pa->SetLineWidth(2); hist06pa->SetLineColor(4);
    hist07pa->SetLineWidth(2); hist07pa->SetLineColor(4);
    hist08pa->SetLineWidth(2); hist08pa->SetLineColor(4);
    hist09pa->SetLineWidth(2); hist09pa->SetLineColor(4);
    hist10pa->SetLineWidth(2); hist10pa->SetLineColor(4);
    hist11pa->SetLineWidth(2); hist11pa->SetLineColor(4);
    hist12pa->SetLineWidth(2); hist12pa->SetLineColor(4);
    hist13pa->SetLineWidth(2); hist13pa->SetLineColor(4);
    hist14pa->SetLineWidth(2); hist14pa->SetLineColor(4);
    hist15pa->SetLineWidth(2); hist15pa->SetLineColor(4);
    hist01pb->SetLineWidth(2); hist01pb->SetLineColor(2);
    hist02pb->SetLineWidth(2); hist02pb->SetLineColor(2);
    hist03pb->SetLineWidth(2); hist03pb->SetLineColor(2);
    hist04pb->SetLineWidth(2); hist04pb->SetLineColor(2);
    hist05pb->SetLineWidth(2); hist05pb->SetLineColor(2);
    hist06pb->SetLineWidth(2); hist06pb->SetLineColor(2);
    hist07pb->SetLineWidth(2); hist07pb->SetLineColor(2);
    hist08pb->SetLineWidth(2); hist08pb->SetLineColor(2);
    hist09pb->SetLineWidth(2); hist09pb->SetLineColor(2);
    hist10pb->SetLineWidth(2); hist10pb->SetLineColor(2);
    hist11pb->SetLineWidth(2); hist11pb->SetLineColor(2);
    hist12pb->SetLineWidth(2); hist12pb->SetLineColor(2);
    hist13pb->SetLineWidth(2); hist13pb->SetLineColor(2);
    hist14pb->SetLineWidth(2); hist14pb->SetLineColor(2);
    hist15pb->SetLineWidth(2); hist15pb->SetLineColor(2);
  
    hist01ja->SetLineWidth(2); hist01ja->SetLineColor(4);
    hist02ja->SetLineWidth(2); hist02ja->SetLineColor(4);
    hist03ja->SetLineWidth(2); hist03ja->SetLineColor(4);
    hist04ja->SetLineWidth(2); hist04ja->SetLineColor(4);
    hist05ja->SetLineWidth(2); hist05ja->SetLineColor(4);
    hist06ja->SetLineWidth(2); hist06ja->SetLineColor(4);
    hist07ja->SetLineWidth(2); hist07ja->SetLineColor(4);
    hist08ja->SetLineWidth(2); hist08ja->SetLineColor(4);
    //hist09ja->SetLineWidth(2); hist09ja->SetLineColor(4);
    hist10ja->SetLineWidth(2); hist10ja->SetLineColor(4);
    hist01jb->SetLineWidth(2); hist01jb->SetLineColor(2);
    hist02jb->SetLineWidth(2); hist02jb->SetLineColor(2);
    hist03jb->SetLineWidth(2); hist03jb->SetLineColor(2);
    hist04jb->SetLineWidth(2); hist04jb->SetLineColor(2);
    hist05jb->SetLineWidth(2); hist05jb->SetLineColor(2);
    hist06jb->SetLineWidth(2); hist06jb->SetLineColor(2);
    hist07jb->SetLineWidth(2); hist07jb->SetLineColor(2);
    hist08jb->SetLineWidth(2); hist08jb->SetLineColor(2);
    //hist09jb->SetLineWidth(2); hist09jb->SetLineColor(2);
    hist10jb->SetLineWidth(2); hist10jb->SetLineColor(2);

    hist01sa->SetLineWidth(2); hist01sa->SetLineColor(4);
    hist01sb->SetLineWidth(2); hist01sb->SetLineColor(2);
    
    TLegend *legend = new TLegend(0.6,0.75,0.85,0.85);
    TAxis *xaxis = NULL;
    TAxis *yaxis = NULL;
    TCanvas *c01p = new TCanvas("c01p","pt_{#gamma}",1200,900);
    xaxis = hist01pa->GetXaxis();
    yaxis = hist01pa->GetYaxis();
    xaxis->SetTitle("pt_{#gamma} (GeV)");
    yaxis->SetTitle("Entries / 50 GeV");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c01p->SetLogy();
    c01p->cd();
    hist01pa->Draw("HIST");
    hist01pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist01pa,histotitle+" Fired trigger","f");
    legend->AddEntry(hist01pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c01p->Print("p_pt.png");

    TCanvas *c02p = new TCanvas("c02p","#eta_{#gamma}",1200,900);
    xaxis = hist02pa->GetXaxis();
    yaxis = hist02pa->GetYaxis();
    xaxis->SetTitle("#eta_{#gamma}");
    yaxis->SetTitle("Entries / 0.2");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c02p->SetLogy();
    c02p->cd();
    hist02pa->Draw("HIST");
    hist02pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist02pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist02pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c02p->Print("p_eta.png");

    TCanvas *c03p = new TCanvas("c03p","#varphi_{#gamma}",1200,900);
    xaxis = hist03pa->GetXaxis();
    yaxis = hist03pa->GetYaxis();
    xaxis->SetTitle("#varphi_{#gamma} (rad)");
    yaxis->SetTitle("Entries / #pi/25");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c03p->SetLogy();
    c03p->cd();
    hist03pa->Draw("HIST");
    hist03pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist03pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist03pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c03p->Print("p_phi.png");

    TCanvas *c04p = new TCanvas("c04p","E_{#gamma}",1200,900);
    xaxis = hist04pa->GetXaxis();
    yaxis = hist04pa->GetYaxis();
    xaxis->SetTitle("E_{#gamma} (GeV)");
    yaxis->SetTitle("Entries / 50 GeV");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c04p->SetLogy();
    c04p->cd();
    hist04pa->Draw("HIST");
    hist04pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist04pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist04pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c04p->Print("p_e.png");

    TCanvas *c05p = new TCanvas("c05p","cos(#theta^*)_{#gamma}",1200,900);
    xaxis = hist05pa->GetXaxis();
    yaxis = hist05pa->GetYaxis();
    xaxis->SetTitle("cos(#theta^*)_{#gamma}");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0,1000);
    //c05p->SetLogy();
    c05p->cd();
    hist05pa->Draw("HIST");
    hist05pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist05pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist05pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c05p->Print("p_costhetastar.png");

    TCanvas *c06p = new TCanvas("c06p","hovere_{#gamma}",1200,900);
    xaxis = hist06pa->GetXaxis();
    yaxis = hist06pa->GetYaxis();
    xaxis->SetTitle("H/E_{#gamma}");
    yaxis->SetTitle("Entries / 0.02");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c06p->SetLogy();
    c06p->cd();
    hist06pa->Draw("HIST");
    hist06pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist06pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist06pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c06p->Print("p_hovere.png");

    TCanvas *c07p = new TCanvas("c07p","pt/M_{#gamma}",1200,900);
    xaxis = hist07pa->GetXaxis();
    yaxis = hist07pa->GetYaxis();
    xaxis->SetTitle("pt_{#gamma}/M");
    yaxis->SetTitle("Entries / 0.1");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c07p->SetLogy();
    c07p->cd();
    hist07pa->Draw("HIST");
    hist07pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist07pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist07pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c07p->Print("p_PtoverM.png");

    TCanvas *c08p = new TCanvas("c08p","Loose Photon ID",1200,900);
    xaxis = hist08pa->GetXaxis();
    yaxis = hist08pa->GetYaxis();
    xaxis->SetTitle("Loose Photon ID");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c08p->SetLogy();
    c08p->cd();
    hist08pa->Draw("HIST");
    hist08pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist08pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist08pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c08p->Print("p_looseID.png");

    TCanvas *c09p = new TCanvas("c09p","Medium Photon ID",1200,900);
    xaxis = hist09pa->GetXaxis();
    yaxis = hist09pa->GetYaxis();
    xaxis->SetTitle("Medium Photon ID");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c09p->SetLogy();
    c09p->cd();
    hist09pa->Draw("HIST");
    hist09pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist09pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist09pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c09p->Print("p_mediumID.png");

    TCanvas *c10p = new TCanvas("c10p","Tight Photon ID",1200,900);
    xaxis = hist10pa->GetXaxis();
    yaxis = hist10pa->GetYaxis();
    xaxis->SetTitle("Tight Photon ID");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c10p->SetLogy();
    c10p->cd();
    hist10pa->Draw("HIST");
    hist10pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist10pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist10pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c10p->Print("p_tightID.png");

    TCanvas *c11p = new TCanvas("c11p","Photon mvaID Value",1200,900);
    xaxis = hist11pa->GetXaxis();
    yaxis = hist11pa->GetYaxis();
    xaxis->SetTitle("Photon mvaID Value");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c11p->SetLogy();
    c11p->cd();
    hist11pa->Draw("HIST");
    hist11pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist11pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist11pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c11p->Print("p_mvaVal.png");

    TCanvas *c12p = new TCanvas("c12p","Photon mvaID Category",1200,900);
    xaxis = hist12pa->GetXaxis();
    yaxis = hist12pa->GetYaxis();
    xaxis->SetTitle("Photon mvaID Category");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c12p->SetLogy();
    c12p->cd();
    hist12pa->Draw("HIST");
    hist12pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist12pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist12pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c12p->Print("p_mvaCat.png");

    TCanvas *c13p = new TCanvas("c13p","isoGamma_{#gamma}",1200,900);
    xaxis = hist13pa->GetXaxis();
    yaxis = hist13pa->GetYaxis();
    xaxis->SetTitle("Photon mvaID Value");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c13p->SetLogy();
    c13p->cd();
    hist13pa->Draw("HIST");
    hist13pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist13pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist13pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c13p->Print("p_isoGamma.png");

    TCanvas *c14p = new TCanvas("c14p","isoCh_{#gamma}",1200,900);
    xaxis = hist14pa->GetXaxis();
    yaxis = hist14pa->GetYaxis();
    xaxis->SetTitle("Photon mvaID Category");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c14p->SetLogy();
    c14p->cd();
    hist14pa->Draw("HIST");
    hist14pb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist14pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist14pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c14p->Print("p_isoCh.png");

    TCanvas *c15p = new TCanvas("c15p","ElectronVeto_{#gamma}",1200,900);
    xaxis = hist15pa->GetXaxis();
    yaxis = hist15pa->GetYaxis();
    xaxis->SetTitle("ElectronVeto");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c15p->SetLogy();
    c15p->cd();
    hist15pa->Draw("HIST");
    legend->Clear();
    legend->AddEntry(hist15pa,histotitle+" Fired trigger" ,"f");
    legend->AddEntry(hist15pb,histotitle+" Pass Selection","f");
    legend->Draw();
    c15p->Print("p_EleVeto.png");


    TCanvas *c01j = new TCanvas("c01j","pt_{j}",1200,900);
    xaxis = hist01ja->GetXaxis();
    yaxis = hist01ja->GetYaxis();
    xaxis->SetTitle("pt_{j}");
    yaxis->SetTitle("Entries / 50 GeV");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c01j->SetLogy();
    c01j->cd();
    hist01ja->Draw("HIST");
    hist01jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist01ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist01jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c01j->Print("j_pt.png");

    TCanvas *c02j = new TCanvas("c02j","#eta_{j}",1200,900);
    xaxis = hist02ja->GetXaxis();
    yaxis = hist02ja->GetYaxis();
    xaxis->SetTitle("#eta_{j}");
    yaxis->SetTitle("Entries / 0.2");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c02j->SetLogy();
    c02j->cd();
    hist02ja->Draw("HIST");
    hist02jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist02ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist02jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c02j->Print("j_eta.png");

    TCanvas *c03j = new TCanvas("c03j","#varphi_{j}",1200,900);
    xaxis = hist03ja->GetXaxis();
    yaxis = hist03ja->GetYaxis();
    xaxis->SetTitle("#varphi_{j}");
    yaxis->SetTitle("Entries / #pi/25");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c03j->SetLogy();
    c03j->cd();
    hist03ja->Draw("HIST");
    hist03jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist03ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist03jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c03j->Print("j_phi.png");

    TCanvas *c04j = new TCanvas("c04j","E_{j}",1200,900);
    xaxis = hist04ja->GetXaxis();
    yaxis = hist04ja->GetYaxis();
    xaxis->SetTitle("E_{j}");
    yaxis->SetTitle("Entries / 60 GeV");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c04j->SetLogy();
    c04j->cd();
    hist04ja->Draw("HIST");
    hist04jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist04ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist04jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c04j->Print("j_e.png");

    TCanvas *c05j = new TCanvas("c05j","m_{j}",1200,900);
    xaxis = hist05ja->GetXaxis();
    yaxis = hist05ja->GetYaxis();
    xaxis->SetTitle("m_{j}");
    yaxis->SetTitle("Entries / 10 GeV");
    yaxis->SetTitleOffset(1.3);
    //yaxis->SetRangeUser(0.1,10000000);
    //c05j->SetLogy();
    c05j->cd();
    hist05ja->Draw("HIST");
    hist05jb->Draw("SAME");
    //g1->Draw("SAME");
    //g2->Draw("SAME");
    //f1->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist05ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist05jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c05j->Print("j_m.png");
  
    TCanvas *c06j = new TCanvas("c06j","#tau1_{j}",1200,900);
    xaxis = hist06ja->GetXaxis();
    yaxis = hist06ja->GetYaxis();
    xaxis->SetTitle("#tau1_{j}");
    yaxis->SetTitle("Entries / 0.02");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c06j->SetLogy();
    c06j->cd();
    hist06ja->Draw("HIST");
    hist06jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist06ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist06jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c06j->Print("j_tau1.png");

    TCanvas *c07j = new TCanvas("c07j","#tau2_{j}",1200,900);
    xaxis = hist07ja->GetXaxis();
    yaxis = hist07ja->GetYaxis();
    xaxis->SetTitle("#tau2_{j}");
    yaxis->SetTitle("Entries / 0.02");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c07j->SetLogy();
    c07j->cd();
    hist07ja->Draw("HIST");
    hist07jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist07ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist07jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c07j->Print("j_tau2.png");

    TCanvas *c08j = new TCanvas("c08j","#tau21_{j}",1200,900);
    xaxis = hist08ja->GetXaxis();
    yaxis = hist08ja->GetYaxis();
    xaxis->SetTitle("#tau21_{j}");
    yaxis->SetTitle("Entries / 0.02");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c08j->SetLogy();
    c08j->cd();
    hist08ja->Draw("HIST");
    hist08jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist08ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist08jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c08j->Print("j_tau21.png");
  */

    /*
      TCanvas *c09j = new TCanvas("c09j","TightLepVeto",1200,900);
      xaxis = hist09ja->GetXaxis();
      yaxis = hist09ja->GetYaxis();
      xaxis->SetTitle("Tight Lepton Veto");
      yaxis->SetTitle("Entries");
      yaxis->SetTitleOffset(1.3);
      yaxis->SetRangeUser(0.1,10000000);
      c09j->SetLogy();
      c09j->cd();
      hist09ja->Draw("HIST");
      hist09jb->Draw("SAME");
      legend->Clear();
      legend->AddEntry(hist09ja,histotitle+" Fired trigger","f");
      legend->AddEntry(hist09jb,histotitle+" Pass Selection","f");
      legend->Draw();
      c09j->Print("j_LepVeto.png");
    */
  /*
    TCanvas *c10j = new TCanvas("c10j","Tight jet ID",1200,900);
    xaxis = hist10ja->GetXaxis();
    yaxis = hist10ja->GetYaxis();
    xaxis->SetTitle("Tight jet ID");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c10j->SetLogy();
    c10j->cd();
    hist10ja->Draw("HIST");
    hist10jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist10ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist10jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c10j->Print("j_tightID.png");

    TCanvas *c11j = new TCanvas("c11j","Jet mass difference",1200,900);
    xaxis = hist11ja->GetXaxis();
    yaxis = hist11ja->GetYaxis();
    xaxis->SetTitle("Mass difference (GeV)");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    yaxis->SetRangeUser(0.1,10000000);
    c11j->SetLogy();
    c11j->cd();
    hist11ja->Draw("HIST");
    hist11jb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist11ja,histotitle+" Fired trigger","f");
    legend->AddEntry(hist11jb,histotitle+" Pass Selection","f");
    legend->Draw();
    c11j->Print("j_massdiff.png");

    TCanvas *c01s = new TCanvas("c01s","invmass",1200,900);
    xaxis = hist01sa->GetXaxis();
    yaxis = hist01sa->GetYaxis();
    xaxis->SetTitle("invariant mass (GeV)");
    yaxis->SetTitle("Entries / 40 GeV");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    //yaxis->SetRangeUser(0.1,10000000);
    //c01s->SetLogy();
    c01s->cd();
    hist01sa->Draw("HIST");
    hist01sb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist01sa,histotitle+" Fired trigger","f");
    legend->AddEntry(hist01sb,histotitle+" Pass Selection","f");
    legend->Draw();
    c01s->Print("invmass.png");

    TCanvas *c02s = new TCanvas("c02s","seperation",1200,900);
    xaxis = hist02sa->GetXaxis();
    yaxis = hist02sa->GetYaxis();
    xaxis->SetTitle("deltaR");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.3);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetRangeUser(0.1,10000000);
    c02s->SetLogy();
    c02s->cd();
    hist02sa->Draw("HIST");
    hist02sb->Draw("SAME");
    legend->Clear();
    legend->AddEntry(hist02sa,histotitle+" Fired trigger","f");
    legend->AddEntry(hist02sb,histotitle+" Pass Selection","f");
    legend->Draw();
    c02s->Print("deltaR.png");
    //Non stacked hist END-------------------------------------------------------------------

  }
*/
  
  gBenchmark->Show("selectWG");

}
