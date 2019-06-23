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
#include <TStopwatch.h>
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

void select(const TString conf="samples.conf", // input file
	      const TString outputDir=".",  // output directory
	      const Bool_t  doScaleCorr=0   // apply energy scale corrections?
	       ) {
  gBenchmark->Start("selectWG");
  gSystem->Load("lib/libMylib.so");
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //=============================================================================================================
  TString type = "DATA"; //DATA, SIGNAL, BKG, CONTROL
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

  //-----------------------------------CODE MODIFICATION INSTRUCTION FOR BKG/SIGNAL/DATA PROCESS---------------------------
  /*
  for QCD bkg, search for "SPECIFIC", uncomment lines with "STACK SPECIFIC" and "MC SPECIFIC"
  */

  //---------------------MC STACK HIST SPECIFIC-----------------------------
  /*
  //Initialize full hist
  TH1F *hs01 = new TH1F("hs01",histotitle+", pt_{#gamma}",48,0,2400);
  TH1F *hs02 = new TH1F("hs02",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F *hs03 = new TH1F("hs03",histotitle+", #varphi_{#gamma}",50,-3.14,3.14);
  TH1F *hs04 = new TH1F("hs04",histotitle+", E_{#gamma}",48,0,2400);
  TH1F *hs05 = new TH1F("hs05",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F *hs06 = new TH1F("hs06",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hs07 = new TH1F("hs07",histotitle+", cos(#theta^{*})_{#gamma}",50,0,1.1);
  TH1F *hs11 = new TH1F("hs11",histotitle+", pt_{j} (AK8)",48,0,2400);
  TH1F *hs12 = new TH1F("hs12",histotitle+", #eta_{j} (AK8)",50,-5,5);
  TH1F *hs13 = new TH1F("hs13",histotitle+", #varphi_{j} (AK8)",50,-3.14,3.14);
  TH1F *hs14 = new TH1F("hs14",histotitle+", E_{j} (AK8)",48,0,2400);
  TH1F *hs15 = new TH1F("hs15",histotitle+", m_{j} (AK8)",50,0,250);
  TH1F *hs16 = new TH1F("hs16",histotitle+", softdrop m_{j} (AK8)",50,0,250);
  TH1F *hs17 = new TH1F("hs17",histotitle+", HT",100,0,5000);
  
  
  //Initialize histograms for different MC samples
  //a: photon variable b: jet variable
  TH1F* histMC01a = new TH1F(); TH1F* histMC11a = new TH1F();//file1
  TH1F* histMC02a = new TH1F(); TH1F* histMC12a = new TH1F();
  TH1F* histMC03a = new TH1F(); TH1F* histMC13a = new TH1F();
  TH1F* histMC04a = new TH1F(); TH1F* histMC14a = new TH1F();
  TH1F* histMC05a = new TH1F(); TH1F* histMC15a = new TH1F();
  TH1F* histMC06a = new TH1F(); TH1F* histMC16a = new TH1F();
  TH1F* histMC07a = new TH1F(); TH1F* histMC17a = new TH1F();
  TH1F* histMC01b = new TH1F(); TH1F* histMC11b = new TH1F();//file2
  TH1F* histMC02b = new TH1F(); TH1F* histMC12b = new TH1F();
  TH1F* histMC03b = new TH1F(); TH1F* histMC13b = new TH1F();
  TH1F* histMC04b = new TH1F(); TH1F* histMC14b = new TH1F();
  TH1F* histMC05b = new TH1F(); TH1F* histMC15b = new TH1F();
  TH1F* histMC06b = new TH1F(); TH1F* histMC16b = new TH1F();
  TH1F* histMC07b = new TH1F(); TH1F* histMC17b = new TH1F();
  TH1F* histMC01c = new TH1F(); TH1F* histMC11c = new TH1F();//file3
  TH1F* histMC02c = new TH1F(); TH1F* histMC12c = new TH1F();
  TH1F* histMC03c = new TH1F(); TH1F* histMC13c = new TH1F();
  TH1F* histMC04c = new TH1F(); TH1F* histMC14c = new TH1F();
  TH1F* histMC05c = new TH1F(); TH1F* histMC15c = new TH1F();
  TH1F* histMC06c = new TH1F(); TH1F* histMC16c = new TH1F();
  TH1F* histMC07c = new TH1F(); TH1F* histMC17c = new TH1F();
  TH1F* histMC01d = new TH1F(); TH1F* histMC11d = new TH1F();//file4
  TH1F* histMC02d = new TH1F(); TH1F* histMC12d = new TH1F();
  TH1F* histMC03d = new TH1F(); TH1F* histMC13d = new TH1F();
  TH1F* histMC04d = new TH1F(); TH1F* histMC14d = new TH1F();
  TH1F* histMC05d = new TH1F(); TH1F* histMC15d = new TH1F();
  TH1F* histMC06d = new TH1F(); TH1F* histMC16d = new TH1F();
  TH1F* histMC07d = new TH1F(); TH1F* histMC17d = new TH1F();
  TH1F* histMC01e = new TH1F(); TH1F* histMC11e = new TH1F();//file5
  TH1F* histMC02e = new TH1F(); TH1F* histMC12e = new TH1F();
  TH1F* histMC03e = new TH1F(); TH1F* histMC13e = new TH1F();
  TH1F* histMC04e = new TH1F(); TH1F* histMC14e = new TH1F();
  TH1F* histMC05e = new TH1F(); TH1F* histMC15e = new TH1F();
  TH1F* histMC06e = new TH1F(); TH1F* histMC16e = new TH1F();
  TH1F* histMC07e = new TH1F(); TH1F* histMC17e = new TH1F();
  TH1F* histMC01f = new TH1F(); TH1F* histMC11f = new TH1F();//file6
  TH1F* histMC02f = new TH1F(); TH1F* histMC12f = new TH1F();
  TH1F* histMC03f = new TH1F(); TH1F* histMC13f = new TH1F();
  TH1F* histMC04f = new TH1F(); TH1F* histMC14f = new TH1F();
  TH1F* histMC05f = new TH1F(); TH1F* histMC15f = new TH1F();
  TH1F* histMC06f = new TH1F(); TH1F* histMC16f = new TH1F();
  TH1F* histMC07f = new TH1F(); TH1F* histMC17f = new TH1F();
  TH1F* histMC01g = new TH1F(); TH1F* histMC11g = new TH1F();//file7
  TH1F* histMC02g = new TH1F(); TH1F* histMC12g = new TH1F();
  TH1F* histMC03g = new TH1F(); TH1F* histMC13g = new TH1F();
  TH1F* histMC04g = new TH1F(); TH1F* histMC14g = new TH1F();
  TH1F* histMC05g = new TH1F(); TH1F* histMC15g = new TH1F();
  TH1F* histMC06g = new TH1F(); TH1F* histMC16g = new TH1F();
  TH1F* histMC07g = new TH1F(); TH1F* histMC17g = new TH1F();
  //--------------------MC STACK HIST SPECIFIC-----------------------------
  */

  // ------------------------Photons------------------------
  TH1F* hist01pa = new TH1F("WGamma01pa",histotitle+", pt_{#gamma}",48,0,2400); // Pass pre-selection and trigger selection
  TH1F* hist02pa = new TH1F("WGamma02pa",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pa = new TH1F("WGamma03pa",histotitle+", #varphi_{#gamma}",50,-3.14,3.14);
  TH1F* hist04pa = new TH1F("WGamma04pa",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist05pa = new TH1F("WGamma05pa",histotitle+", cos(#theta^*)_{#gamma}",50,0,1);
  TH1F* hist06pa = new TH1F("WGamma06pa",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist07pa = new TH1F("WGamma07pa",histotitle+", pt/M",50,0,1000);
  TH1F *hist08pa = new TH1F("WGamma08pa",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pa = new TH1F("WGamma09pa",histotitle+", Medium Photon ID",10,-1,3);
  TH1F *hist10pa = new TH1F("WGamma10pa",histotitle+", Tight Photon ID",10,-1,3);
  TH1F *hist11pa = new TH1F("WGamma11pa",histotitle+", Photon mvaID Value",50,-1.5,1.5);
  TH1F *hist12pa = new TH1F("WGamma12pa",histotitle+", Photon mvaID Category",10,-1,3);
  TH1F *hist13pa = new TH1F("WGamma13pa",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist14pa = new TH1F("WGamma14pa",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist15pa = new TH1F("WGamma15pa",histotitle+", Eletron Veto",10,-1,3);

  TH1F* hist01pb = new TH1F("WGamma01pb",histotitle+", pt_{#gamma}",48,0,2400); // Pass photon mvaID
  TH1F* hist02pb = new TH1F("WGamma02pb",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pb = new TH1F("WGamma03pb",histotitle+", #varphi_{#gamma}",50,-3.14,3.14);
  TH1F* hist04pb = new TH1F("WGamma04pb",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist05pb = new TH1F("WGamma05pb",histotitle+", cos(#theta^*)_{#gamma}",50,0,1);
  TH1F* hist06pb = new TH1F("WGamma06pb",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist07pb = new TH1F("WGamma07pb",histotitle+", pt/M",50,0,1000);
  TH1F *hist08pb = new TH1F("WGamma08pb",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pb = new TH1F("WGamma09pb",histotitle+", Medium Photon ID",10,-1,3);
  TH1F *hist10pb = new TH1F("WGamma10pb",histotitle+", Tight Photon ID",10,-1,3);
  TH1F *hist11pb = new TH1F("WGamma11pb",histotitle+", Photon mvaID Value",50,-1.5,1.5);
  TH1F *hist12pb = new TH1F("WGamma12pb",histotitle+", Photon mvaID Category",10,-1,3);
  TH1F *hist13pb = new TH1F("WGamma13pb",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist14pb = new TH1F("WGamma14pb",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist15pb = new TH1F("WGamma15pb",histotitle+", Eletron Veto",10,-1,3);
  //Jets
  TH1F* hist01ja = new TH1F("WGamma01ja",histotitle+", pt_{j} (AK8)",48,0,2400);
  TH1F* hist02ja = new TH1F("WGamma02ja",histotitle+", #eta_{j} (AK8)",50,-5,5);
  TH1F* hist03ja = new TH1F("WGamma03ja",histotitle+", #varphi_{j} (AK8)",50,-3.14,3.14);
  TH1F* hist04ja = new TH1F("WGamma04ja",histotitle+", E_{j} (AK8)",48,0,2400);
  TH1F* hist05ja = new TH1F("WGamma05ja",histotitle+", m_{j} (AK8)",50,0,250);
  TH1F* hist06ja = new TH1F("WGamma06ja",histotitle+", softdrop m_{j} (AK8)",50,0,250);
  TH1F* hist07ja = new TH1F("WGamma07ja",histotitle+", HT",100,0,5000);
  TH1F* hist08ja = new TH1F("WGamma08ja",histotitle+", Loose jet ID",50,-2,2);
  TH1F* hist09ja = new TH1F("WGamma09ja",histotitle+", Tight jet ID",50,-2,2);
  
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;
  gStyle->SetOptStat(0);

  
  // load trigger menu
  // const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  // load pileup reweighting file
  /*
  TFile *f_rw = TFile::Open("../Utils/data/puWeights_76x.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("puWeights");
  TH1D *h_rw_up = (TH1D*) f_rw->Get("puWeightsUp");
  TH1D *h_rw_down = (TH1D*) f_rw->Get("puWeightsDown");
  */


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
  //--Jets(AK8)--
  Int_t jetAK8_N = -99;
  std::vector<float> *jetAK8_pt = new std::vector<float>();
  std::vector<float> *jetAK8_eta = new std::vector<float>();
  std::vector<float> *jetAK8_phi = new std::vector<float>();
  std::vector<float> *jetAK8_mass = new std::vector<float>();
  std::vector<float> *jetAK8_E = new std::vector<float>();
  std::vector<float> *jetAK8_jec = new std::vector<float>();
  std::vector<float> *jetAK8_jecUp = new std::vector<float>();
  std::vector<float> *jetAK8_jecDown = new std::vector<float>();
  std::vector<bool>  *jetAK8_IDLoose = new std::vector<bool>();
  std::vector<bool>  *jetAK8_IDTight = new std::vector<bool>();
  std::vector<bool>  *jetAK8_IDTightLepVeto = new std::vector<bool>();
  std::vector<int>   *jetAK8_charge = new std::vector<int>();
  std::vector<int>   *jetAK8_partonFlavour = new std::vector<int>();
  std::vector<int>   *jetAK8_hadronFlavour = new std::vector<int>();
  std::vector<int>   *jetAK8_genParton_pdgID = new std::vector<int>();
  std::vector<int>   *jetAK8_nbHadrons = new std::vector<int>();
  std::vector<int>   *jetAK8_ncHadrons  = new std::vector<int>();
  std::vector<float> *jetAK8_jer_sf = new std::vector<float>();
  std::vector<float> *jetAK8_jer_sf_up = new std::vector<float>();
  std::vector<float> *jetAK8_jer_sf_down = new std::vector<float>();
  std::vector<float> *jetAK8_jer_sigma_pt = new std::vector<float>();
  //--Jets(AK8 pruned)
  /*
  std::vector<float> *jetAK8_chs_pruned_mass = new std::vector<float>();
  std::vector<float> *jetAK8_chs_pruned_massCorr = new std::vector<float>();
  std::vector<float> *jetAK8_chs_pruned_jec = new std::vector<float>();
  std::vector<float> *jetAK8_chs_pruned_jecUp = new std::vector<float>();
  std::vector<float> *jetAK8_chs_pruned_jecDown = new std::vector<float>();
  */
  //--Jets(AK8 softdrop)
  std::vector<float> *jetAK8_softdrop_mass = new std::vector<float>();
  std::vector<float> *jetAK8_softdrop_massCorr = new std::vector<float>();
  std::vector<float> *jetAK8_softdrop_jec = new std::vector<float>();
  std::vector<float> *jetAK8_softdrop_jecUp = new std::vector<float>();
  std::vector<float> *jetAK8_softdrop_jecDown = new std::vector<float>();
  //--Jets(AK8 puppi)
  /*
  UInt_t* jetAK8_puppi_N = new UInt_t();
  std::vector<float> *jetAK8_puppi_pt = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_eta = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_phi = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_mass = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_E = new std::vector<float>();
  //Jets(AK8 puppi pruned)(not avaliable in 2017 codes)
  //Jet(AK8 puppi softdrop)
  std::vector<float> *jetAK8_puppi_softdrop_mass = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_massCorr = new std::vector<float>();
  std::vector<float> *jetAK8_puppi_softdrop_jec = new std::vector<float>();
  */
  //Trigger
  std::map<std::string,bool> *HLT_isFired = new std::map<std::string,bool>();
  //std::vector<TriggerObjects>   *muon_hltMatchBits = new std::vector<std::bitset<256> >();
  
  // loop over samples
  TTree* eventTree = 0;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    CSample* samp = samplev[isam];
    cout<<"begin loop over files"<<endl;
    TStopwatch stopwatch;

    //Lumi-section check
    vector<RunLumiRangeMap*> Lumi_Photon175;
    Int_t nprescales = samp->prescaleJSONv.size();
    for(int iprescale=0; iprescale<nprescales; iprescale++){
      Lumi_Photon175.push_back(new RunLumiRangeMap());
      Lumi_Photon175.back()->addJSONFile(samp->prescaleJSONv[iprescale]);
    }

    // loop through files
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {  
      
      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile]; cout.flush();
      TFile *infile = TFile::Open(samp->fnamev[ifile]);
      assert(infile);

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
      eventTree->SetBranchAddress("ph_eta", &ph_eta);                           TBranch *photonEtaBr = eventTree->GetBranch("ph_eta");
      eventTree->SetBranchAddress("ph_phi", &ph_phi);                           TBranch *photonPhiBr = eventTree->GetBranch("ph_phi");
      eventTree->SetBranchAddress("ph_e", &ph_E);                               TBranch *photonEBr = eventTree->GetBranch("ph_e");
      eventTree->SetBranchAddress("ph_et", &ph_Et);                             TBranch *photonEtBr = eventTree->GetBranch("ph_et");
      eventTree->SetBranchAddress("ph_mass", &ph_m);                            TBranch *photonMBr = eventTree->GetBranch("ph_mass");
      eventTree->SetBranchAddress("ph_superCluster_eta", &ph_superCluster_eta); TBranch *photonSupCluEtaBr = eventTree->GetBranch("ph_superCluster_eta");
      eventTree->SetBranchAddress("ph_superCluster_phi", &ph_superCluster_phi); TBranch *photonSupCluPhiBr = eventTree->GetBranch("ph_superCluster_phi");
      eventTree->SetBranchAddress("ph_sigmaIetaIeta", &ph_sigmaIetaIeta);       TBranch *photonSigmaIEtaIEtaBr = eventTree->GetBranch("ph_sigmaIetaIeta");
      eventTree->SetBranchAddress("ph_hOverE", &ph_hOverE);                     TBranch *photonHoverEBr = eventTree->GetBranch("ph_hOverE");
      eventTree->SetBranchAddress("ph_isoGamma", &ph_isoGamma);                 TBranch *photonIsoGammaBr = eventTree->GetBranch("ph_isoGamma");
      eventTree->SetBranchAddress("ph_isoCh", &ph_isoCh);                       TBranch *photonIsoChBr = eventTree->GetBranch("ph_isoCh");
      eventTree->SetBranchAddress("ph_passEleVeto", &ph_passEleVeto);           TBranch *photonPassEleVetoBr = eventTree->GetBranch("ph_passEleVeto");
      eventTree->SetBranchAddress("ph_passLooseId", &ph_passLooseId);           TBranch *photonPassLooseIdBr = eventTree->GetBranch("ph_passLooseId");
      eventTree->SetBranchAddress("ph_passMediumId", &ph_passMediumId);         TBranch *photonPassMediumIdBr = eventTree->GetBranch("ph_passMediumId");
      eventTree->SetBranchAddress("ph_passTightId", &ph_passTightId);           TBranch *photonPassTightIdBr = eventTree->GetBranch("ph_passTightId");
      eventTree->SetBranchAddress("ph_mvaVal", &ph_mvaVal);                     TBranch *photonMvaValBr = eventTree->GetBranch("ph_mvaVal");
      eventTree->SetBranchAddress("ph_mvaCat", &ph_mvaCat);                     TBranch *photonMvaCatBr = eventTree->GetBranch("ph_mvaCat");
      //--Jets(AK8 softdrop)
      eventTree->SetBranchAddress("jetAK8_N", &jetAK8_N);                                 TBranch *jetAK8NBr = eventTree->GetBranch("jetAK8_N");
      eventTree->SetBranchAddress("jetAK8_pt", &jetAK8_pt);                               TBranch *jetAK8PtBr = eventTree->GetBranch("jetAK8_pt");
      eventTree->SetBranchAddress("jetAK8_eta", &jetAK8_eta);                             TBranch *jetAK8EtaBr = eventTree->GetBranch("jetAK8_eta");
      eventTree->SetBranchAddress("jetAK8_phi", &jetAK8_phi);                             TBranch *jetAK8PhiBr = eventTree->GetBranch("jetAK8_phi");
      eventTree->SetBranchAddress("jetAK8_e", &jetAK8_E);                                 TBranch *jetAK8EBr = eventTree->GetBranch("jetAK8_e");
      eventTree->SetBranchAddress("jetAK8_mass", &jetAK8_mass);                           TBranch *jetAK8MBr = eventTree->GetBranch("jetAK8_mass");
      eventTree->SetBranchAddress("jetAK8_jec", &jetAK8_jec);                             TBranch *jetAK8JecBr = eventTree->GetBranch("jetAK8_jec");
      eventTree->SetBranchAddress("jetAK8_jecUp", &jetAK8_jecUp);                         TBranch *jetAK8JecUpBr = eventTree->GetBranch("jetAK8_jecUp");
      eventTree->SetBranchAddress("jetAK8_jecDown", &jetAK8_jecDown);                     TBranch *jetAK8JecDownBr = eventTree->GetBranch("jetAK8_jecDown");
      eventTree->SetBranchAddress("jetAK8_IDLoose", &jetAK8_IDLoose);                     TBranch *jetAK8IdLooseBr = eventTree->GetBranch("jetAK8_IDLoose");
      eventTree->SetBranchAddress("jetAK8_IDTight", &jetAK8_IDTight);                     TBranch *jetAK8IdTightBr = eventTree->GetBranch("jetAK8_IDTight");
      eventTree->SetBranchAddress("jetAK8_IDTightLepVeto", &jetAK8_IDTightLepVeto);       TBranch *jetAK8IDTightLepVetoBr = eventTree->GetBranch("jetAK8_IDTightLepVeto");
      eventTree->SetBranchAddress("jetAK8_charge", &jetAK8_charge);                       TBranch *jetAK8ChargeBr = eventTree->GetBranch("jetAK8_charge");

      //MC specific variables------------------------------------------------
      /*
      eventTree->SetBranchAddress("jetAK8_partonFlavour", &jetAK8_partonFlavour);         TBranch *jetAK8PartonFlavourBr = eventTree->GetBranch("jetAK8_partonFlavour");
      eventTree->SetBranchAddress("jetAK8_hadronFlavour", &jetAK8_hadronFlavour);         TBranch *jetAK8HadronFlavourBr = eventTree->GetBranch("jetAK8_hadronFlavour");
      eventTree->SetBranchAddress("jetAK8_genParton_pdgID", &jetAK8_genParton_pdgID);     TBranch *jetAK8GenParton_pdgIDBr = eventTree->GetBranch("jetAK8_genParton_pdgID");
      eventTree->SetBranchAddress("jetAK8_nbHadrons", &jetAK8_nbHadrons);                 TBranch *jetAK8NbHadronsBr = eventTree->GetBranch("jetAK8_nbHadrons");
      eventTree->SetBranchAddress("jetAK8_ncHadrons", &jetAK8_ncHadrons);                 TBranch *jetAK8NcHadronsBr = eventTree->GetBranch("jetAK8_ncHadrons");
      eventTree->SetBranchAddress("jetAK8_jer_sf", &jetAK8_jer_sf);                       TBranch *jetAK8Jer_sfBr = eventTree->GetBranch("jetAK8_jer_sf");
      eventTree->SetBranchAddress("jetAK8_jer_sf_up", &jetAK8_jer_sf_up);                 TBranch *jetAK8Jer_sf_upBr = eventTree->GetBranch("jetAK8_jer_sf_up");
      eventTree->SetBranchAddress("jetAK8_jer_sf_down", &jetAK8_jer_sf_down);             TBranch *jetAK8Jer_sf_downBr = eventTree->GetBranch("jetAK8_jer_sf_down");
      eventTree->SetBranchAddress("jetAK8_jer_sigma_pt", &jetAK8_jer_sigma_pt);           TBranch *jetAK8Jer_sigma_ptBr = eventTree->GetBranch("jetAK8_jer_sigma_pt");
      */
      //----------------------------------------------------------------------
      eventTree->SetBranchAddress("jetAK8_softdrop_mass", &jetAK8_softdrop_mass);         TBranch *jetAK8SoftdropMBr = eventTree->GetBranch("jetAK8_softdrop_mass");
      eventTree->SetBranchAddress("jetAK8_softdrop_massCorr", &jetAK8_softdrop_massCorr); TBranch *jetAK8SoftdropMCorrBr = eventTree->GetBranch("jetAK8_softdrop_massCorr");
      eventTree->SetBranchAddress("jetAK8_softdrop_jec", &jetAK8_softdrop_jec);           TBranch *jetAK8SoftdropJecBr = eventTree->GetBranch("jetAK8_softdrop_jec");
      eventTree->SetBranchAddress("jetAK8_softdrop_jecUp", &jetAK8_softdrop_jecUp);       TBranch *jetAK8SoftdropJecUpBr = eventTree->GetBranch("jetAK8_softdrop_jecUp");
      eventTree->SetBranchAddress("jetAK8_softdrop_jecDown", &jetAK8_softdrop_jecDown);   TBranch *jetAK8SoftdropJecDownBr = eventTree->GetBranch("jetAK8_softdrop_jecDown");
      eventTree->SetBranchAddress("HLT_isFired", &HLT_isFired);                           TBranch *HLTisFiredBr = eventTree->GetBranch("HLT_isFired");

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {

	// Get Events
	runNumBr->GetEntry(ientry);
	evtNumBr->GetEntry(ientry);
	photonNBr->GetEntry(ientry);
	ph_pt->clear();                photonPtBr->GetEntry(ientry);
	ph_eta->clear();               photonEtaBr->GetEntry(ientry);
	ph_phi->clear();               photonPhiBr->GetEntry(ientry);
	ph_E->clear();                 photonEBr->GetEntry(ientry);
	ph_Et->clear();                photonEtBr->GetEntry(ientry);
	ph_m->clear();                 photonMBr->GetEntry(ientry);
	ph_superCluster_eta->clear();  photonSupCluEtaBr->GetEntry(ientry);
	ph_superCluster_phi->clear();  photonSupCluPhiBr->GetEntry(ientry);
	ph_sigmaIetaIeta->clear();     photonSigmaIEtaIEtaBr->GetEntry(ientry);
	ph_hOverE->clear();            photonHoverEBr->GetEntry(ientry);
	ph_isoGamma->clear();          photonIsoGammaBr->GetEntry(ientry);
	ph_isoCh->clear();             photonIsoChBr->GetEntry(ientry);
	ph_passEleVeto->clear();       photonPassEleVetoBr->GetEntry(ientry);
	ph_passLooseId->clear();       photonPassLooseIdBr->GetEntry(ientry);
	ph_passMediumId->clear();      photonPassMediumIdBr->GetEntry(ientry);
	ph_passTightId->clear();       photonPassTightIdBr->GetEntry(ientry);
	ph_mvaVal->clear();            photonMvaValBr->GetEntry(ientry);
	ph_mvaCat->clear();            photonMvaCatBr->GetEntry(ientry);
	jetAK8NBr->GetEntry(ientry);
	jetAK8_pt->clear();               jetAK8PtBr->GetEntry(ientry);
	jetAK8_eta->clear();              jetAK8EtaBr->GetEntry(ientry);
	jetAK8_phi->clear();              jetAK8PhiBr->GetEntry(ientry);
	jetAK8_mass->clear();             jetAK8MBr->GetEntry(ientry);
	jetAK8_E->clear();                jetAK8EBr->GetEntry(ientry);
	jetAK8_jec->clear();              jetAK8JecBr->GetEntry(ientry);
	jetAK8_jecUp->clear();            jetAK8JecUpBr->GetEntry(ientry);
	jetAK8_jecDown->clear();          jetAK8JecDownBr->GetEntry(ientry);
	jetAK8_IDLoose->clear();          jetAK8IdLooseBr->GetEntry(ientry);
	jetAK8_IDTight->clear();          jetAK8IdTightBr->GetEntry(ientry);
	jetAK8_IDTightLepVeto->clear();   jetAK8IDTightLepVetoBr->GetEntry(ientry);
	jetAK8_charge->clear();           jetAK8ChargeBr->GetEntry(ientry);
	//MC specific variables-----------------------------
	/*
	jetAK8_partonFlavour->clear();    jetAK8PartonFlavourBr->GetEntry(ientry);
	jetAK8_hadronFlavour->clear();    jetAK8HadronFlavourBr->GetEntry(ientry);
	jetAK8_genParton_pdgID->clear();  jetAK8GenParton_pdgIDBr->GetEntry(ientry);
	jetAK8_nbHadrons->clear();        jetAK8NbHadronsBr->GetEntry(ientry);
	jetAK8_ncHadrons->clear();        jetAK8NcHadronsBr->GetEntry(ientry);
	jetAK8_jer_sf->clear();           jetAK8Jer_sfBr->GetEntry(ientry);
	jetAK8_jer_sf_up->clear();        jetAK8Jer_sf_upBr->GetEntry(ientry);
	jetAK8_jer_sf_down->clear();      jetAK8Jer_sf_downBr->GetEntry(ientry);
	jetAK8_jer_sigma_pt->clear();     jetAK8Jer_sigma_ptBr->GetEntry(ientry);
	*/
	//--------------------------------------------------
	jetAK8_softdrop_mass->clear();    jetAK8SoftdropMBr->GetEntry(ientry);
	jetAK8_softdrop_massCorr->clear();jetAK8SoftdropMCorrBr->GetEntry(ientry);
	jetAK8_softdrop_jec->clear();     jetAK8SoftdropJecBr->GetEntry(ientry);
	jetAK8_softdrop_jecUp->clear();   jetAK8SoftdropJecUpBr->GetEntry(ientry);
	jetAK8_softdrop_jecDown->clear(); jetAK8SoftdropJecDownBr->GetEntry(ientry);
	HLTisFiredBr->GetEntry(ientry);

	//Only study events contain photons (EM objects here) and jets -- DATA -- 1st skimming
	if(ph_N < 1 || jetAK8_N <1) continue;

	// Locate the first mvaID photon in EM objects array (highest pt) and kinetic cut
	// See https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2
	Int_t index_p = -99;
	for(int i=0; i<ph_pt->size(); i++){
	  if(ph_passEleVeto->at(i) != true) continue;
	  if(abs(ph_eta->at(i)) > 2.4) continue;
	  if(abs(ph_eta->at(i)) > 1.4442 && abs(ph_eta->at(i))< 1.566) continue;
	  if((ph_mvaCat->at(i) == 0 && ph_mvaVal->at(i)< -0.02) || (ph_mvaCat->at(i) == 1 && ph_mvaVal->at(i)< -0.26)) continue;
	  index_p = i;
	  break;
	} 
	if(index_p == -99) continue;
	count1++;

	// Calculate HT (pt>10GeV, |eta|<3)
	Double_t HT = 0;
	for(int i=0; i<jetAK8_pt->size(); i++){
	  if(jetAK8_pt->at(i)>10 && abs(jetAK8_eta->at(i))<3)
	  HT += jetAK8_pt->at(i);
	}

	// Assigning leading photon and W jet
	TLorentzVector v_j, v_p, v_sys;
	v_j.SetPtEtaPhiE(jetAK8_pt->at(0), jetAK8_eta->at(0), jetAK8_phi->at(0), jetAK8_E->at(0));
	v_p.SetPtEtaPhiE(ph_pt->at(index_p), ph_eta->at(index_p), ph_phi->at(index_p), ph_E->at(index_p));
	v_sys = v_j + v_p;
	

	// Calculate pt/m
	Double_t invmass = 0;
	Double_t ptoverM = 0;
	invmass = (v_j+v_p).M();
	ptoverM = ph_pt->at(index_p)/invmass;

	// Calculate cos(theta*)
	Double_t cosThetaStar = -999;
	TLorentzVector v_boosted_j, v_boosted_p;
	v_boosted_p = v_p;
        v_boosted_p.Boost(-(v_sys.BoostVector()));
        cosThetaStar = abs(v_boosted_p.Pz()/v_boosted_p.P());

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

	// HLT Trigger decision
	bool passTrig = false;
	for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it) {
	  if (it->first.find("HLT_Photon200_v") != std::string::npos && it->second == 1){
	    passTrig = true;
	  }	  
	}
	if (!passTrig) continue;
	count2++;

	// Un-weighted Filling
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
	hist01ja->Fill(jetAK8_pt->at(0));
	hist02ja->Fill(jetAK8_eta->at(0));
	hist03ja->Fill(jetAK8_phi->at(0));
	hist04ja->Fill(jetAK8_E->at(0));
	hist05ja->Fill(jetAK8_mass->at(0));
	hist06ja->Fill(jetAK8_softdrop_mass->at(0));
	hist07ja->Fill(HT);
	hist08ja->Fill(jetAK8_IDLoose->at(0));
	hist09ja->Fill(jetAK8_IDTight->at(0));
      }//end of event loop
      cout<<"Number of events in this file: "<<eventTree->GetEntries()<<endl;
      cout<<"Events passed trigger: "<<count1<<endl;
      cout<<"Events with tight photon: "<<count2<<endl;
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
      //---------------------MC STACK HIST SPECIFIC-----------------------------
      /*
      Double_t scale = 1;//Scale to 1pb^-1
      if(ifile == 0){
	scale = double(1547000)/count1;
	cout<<"Scale is: "<<scale<<endl;
	histMC01a = (TH1F*)hist01->Clone(); histMC01a->SetLineColor(kBlue); histMC01a->SetLineWidth(2); histMC01a->SetLineStyle(2); histMC01a->Scale(scale);
	histMC02a = (TH1F*)hist02->Clone(); histMC02a->SetLineColor(kBlue); histMC02a->SetLineWidth(2);	histMC02a->SetLineStyle(2); histMC02a->Scale(scale);
	histMC03a = (TH1F*)hist03->Clone(); histMC03a->SetLineColor(kBlue); histMC03a->SetLineWidth(2);	histMC03a->SetLineStyle(2); histMC03a->Scale(scale);
	histMC04a = (TH1F*)hist04->Clone(); histMC04a->SetLineColor(kBlue); histMC04a->SetLineWidth(2);	histMC04a->SetLineStyle(2); histMC04a->Scale(scale);
	histMC05a = (TH1F*)hist05->Clone(); histMC05a->SetLineColor(kBlue); histMC05a->SetLineWidth(2);	histMC05a->SetLineStyle(2); histMC05a->Scale(scale);
	histMC06a = (TH1F*)hist06->Clone(); histMC06a->SetLineColor(kBlue); histMC06a->SetLineWidth(2);	histMC06a->
SetLineStyle(2); histMC06a->Scale(scale);
	histMC07a = (TH1F*)hist07->Clone(); histMC07a->SetLineColor(kBlue); histMC07a->SetLineWidth(2);	histMC07a->SetLineStyle(2); histMC07a->Scale(scale);
	histMC11a = (TH1F*)hist11->Clone(); histMC11a->SetLineColor(kBlue); histMC11a->SetLineWidth(2); histMC11a->SetLineStyle(2); histMC11a->Scale(scale);
	histMC12a = (TH1F*)hist12->Clone(); histMC12a->SetLineColor(kBlue); histMC12a->SetLineWidth(2); histMC12a->SetLineStyle(2); histMC12a->Scale(scale);
	histMC13a = (TH1F*)hist13->Clone(); histMC13a->SetLineColor(kBlue); histMC13a->SetLineWidth(2); histMC13a->SetLineStyle(2); histMC13a->Scale(scale);
	histMC14a = (TH1F*)hist14->Clone(); histMC14a->SetLineColor(kBlue); histMC14a->SetLineWidth(2); histMC14a->SetLineStyle(2); histMC14a->Scale(scale);
	histMC15a = (TH1F*)hist15->Clone(); histMC15a->SetLineColor(kBlue); histMC15a->SetLineWidth(2); histMC15a->SetLineStyle(2); histMC15a->Scale(scale);
	histMC16a = (TH1F*)hist16->Clone(); histMC16a->SetLineColor(kBlue); histMC16a->SetLineWidth(2); histMC16a->SetLineStyle(2); histMC16a->Scale(scale);
	histMC17a = (TH1F*)hist17->Clone(); histMC17a->SetLineColor(kBlue); histMC17a->SetLineWidth(2); histMC17a->SetLineStyle(2); histMC17a->Scale(scale);
      }
      else if(ifile == 1){
	scale = double(322600)/count1;
	cout<<"Scale is: "<<scale<<endl;
	histMC01b = (TH1F*)hist01->Clone(); histMC01b->SetLineColor(kRed); histMC01b->SetLineWidth(2); histMC01b->SetLineStyle(2); histMC01b->Scale(scale);
	histMC02b = (TH1F*)hist02->Clone(); histMC02b->SetLineColor(kRed); histMC02b->SetLineWidth(2); histMC02b->SetLineStyle(2); histMC02b->Scale(scale);
	histMC03b = (TH1F*)hist03->Clone(); histMC03b->SetLineColor(kRed); histMC03b->SetLineWidth(2); histMC03b->SetLineStyle(2); histMC03b->Scale(scale);
	histMC04b = (TH1F*)hist04->Clone(); histMC04b->SetLineColor(kRed); histMC04b->SetLineWidth(2); histMC04b->SetLineStyle(2); histMC04b->Scale(scale);
	histMC05b = (TH1F*)hist05->Clone(); histMC05b->SetLineColor(kRed); histMC05b->SetLineWidth(2); histMC05b->SetLineStyle(2); histMC05b->Scale(scale);
	histMC06b = (TH1F*)hist06->Clone(); histMC06b->SetLineColor(kRed); histMC06b->SetLineWidth(2); histMC06b->SetLineStyle(2); histMC06b->Scale(scale);
	histMC07b = (TH1F*)hist07->Clone(); histMC07b->SetLineColor(kRed); histMC07b->SetLineWidth(2); histMC07b->SetLineStyle(2); histMC07b->Scale(scale);
	histMC11b = (TH1F*)hist11->Clone(); histMC11b->SetLineColor(kRed); histMC11b->SetLineWidth(2); histMC11b->SetLineStyle(2); histMC11b->Scale(scale);
	histMC12b = (TH1F*)hist12->Clone(); histMC12b->SetLineColor(kRed); histMC12b->SetLineWidth(2); histMC12b->SetLineStyle(2); histMC12b->Scale(scale);
	histMC13b = (TH1F*)hist13->Clone(); histMC13b->SetLineColor(kRed); histMC13b->SetLineWidth(2); histMC13b->SetLineStyle(2); histMC13b->Scale(scale);
	histMC14b = (TH1F*)hist14->Clone(); histMC14b->SetLineColor(kRed); histMC14b->SetLineWidth(2); histMC14b->SetLineStyle(2); histMC14b->Scale(scale);
	histMC15b = (TH1F*)hist15->Clone(); histMC15b->SetLineColor(kRed); histMC15b->SetLineWidth(2); histMC15b->SetLineStyle(2); histMC15b->Scale(scale);
	histMC16b = (TH1F*)hist16->Clone(); histMC16b->SetLineColor(kRed); histMC16b->SetLineWidth(2); histMC16b->SetLineStyle(2); histMC16b->Scale(scale);
	histMC17b = (TH1F*)hist17->Clone(); histMC17b->SetLineColor(kRed); histMC17b->SetLineWidth(2); histMC17b->SetLineStyle(2); histMC17b->Scale(scale);
      }
      else if(ifile == 2){
	scale = double(29980)/count1;
	cout<<"Scale is: "<<scale<<endl;
	histMC01c = (TH1F*)hist01->Clone(); histMC01c->SetLineColor(kGreen); histMC01c->SetLineWidth(2); histMC01c->SetLineStyle(2); histMC01c->Scale(scale);
	histMC02c = (TH1F*)hist02->Clone(); histMC02c->SetLineColor(kGreen); histMC02c->SetLineWidth(2); histMC02c->SetLineStyle(2); histMC02c->Scale(scale);
	histMC03c = (TH1F*)hist03->Clone(); histMC03c->SetLineColor(kGreen); histMC03c->SetLineWidth(2); histMC03c->SetLineStyle(2); histMC03c->Scale(scale);
	histMC04c = (TH1F*)hist04->Clone(); histMC04c->SetLineColor(kGreen); histMC04c->SetLineWidth(2); histMC04c->SetLineStyle(2); histMC04c->Scale(scale);
	histMC05c = (TH1F*)hist05->Clone(); histMC05c->SetLineColor(kGreen); histMC05c->SetLineWidth(2); histMC05c->SetLineStyle(2); histMC05c->Scale(scale);
	histMC06c = (TH1F*)hist06->Clone(); histMC06c->SetLineColor(kGreen); histMC06c->SetLineWidth(2); histMC06c->SetLineStyle(2); histMC06c->Scale(scale);
	histMC07c = (TH1F*)hist07->Clone(); histMC07c->SetLineColor(kGreen); histMC07c->SetLineWidth(2); histMC07c->SetLineStyle(2); histMC07c->Scale(scale);
	histMC11c = (TH1F*)hist11->Clone(); histMC11c->SetLineColor(kGreen); histMC11c->SetLineWidth(2); histMC11c->SetLineStyle(2); histMC11c->Scale(scale);
	histMC12c = (TH1F*)hist12->Clone(); histMC12c->SetLineColor(kGreen); histMC12c->SetLineWidth(2); histMC12c->SetLineStyle(2); histMC12c->Scale(scale);
	histMC13c = (TH1F*)hist13->Clone(); histMC13c->SetLineColor(kGreen); histMC13c->SetLineWidth(2); histMC13c->SetLineStyle(2); histMC13c->Scale(scale);
	histMC14c = (TH1F*)hist14->Clone(); histMC14c->SetLineColor(kGreen); histMC14c->SetLineWidth(2); histMC14c->SetLineStyle(2); histMC14c->Scale(scale);
	histMC15c = (TH1F*)hist15->Clone(); histMC15c->SetLineColor(kGreen); histMC15c->SetLineWidth(2); histMC15c->SetLineStyle(2); histMC15c->Scale(scale);
	histMC16c = (TH1F*)hist16->Clone(); histMC16c->SetLineColor(kGreen); histMC16c->SetLineWidth(2); histMC16c->SetLineStyle(2); histMC16c->Scale(scale);
	histMC17c = (TH1F*)hist17->Clone(); histMC17c->SetLineColor(kGreen); histMC17c->SetLineWidth(2); histMC17c->SetLineStyle(2); histMC17c->Scale(scale);
      }
      else if(ifile == 3){
	scale = double(6334)/count1;
	cout<<"Scale is: "<<scale<<endl;
	histMC01d = (TH1F*)hist01->Clone(); histMC01d->SetLineColor(6); histMC01d->SetLineWidth(2); histMC01d->SetLineStyle(2); histMC01d->Scale(scale);
	histMC02d = (TH1F*)hist02->Clone(); histMC02d->SetLineColor(6); histMC02d->SetLineWidth(2); histMC02d->SetLineStyle(2); histMC02d->Scale(scale);
	histMC03d = (TH1F*)hist03->Clone(); histMC03d->SetLineColor(6); histMC03d->SetLineWidth(2); histMC03d->SetLineStyle(2); histMC03d->Scale(scale);
	histMC04d = (TH1F*)hist04->Clone(); histMC04d->SetLineColor(6); histMC04d->SetLineWidth(2); histMC04d->SetLineStyle(2); histMC04d->Scale(scale);
	histMC05d = (TH1F*)hist05->Clone(); histMC05d->SetLineColor(6); histMC05d->SetLineWidth(2); histMC05d->SetLineStyle(2); histMC05d->Scale(scale);
	histMC06d = (TH1F*)hist06->Clone(); histMC06d->SetLineColor(6); histMC06d->SetLineWidth(2); histMC06d->SetLineStyle(2); histMC06d->Scale(scale);
	histMC07d = (TH1F*)hist07->Clone(); histMC07d->SetLineColor(6); histMC07d->SetLineWidth(2); histMC07d->SetLineStyle(2); histMC07d->Scale(scale);
	histMC11d = (TH1F*)hist11->Clone(); histMC11d->SetLineColor(6); histMC11d->SetLineWidth(2); histMC11d->SetLineStyle(2); histMC11d->Scale(scale);
	histMC12d = (TH1F*)hist12->Clone(); histMC12d->SetLineColor(6); histMC12d->SetLineWidth(2); histMC12d->SetLineStyle(2); histMC12d->Scale(scale);
	histMC13d = (TH1F*)hist13->Clone(); histMC13d->SetLineColor(6); histMC13d->SetLineWidth(2); histMC13d->SetLineStyle(2); histMC13d->Scale(scale);
	histMC14d = (TH1F*)hist14->Clone(); histMC14d->SetLineColor(6); histMC14d->SetLineWidth(2); histMC14d->SetLineStyle(2); histMC14d->Scale(scale);
	histMC15d = (TH1F*)hist15->Clone(); histMC15d->SetLineColor(6); histMC15d->SetLineWidth(2); histMC15d->SetLineStyle(2); histMC15d->Scale(scale);
	histMC16d = (TH1F*)hist16->Clone(); histMC16d->SetLineColor(6); histMC16d->SetLineWidth(2); histMC16d->SetLineStyle(2); histMC16d->Scale(scale);
	histMC17d = (TH1F*)hist17->Clone(); histMC17d->SetLineColor(6); histMC17d->SetLineWidth(2); histMC17d->SetLineStyle(2); histMC17d->Scale(scale);
      }
      else if(ifile == 4){
	scale = double(1088)/count1;
	cout<<"Scale is: "<<scale<<endl;
	histMC01e = (TH1F*)hist01->Clone(); histMC01e->SetLineColor(kOrange-3); histMC01e->SetLineWidth(2); histMC01e->SetLineStyle(2); histMC01e->Scale(scale);
	histMC02e = (TH1F*)hist02->Clone(); histMC02e->SetLineColor(kOrange-3); histMC02e->SetLineWidth(2); histMC02e->SetLineStyle(2); histMC02e->Scale(scale);
	histMC03e = (TH1F*)hist03->Clone(); histMC03e->SetLineColor(kOrange-3); histMC03e->SetLineWidth(2); histMC03e->SetLineStyle(2); histMC03e->Scale(scale);
	histMC04e = (TH1F*)hist04->Clone(); histMC04e->SetLineColor(kOrange-3); histMC04e->SetLineWidth(2); histMC04e->SetLineStyle(2); histMC04e->Scale(scale);
	histMC05e = (TH1F*)hist05->Clone(); histMC05e->SetLineColor(kOrange-3); histMC05e->SetLineWidth(2); histMC05e->SetLineStyle(2); histMC05e->Scale(scale);
	histMC06e = (TH1F*)hist06->Clone(); histMC06e->SetLineColor(kOrange-3); histMC06e->SetLineWidth(2); histMC06e->SetLineStyle(2); histMC06e->Scale(scale);
	histMC07e = (TH1F*)hist07->Clone(); histMC07e->SetLineColor(kOrange-3); histMC07e->SetLineWidth(2); histMC07e->SetLineStyle(2); histMC07e->Scale(scale);
	histMC11e = (TH1F*)hist11->Clone(); histMC11e->SetLineColor(kOrange-3); histMC11e->SetLineWidth(2); histMC11e->SetLineStyle(2); histMC11e->Scale(scale);
	histMC12e = (TH1F*)hist12->Clone(); histMC12e->SetLineColor(kOrange-3); histMC12e->SetLineWidth(2); histMC12e->SetLineStyle(2); histMC12e->Scale(scale);
	histMC13e = (TH1F*)hist13->Clone(); histMC13e->SetLineColor(kOrange-3); histMC13e->SetLineWidth(2); histMC13e->SetLineStyle(2); histMC13e->Scale(scale);
	histMC14e = (TH1F*)hist14->Clone(); histMC14e->SetLineColor(kOrange-3); histMC14e->SetLineWidth(2); histMC14e->SetLineStyle(2); histMC14e->Scale(scale);
	histMC15e = (TH1F*)hist15->Clone(); histMC15e->SetLineColor(kOrange-3); histMC15e->SetLineWidth(2); histMC15e->SetLineStyle(2); histMC15e->Scale(scale);
	histMC16e = (TH1F*)hist16->Clone(); histMC16e->SetLineColor(kOrange-3); histMC16e->SetLineWidth(2); histMC16e->SetLineStyle(2); histMC16e->Scale(scale);
	histMC17e = (TH1F*)hist17->Clone(); histMC17e->SetLineColor(kOrange-3); histMC17e->SetLineWidth(2); histMC17e->SetLineStyle(2); histMC17e->Scale(scale);
      }
      else if(ifile == 5){
	scale = double(99.11)/count1;
	cout<<"Scale is: "<<scale<<endl;
	histMC01f = (TH1F*)hist01->Clone(); histMC01f->SetLineColor(7); histMC01f->SetLineWidth(2); histMC01f->SetLineStyle(2); histMC01f->Scale(scale);
	histMC02f = (TH1F*)hist02->Clone(); histMC02f->SetLineColor(7); histMC02f->SetLineWidth(2); histMC02f->SetLineStyle(2); histMC02f->Scale(scale);
	histMC03f = (TH1F*)hist03->Clone(); histMC03f->SetLineColor(7); histMC03f->SetLineWidth(2); histMC03f->SetLineStyle(2); histMC03f->Scale(scale);
	histMC04f = (TH1F*)hist04->Clone(); histMC04f->SetLineColor(7); histMC04f->SetLineWidth(2); histMC04f->SetLineStyle(2); histMC04f->Scale(scale);
	histMC05f = (TH1F*)hist05->Clone(); histMC05f->SetLineColor(7); histMC05f->SetLineWidth(2); histMC05f->SetLineStyle(2); histMC05f->Scale(scale);
	histMC06f = (TH1F*)hist06->Clone(); histMC06f->SetLineColor(7); histMC06f->SetLineWidth(2); histMC06f->SetLineStyle(2); histMC06f->Scale(scale);
	histMC07f = (TH1F*)hist07->Clone(); histMC07f->SetLineColor(7); histMC07f->SetLineWidth(2); histMC07f->SetLineStyle(2); histMC07f->Scale(scale);
	histMC11f = (TH1F*)hist11->Clone(); histMC11f->SetLineColor(7); histMC11f->SetLineWidth(2); histMC11f->SetLineStyle(2); histMC11f->Scale(scale);
	histMC12f = (TH1F*)hist12->Clone(); histMC12f->SetLineColor(7); histMC12f->SetLineWidth(2); histMC12f->SetLineStyle(2); histMC12f->Scale(scale);
	histMC13f = (TH1F*)hist13->Clone(); histMC13f->SetLineColor(7); histMC13f->SetLineWidth(2); histMC13f->SetLineStyle(2); histMC13f->Scale(scale);
	histMC14f = (TH1F*)hist14->Clone(); histMC14f->SetLineColor(7); histMC14f->SetLineWidth(2); histMC14f->SetLineStyle(2); histMC14f->Scale(scale);
	histMC15f = (TH1F*)hist15->Clone(); histMC15f->SetLineColor(7); histMC15f->SetLineWidth(2); histMC15f->SetLineStyle(2); histMC15f->Scale(scale);
	histMC16f = (TH1F*)hist16->Clone(); histMC16f->SetLineColor(7); histMC16f->SetLineWidth(2); histMC16f->SetLineStyle(2); histMC16f->Scale(scale);
	histMC17f = (TH1F*)hist17->Clone(); histMC17f->SetLineColor(7); histMC17f->SetLineWidth(2); histMC17f->SetLineStyle(2); histMC17f->Scale(scale);
      }
      else if(ifile == 6){
	scale = double(20.23)/count1;
	cout<<"Scale is: "<<scale<<endl;
	histMC01g = (TH1F*)hist01->Clone(); histMC01g->SetLineColor(8); histMC01g->SetLineWidth(2); histMC01g->SetLineStyle(2); histMC01g->Scale(scale);
	histMC02g = (TH1F*)hist02->Clone(); histMC02g->SetLineColor(8); histMC02g->SetLineWidth(2); histMC02g->SetLineStyle(2); histMC02g->Scale(scale);
	histMC03g = (TH1F*)hist03->Clone(); histMC03g->SetLineColor(8); histMC03g->SetLineWidth(2); histMC03g->SetLineStyle(2); histMC03g->Scale(scale);
	histMC04g = (TH1F*)hist04->Clone(); histMC04g->SetLineColor(8); histMC04g->SetLineWidth(2); histMC04g->SetLineStyle(2); histMC04g->Scale(scale);
	histMC05g = (TH1F*)hist05->Clone(); histMC05g->SetLineColor(8); histMC05g->SetLineWidth(2); histMC05g->SetLineStyle(2); histMC05g->Scale(scale);
	histMC06g = (TH1F*)hist06->Clone(); histMC06g->SetLineColor(8); histMC06g->SetLineWidth(2); histMC06g->SetLineStyle(2); histMC06g->Scale(scale);
	histMC07g = (TH1F*)hist07->Clone(); histMC07g->SetLineColor(8); histMC07g->SetLineWidth(2); histMC07g->SetLineStyle(2); histMC07g->Scale(scale);
	histMC11g = (TH1F*)hist11->Clone(); histMC11g->SetLineColor(8); histMC11g->SetLineWidth(2); histMC11g->SetLineStyle(2); histMC11g->Scale(scale);
	histMC12g = (TH1F*)hist12->Clone(); histMC12g->SetLineColor(8); histMC12g->SetLineWidth(2); histMC12g->SetLineStyle(2); histMC12g->Scale(scale);
	histMC13g = (TH1F*)hist13->Clone(); histMC13g->SetLineColor(8); histMC13g->SetLineWidth(2); histMC13g->SetLineStyle(2); histMC13g->Scale(scale);
	histMC14g = (TH1F*)hist14->Clone(); histMC14g->SetLineColor(8); histMC14g->SetLineWidth(2); histMC14g->SetLineStyle(2); histMC14g->Scale(scale);
	histMC15g = (TH1F*)hist15->Clone(); histMC15g->SetLineColor(8); histMC15g->SetLineWidth(2); histMC15g->SetLineStyle(2); histMC15g->Scale(scale);
	histMC16g = (TH1F*)hist16->Clone(); histMC16g->SetLineColor(8); histMC16g->SetLineWidth(2); histMC16g->SetLineStyle(2); histMC16g->Scale(scale);
	histMC17g = (TH1F*)hist17->Clone(); histMC17g->SetLineColor(8); histMC17g->SetLineWidth(2); histMC17g->SetLineStyle(2); histMC17g->Scale(scale);
      }  
      hist01->Reset();
      hist02->Reset();
      hist03->Reset();
      hist04->Reset();
      hist05->Reset();
      hist06->Reset();
      hist07->Reset();
      hist11->Reset();
      hist12->Reset();
      hist13->Reset();
      hist14->Reset();
      hist15->Reset();
      hist16->Reset();
      hist17->Reset();
      //-------------------------------MC STACK HIST SPECIFIC-------------------------------------------
      */
    }//end of file loop
  }//end of sample loop

  //Plots
  /*
  //---------------------MC STACK HIST SPECIFIC-----------------------------
  hs01->Add(histMC01a); hs01->Add(histMC01b); hs01->Add(histMC01c); hs01->Add(histMC01d); hs01->Add(histMC01e); hs01->Add(histMC01f); hs01->Add(histMC01g);
  hs02->Add(histMC02a); hs02->Add(histMC02b); hs02->Add(histMC02c); hs02->Add(histMC02d); hs02->Add(histMC02e); hs02->Add(histMC02f); hs02->Add(histMC02g);
  hs03->Add(histMC03a); hs03->Add(histMC03b); hs03->Add(histMC03c); hs03->Add(histMC03d); hs03->Add(histMC03e); hs03->Add(histMC03f); hs03->Add(histMC03g);
  hs04->Add(histMC04a); hs04->Add(histMC04b); hs04->Add(histMC04c); hs04->Add(histMC04d); hs04->Add(histMC04e); hs04->Add(histMC04f); hs04->Add(histMC04g);
  hs05->Add(histMC05a); hs05->Add(histMC05b); hs05->Add(histMC05c); hs05->Add(histMC05d); hs05->Add(histMC05e); hs05->Add(histMC05f); hs05->Add(histMC05g);
  hs06->Add(histMC06a); hs06->Add(histMC06b); hs06->Add(histMC06c); hs06->Add(histMC06d); hs06->Add(histMC06e); hs06->Add(histMC06f); hs06->Add(histMC06g);
  hs07->Add(histMC07a); hs07->Add(histMC07b); hs07->Add(histMC07c); hs07->Add(histMC07d); hs07->Add(histMC07e); hs07->Add(histMC07f); hs07->Add(histMC07g);
  hs11->Add(histMC11a); hs11->Add(histMC11b); hs11->Add(histMC11c); hs11->Add(histMC11d); hs11->Add(histMC11e); hs11->Add(histMC11f); hs11->Add(histMC11g);
  hs12->Add(histMC12a); hs12->Add(histMC12b); hs12->Add(histMC12c); hs12->Add(histMC12d); hs12->Add(histMC12e); hs12->Add(histMC12f); hs12->Add(histMC12g);
  hs13->Add(histMC13a); hs13->Add(histMC13b); hs13->Add(histMC13c); hs13->Add(histMC13d); hs13->Add(histMC13e); hs13->Add(histMC13f); hs13->Add(histMC13g);
  hs14->Add(histMC14a); hs14->Add(histMC14b); hs14->Add(histMC14c); hs14->Add(histMC14d); hs14->Add(histMC14e); hs14->Add(histMC14f); hs14->Add(histMC14g);
  hs15->Add(histMC15a); hs15->Add(histMC15b); hs15->Add(histMC15c); hs15->Add(histMC15d); hs15->Add(histMC15e); hs15->Add(histMC15f); hs15->Add(histMC15g);
  hs16->Add(histMC16a); hs16->Add(histMC16b); hs16->Add(histMC16c); hs16->Add(histMC16d); hs16->Add(histMC16e); hs16->Add(histMC16f); hs16->Add(histMC16g);
  hs17->Add(histMC17a); hs17->Add(histMC17b); hs17->Add(histMC17c); hs17->Add(histMC17d); hs17->Add(histMC17e); hs17->Add(histMC17f); hs17->Add(histMC17g);
  hist01->SetLineWidth(10);
  hist02->SetLineWidth(10);
  hist03->SetLineWidth(10);
  hist04->SetLineWidth(10);
  hist05->SetLineWidth(10);
  hist06->SetLineWidth(10);
  hist07->SetLineWidth(10);
  hist11->SetLineWidth(10);
  hist12->SetLineWidth(10);
  hist13->SetLineWidth(10);
  hist14->SetLineWidth(10);
  hist15->SetLineWidth(10);
  hist16->SetLineWidth(10);
  hist17->SetLineWidth(10);

  TLegend *legend = new TLegend(0.7,0.75,0.9,0.9);
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  TCanvas *c01 = new TCanvas("c01","pt_{#gamma}",1200,900);
  c01->cd();
  c01->SetLogy();
  hs01->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs01->GetXaxis();
  yaxis = hs01->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000); //NEED TO SPECIFY RANGE EXPLICITLY
  histMC01a->Draw("SAMEHIST");
  histMC01b->Draw("SAMEHIST");
  histMC01c->Draw("SAMEHIST");
  histMC01d->Draw("SAMEHIST");
  histMC01e->Draw("SAMEHIST");
  histMC01f->Draw("SAMEHIST");
  histMC01g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs01,"2017 QCD","f");
  legend->AddEntry(histMC01a,"HT 200-300","f");
  legend->AddEntry(histMC01b,"HT 300-500","f");
  legend->AddEntry(histMC01c,"HT 500-700","f");
  legend->AddEntry(histMC01d,"HT 700-1000","f");
  legend->AddEntry(histMC01e,"HT 1000-1500","f");
  legend->AddEntry(histMC01f,"HT 1500-2000","f");
  legend->AddEntry(histMC01g,"HT 2000-inf","f");
  legend->Draw();
  c01->Print("p_pt.png");

  TCanvas *c02 = new TCanvas("c02","#eta_{#gamma}",1200,900);
  c02->cd();
  c02->SetLogy();
  hs02->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs02->GetXaxis();
  yaxis = hs02->GetYaxis();
  xaxis->SetTitle("#eta_{#gamma}");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 0.2");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC02a->Draw("SAMEHIST");
  histMC02b->Draw("SAMEHIST");
  histMC02c->Draw("SAMEHIST");
  histMC02d->Draw("SAMEHIST");
  histMC02e->Draw("SAMEHIST");
  histMC02f->Draw("SAMEHIST");
  histMC02g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs02,"2017 QCD","f");
  legend->AddEntry(histMC02a,"HT 200-300","f");
  legend->AddEntry(histMC02b,"HT 300-500","f");
  legend->AddEntry(histMC02c,"HT 500-700","f");
  legend->AddEntry(histMC02d,"HT 700-1000","f");
  legend->AddEntry(histMC02e,"HT 1000-1500","f");
  legend->AddEntry(histMC02f,"HT 1500-2000","f");
  legend->AddEntry(histMC02g,"HT 2000-inf","f");
  legend->Draw();
  c02->Print("p_eta.png");

  TCanvas *c03 = new TCanvas("c03","#varphi_{#gamma}",1200,900);
  c03->cd();
  c03->SetLogy();
  hs03->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs03->GetXaxis();
  yaxis = hs03->GetYaxis();
  xaxis->SetTitle("#varphi_{#gamma} (rad)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / #pi/25");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC03a->Draw("SAMEHIST");
  histMC03b->Draw("SAMEHIST");
  histMC03c->Draw("SAMEHIST");
  histMC03d->Draw("SAMEHIST");
  histMC03e->Draw("SAMEHIST");
  histMC03f->Draw("SAMEHIST");
  histMC03g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs03,"2017 QCD","f");
  legend->AddEntry(histMC03a,"HT 200-300","f");
  legend->AddEntry(histMC03b,"HT 300-500","f");
  legend->AddEntry(histMC03c,"HT 500-700","f");
  legend->AddEntry(histMC03d,"HT 700-1000","f");
  legend->AddEntry(histMC03e,"HT 1000-1500","f");
  legend->AddEntry(histMC03f,"HT 1500-2000","f");
  legend->AddEntry(histMC03g,"HT 2000-inf","f");
  legend->Draw();
  c03->Print("p_phi.png");

  TCanvas *c04 = new TCanvas("c04","E_{#gamma}",1200,900);
  c04->cd();
  c04->SetLogy();
  hs04->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs04->GetXaxis();
  yaxis = hs04->GetYaxis();
  xaxis->SetTitle("E_{#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC04a->Draw("SAMEHIST");
  histMC04b->Draw("SAMEHIST");
  histMC04c->Draw("SAMEHIST");
  histMC04d->Draw("SAMEHIST");
  histMC04e->Draw("SAMEHIST");
  histMC04f->Draw("SAMEHIST");
  histMC04g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs04,"2017 QCD","f");
  legend->AddEntry(histMC04a,"HT 200-300","f");
  legend->AddEntry(histMC04b,"HT 300-500","f");
  legend->AddEntry(histMC04c,"HT 500-700","f");
  legend->AddEntry(histMC04d,"HT 700-1000","f");
  legend->AddEntry(histMC04e,"HT 1000-1500","f");
  legend->AddEntry(histMC04f,"HT 1500-2000","f");
  legend->AddEntry(histMC04g,"HT 2000-inf","f");
  legend->Draw();
  c04->Print("p_E.png");

  TCanvas *c05 = new TCanvas("c05","cos(#theta^*)_{#gamma}",1200,900);
  c05->cd();
  c05->SetLogy();
  hs05->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs05->GetXaxis();
  yaxis = hs05->GetYaxis();
  xaxis->SetTitle("cos(#theta^*)_{#gamma}");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC05a->Draw("SAMEHIST");
  histMC05b->Draw("SAMEHIST");
  histMC05c->Draw("SAMEHIST");
  histMC05d->Draw("SAMEHIST");
  histMC05e->Draw("SAMEHIST");
  histMC05f->Draw("SAMEHIST");
  histMC05g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs05,"2017 QCD","f");
  legend->AddEntry(histMC05a,"HT 200-300","f");
  legend->AddEntry(histMC05b,"HT 300-500","f");
  legend->AddEntry(histMC05c,"HT 500-700","f");
  legend->AddEntry(histMC05d,"HT 700-1000","f");
  legend->AddEntry(histMC05e,"HT 1000-1500","f");
  legend->AddEntry(histMC05f,"HT 1500-2000","f");
  legend->AddEntry(histMC05g,"HT 2000-inf","f");
  legend->Draw();
  c05->Print("p_costhetastar.png");

  TCanvas *c06 = new TCanvas("c06","H/E_{#gamma}",1200,900);
  c06->cd();
  c06->SetLogy();
  hs06->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs06->GetXaxis();
  yaxis = hs06->GetYaxis();
  xaxis->SetTitle("H/E_{#gamma}");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 0.02");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC06a->Draw("SAMEHIST");
  histMC06b->Draw("SAMEHIST");
  histMC06c->Draw("SAMEHIST");
  histMC06d->Draw("SAMEHIST");
  histMC06e->Draw("SAMEHIST");
  histMC06f->Draw("SAMEHIST");
  histMC06g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs06,"2017 QCD","f");
  legend->AddEntry(histMC06a,"HT 200-300","f");
  legend->AddEntry(histMC06b,"HT 300-500","f");
  legend->AddEntry(histMC06c,"HT 500-700","f");
  legend->AddEntry(histMC06d,"HT 700-1000","f");
  legend->AddEntry(histMC06e,"HT 1000-1500","f");
  legend->AddEntry(histMC06f,"HT 1500-2000","f");
  legend->AddEntry(histMC06g,"HT 2000-inf","f");
  legend->Draw();
  c06->Print("p_HoverE.png");

  TCanvas *c07 = new TCanvas("c07","cos(#theta^{*})_{#gamma}",1200,900);
  c07->cd();
  c07->SetLogy();
  hs07->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs07->GetXaxis();
  yaxis = hs07->GetYaxis();
  xaxis->SetTitle("cos(#theta^{*})_{#gamma}");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 0.04");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC07a->Draw("SAMEHIST");
  histMC07b->Draw("SAMEHIST");
  histMC07c->Draw("SAMEHIST");
  histMC07d->Draw("SAMEHIST");
  histMC07e->Draw("SAMEHIST");
  histMC07f->Draw("SAMEHIST");
  histMC07g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs07,"2017 QCD","f");
  legend->AddEntry(histMC07a,"HT 200-300","f");
  legend->AddEntry(histMC07b,"HT 300-500","f");
  legend->AddEntry(histMC07c,"HT 500-700","f");
  legend->AddEntry(histMC07d,"HT 700-1000","f");
  legend->AddEntry(histMC07e,"HT 1000-1500","f");
  legend->AddEntry(histMC07f,"HT 1500-2000","f");
  legend->AddEntry(histMC07g,"HT 2000-inf","f");
  legend->Draw();
  c07->Print("costheta.png");

  TCanvas *c11 = new TCanvas("c11","pt_{j}",1200,900);
  c11->cd();
  c11->SetLogy();
  hs11->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs11->GetXaxis();
  yaxis = hs11->GetYaxis();
  xaxis->SetTitle("pt_{j} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC11a->Draw("SAMEHIST");
  histMC11b->Draw("SAMEHIST");
  histMC11c->Draw("SAMEHIST");
  histMC11d->Draw("SAMEHIST");
  histMC11e->Draw("SAMEHIST");
  histMC11f->Draw("SAMEHIST");
  histMC11g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs11,"2017 QCD","f");
  legend->AddEntry(histMC11a,"HT 200-300","f");
  legend->AddEntry(histMC11b,"HT 300-500","f");
  legend->AddEntry(histMC11c,"HT 500-700","f");
  legend->AddEntry(histMC11d,"HT 700-1000","f");
  legend->AddEntry(histMC11e,"HT 1000-1500","f");
  legend->AddEntry(histMC11f,"HT 1500-2000","f");
  legend->AddEntry(histMC11g,"HT 2000-inf","f");
  legend->Draw();
  c11->Print("j_pt.png");

  TCanvas *c12 = new TCanvas("c12","#eta_{j}",1200,900);
  c12->cd();
  c12->SetLogy();
  hs12->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs12->GetXaxis();
  yaxis = hs12->GetYaxis();
  xaxis->SetTitle("#eta_{j}");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 0.2");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC12a->Draw("SAMEHIST");
  histMC12b->Draw("SAMEHIST");
  histMC12c->Draw("SAMEHIST");
  histMC12d->Draw("SAMEHIST");
  histMC12e->Draw("SAMEHIST");
  histMC12f->Draw("SAMEHIST");
  histMC12g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs12,"2017 QCD","f");
  legend->AddEntry(histMC12a,"HT 200-300","f");
  legend->AddEntry(histMC12b,"HT 300-500","f");
  legend->AddEntry(histMC12c,"HT 500-700","f");
  legend->AddEntry(histMC12d,"HT 700-1000","f");
  legend->AddEntry(histMC12e,"HT 1000-1500","f");
  legend->AddEntry(histMC12f,"HT 1500-2000","f");
  legend->AddEntry(histMC12g,"HT 2000-inf","f");
  legend->Draw();
  c12->Print("j_eta.png");

  TCanvas *c13 = new TCanvas("c13","#varphi_{j}",1200,900);
  c13->cd();
  c13->SetLogy();
  hs13->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs13->GetXaxis();
  yaxis = hs13->GetYaxis();
  xaxis->SetTitle("#varphi_{j} (rad)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / #pi/25");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC13a->Draw("SAMEHIST");
  histMC13b->Draw("SAMEHIST");
  histMC13c->Draw("SAMEHIST");
  histMC13d->Draw("SAMEHIST");
  histMC13e->Draw("SAMEHIST");
  histMC13f->Draw("SAMEHIST");
  histMC13g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs13,"2017 QCD","f");
  legend->AddEntry(histMC13a,"HT 200-300","f");
  legend->AddEntry(histMC13b,"HT 300-500","f");
  legend->AddEntry(histMC13c,"HT 500-700","f");
  legend->AddEntry(histMC13d,"HT 700-1000","f");
  legend->AddEntry(histMC13e,"HT 1000-1500","f");
  legend->AddEntry(histMC13f,"HT 1500-2000","f");
  legend->AddEntry(histMC13g,"HT 2000-inf","f");
  legend->Draw();
  c13->Print("j_phi.png");

  TCanvas *c14 = new TCanvas("c14","#E_{j}",1200,900);
  c14->cd();
  c14->SetLogy();
  hs14->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs14->GetXaxis();
  yaxis = hs14->GetYaxis();
  xaxis->SetTitle("E_{j} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC14a->Draw("SAMEHIST");
  histMC14b->Draw("SAMEHIST");
  histMC14c->Draw("SAMEHIST");
  histMC14d->Draw("SAMEHIST");
  histMC14e->Draw("SAMEHIST");
  histMC14f->Draw("SAMEHIST");
  histMC14g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs14,"2017 QCD","f");
  legend->AddEntry(histMC14a,"HT 200-300","f");
  legend->AddEntry(histMC14b,"HT 300-500","f");
  legend->AddEntry(histMC14c,"HT 500-700","f");
  legend->AddEntry(histMC14d,"HT 700-1000","f");
  legend->AddEntry(histMC14e,"HT 1000-1500","f");
  legend->AddEntry(histMC14f,"HT 1500-2000","f");
  legend->AddEntry(histMC14g,"HT 2000-inf","f");
  legend->Draw();
  c14->Print("j_E.png");

  TCanvas *c15 = new TCanvas("c15","#m_{j}",1200,900);
  c15->cd();
  c15->SetLogy();
  hs15->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs15->GetXaxis();
  yaxis = hs15->GetYaxis();
  xaxis->SetTitle("m_{j} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 5 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC15a->Draw("SAMEHIST");
  histMC15b->Draw("SAMEHIST");
  histMC15c->Draw("SAMEHIST");
  histMC15d->Draw("SAMEHIST");
  histMC15e->Draw("SAMEHIST");
  histMC15f->Draw("SAMEHIST");
  histMC15g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs15,"2017 QCD","f");
  legend->AddEntry(histMC15a,"HT 200-300","f");
  legend->AddEntry(histMC15b,"HT 300-500","f");
  legend->AddEntry(histMC15c,"HT 500-700","f");
  legend->AddEntry(histMC15d,"HT 700-1000","f");
  legend->AddEntry(histMC15e,"HT 1000-1500","f");
  legend->AddEntry(histMC15f,"HT 1500-2000","f");
  legend->AddEntry(histMC15g,"HT 2000-inf","f");
  legend->Draw();
  c15->Print("j_m.png");

  TCanvas *c16 = new TCanvas("c16","#m_softdrop_{j}",1200,900);
  c16->cd();
  c16->SetLogy();
  hs16->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs16->GetXaxis();
  yaxis = hs16->GetYaxis();
  xaxis->SetTitle("m_softdrop_{j} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 5 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC16a->Draw("SAMEHIST");
  histMC16b->Draw("SAMEHIST");
  histMC16c->Draw("SAMEHIST");
  histMC16d->Draw("SAMEHIST");
  histMC16e->Draw("SAMEHIST");
  histMC16f->Draw("SAMEHIST");
  histMC16g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs16,"2017 QCD","f");
  legend->AddEntry(histMC16a,"HT 200-300","f");
  legend->AddEntry(histMC16b,"HT 300-500","f");
  legend->AddEntry(histMC16c,"HT 500-700","f");
  legend->AddEntry(histMC16d,"HT 700-1000","f");
  legend->AddEntry(histMC16e,"HT 1000-1500","f");
  legend->AddEntry(histMC16f,"HT 1500-2000","f");
  legend->AddEntry(histMC16g,"HT 2000-inf","f");
  legend->Draw();
  c16->Print("j_m_softdrop.png");

  TCanvas *c17 = new TCanvas("c17","HT",1200,900);
  c17->cd();
  c17->SetLogy();
  hs17->Draw("HIST"); //YOU NEED TO DRAW FIRST TO GET AXIS
  xaxis = hs17->GetXaxis();
  yaxis = hs17->GetYaxis();
  xaxis->SetTitle("HT (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.01,10000000);
  histMC17a->Draw("SAMEHIST");
  histMC17b->Draw("SAMEHIST");
  histMC17c->Draw("SAMEHIST");
  histMC17d->Draw("SAMEHIST");
  histMC17e->Draw("SAMEHIST");
  histMC17f->Draw("SAMEHIST");
  histMC17g->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hs17,"2017 QCD","f");
  legend->AddEntry(histMC17a,"HT 200-300","f");
  legend->AddEntry(histMC17b,"HT 300-500","f");
  legend->AddEntry(histMC17c,"HT 500-700","f");
  legend->AddEntry(histMC17d,"HT 700-1000","f");
  legend->AddEntry(histMC17e,"HT 1000-1500","f");
  legend->AddEntry(histMC17f,"HT 1500-2000","f");
  legend->AddEntry(histMC17g,"HT 2000-inf","f");
  legend->Draw();
  c17->Print("HT.png");
  */

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
  hist01ja->SetLineWidth(2);
  hist02ja->SetLineWidth(2);
  hist03ja->SetLineWidth(2);
  hist04ja->SetLineWidth(2);
  hist05ja->SetLineWidth(2);
  hist06ja->SetLineWidth(2);
  hist07ja->SetLineWidth(2);
  hist08ja->SetLineWidth(2);
  hist09ja->SetLineWidth(2);
  
  TLegend *legend = new TLegend(0.7,0.85,0.9,0.9);
  
  TCanvas *c01p = new TCanvas("c01p","pt_{#gamma}",1200,900);
  TAxis *xaxis = hist01pa->GetXaxis();
  TAxis *yaxis = hist01pa->GetYaxis();
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
  legend->AddEntry(hist01pa,histotitle+" mvaID","f");
  legend->AddEntry(hist01pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist02pa,histotitle+" mvaID","f");
  legend->AddEntry(hist02pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist03pa,histotitle+" mvaID","f");
  legend->AddEntry(hist03pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist04pa,histotitle+" mvaID","f");
  legend->AddEntry(hist04pb,histotitle+" HLT_Photon200","f");
  legend->Draw();
  c04p->Print("p_e.png");

  TCanvas *c05p = new TCanvas("c05p","cos(#theta^*)_{#gamma}",1200,900);
  xaxis = hist05pa->GetXaxis();
  yaxis = hist05pa->GetYaxis();
  xaxis->SetTitle("cos(#theta^*)_{#gamma}");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.1,10000000);
  c05p->SetLogy();
  c05p->cd();
  hist05pa->Draw("HIST");
  hist05pb->Draw("SAME");
  legend->Clear();
  legend->AddEntry(hist05pa,histotitle+" mvaID","f");
  legend->AddEntry(hist05pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist06pa,histotitle+" mvaID","f");
  legend->AddEntry(hist06pb,histotitle+" HLT_Photon200","f");
  legend->Draw();
  c06p->Print("p_hovere.png");

  TCanvas *c07p = new TCanvas("c07p","pt/M_{#gamma}",1200,900);
  xaxis = hist07pa->GetXaxis();
  yaxis = hist07pa->GetYaxis();
  xaxis->SetTitle("pt/M_{#gamma}");
  yaxis->SetTitle("Entries / 20");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.1,10000000);
  c07p->SetLogy();
  c07p->cd();
  hist07pa->Draw("HIST");
  hist07pb->Draw("SAME");
  legend->Clear();
  legend->AddEntry(hist07pa,histotitle+" mvaID","f");
  legend->AddEntry(hist07pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist08pa,histotitle+" mvaID","f");
  legend->AddEntry(hist08pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist09pa,histotitle+" mvaID","f");
  legend->AddEntry(hist09pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist10pa,histotitle+" mvaID","f");
  legend->AddEntry(hist10pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist11pa,histotitle+" mvaID","f");
  legend->AddEntry(hist11pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist12pa,histotitle+" mvaID","f");
  legend->AddEntry(hist12pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist13pa,histotitle+" mvaID","f");
  legend->AddEntry(hist13pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist14pa,histotitle+" mvaID","f");
  legend->AddEntry(hist14pb,histotitle+" HLT_Photon200","f");
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
  legend->AddEntry(hist15pa,histotitle+" mvaID","f");
  legend->AddEntry(hist15pb,histotitle+" HLT_Photon200","f");
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
  legend->Clear();
  legend->AddEntry(hist01ja,histotitle,"f");
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
  legend->Clear();
  legend->AddEntry(hist02ja,histotitle,"f");
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
  legend->Clear();
  legend->AddEntry(hist03ja,histotitle,"f");
  legend->Draw();
  c03j->Print("j_phi.png");

  TCanvas *c04j = new TCanvas("c04j","E_{j}",1200,900);
  xaxis = hist04ja->GetXaxis();
  yaxis = hist04ja->GetYaxis();
  xaxis->SetTitle("E_{j}");
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,10000000);
  c04j->SetLogy();
  c04j->cd();
  hist04ja->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist04ja,histotitle,"f");
  legend->Draw();
  c04j->Print("j_e.png");

  TCanvas *c05j = new TCanvas("c05j","m_{j}",1200,900);
  xaxis = hist05ja->GetXaxis();
  yaxis = hist05ja->GetYaxis();
  xaxis->SetTitle("m_{j}");
  yaxis->SetTitle("Entries / 5 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,10000000);
  c05j->SetLogy();
  c05j->cd();
  hist05ja->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist05ja,histotitle,"f");
  legend->Draw();
  c05j->Print("j_m.png");

  TCanvas *c06j = new TCanvas("c06j","m_softdrop_{j}",1200,900);
  xaxis = hist06ja->GetXaxis();
  yaxis = hist06ja->GetYaxis();
  xaxis->SetTitle("m_softdrop_{j}");
  yaxis->SetTitle("Entries / 5 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,10000000);
  c06j->SetLogy();
  c06j->cd();
  hist06ja->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist06ja,histotitle,"f");
  legend->Draw();
  c06j->Print("j_m_softdrop.png");

  TCanvas *c07j = new TCanvas("c07j","HT",1200,900);
  xaxis = hist07ja->GetXaxis();
  yaxis = hist07ja->GetYaxis();
  xaxis->SetTitle("HT (GeV)");
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,10000000);
  c07j->SetLogy();
  c07j->cd();
  hist07ja->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist07ja,histotitle,"f");
  legend->Draw();
  c07j->Print("HT.png");

  TCanvas *c08j = new TCanvas("c08j","Loose jet ID",1200,900);
  xaxis = hist08ja->GetXaxis();
  yaxis = hist08ja->GetYaxis();
  xaxis->SetTitle("Loose jet ID");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,10000000);
  c08j->SetLogy();
  c08j->cd();
  hist08ja->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist08ja,histotitle,"f");
  legend->Draw();
  c08j->Print("j_looseID.png");

  TCanvas *c09j = new TCanvas("c09j","Tight jet ID",1200,900);
  xaxis = hist09ja->GetXaxis();
  yaxis = hist09ja->GetYaxis();
  xaxis->SetTitle("Tight jet ID");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,10000000);
  c09j->SetLogy();
  c09j->cd();
  hist09ja->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist09ja,histotitle,"f");
  legend->Draw();
  c09j->Print("j_tightID.png");
  //Non stacked hist END-------------------------------------------------------------------

  gBenchmark->Show("selectWG");

}
