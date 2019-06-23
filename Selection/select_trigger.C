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
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"
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
  TH1F* hist01pa = new TH1F("WGamma01pa",histotitle+", pt_{#gamma}",60,0,1200); //Pass: Loose ID and EleVeto
  TH1F* hist02pa = new TH1F("WGamma02pa",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pa = new TH1F("WGamma03pa",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pa = new TH1F("WGamma04pa",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pa = new TH1F("WGamma05pa",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pa = new TH1F("WGamma06pa",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pa = new TH1F("WGamma07pa",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pa = new TH1F("WGamma08pa",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pa = new TH1F("WGamma09pa",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pa = new TH1F("WGamma10pa",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pa = new TH1F("WGamma11pa",histotitle+", Photon mvaID Category",10,-1,3);

  TH1F* hist01pb = new TH1F("WGamma01pb",histotitle+", pt_{#gamma}",60,0,1200); //Fired trigger 175
  TH1F* hist02pb = new TH1F("WGamma02pb",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pb = new TH1F("WGamma03pb",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pb = new TH1F("WGamma04pb",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pb = new TH1F("WGamma05pb",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pb = new TH1F("WGamma06pb",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pb = new TH1F("WGamma07pb",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pb = new TH1F("WGamma08pb",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pb = new TH1F("WGamma09pb",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pb = new TH1F("WGamma10pb",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pb = new TH1F("WGamma11pb",histotitle+", Photon mvaID Category",10,-1,3);

  TH1F* hist01pc = new TH1F("WGamma01pc",histotitle+", pt_{#gamma}",60,0,1200); //Fired trigger 200
  TH1F* hist02pc = new TH1F("WGamma02pc",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pc = new TH1F("WGamma03pc",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pc = new TH1F("WGamma04pc",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pc = new TH1F("WGamma05pc",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pc = new TH1F("WGamma06pc",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pc = new TH1F("WGamma07pc",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pc = new TH1F("WGamma08pc",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pc = new TH1F("WGamma09pc",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pc = new TH1F("WGamma10pc",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pc = new TH1F("WGamma11pc",histotitle+", Photon mvaID Category",10,-1,3);

  TH1F* hist01pd = new TH1F("WGamma01pd",histotitle+", pt_{#gamma}",60,0,1200); //Fired trigger 150
  TH1F* hist02pd = new TH1F("WGamma02pd",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pd = new TH1F("WGamma03pd",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pd = new TH1F("WGamma04pd",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pd = new TH1F("WGamma05pd",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pd = new TH1F("WGamma06pd",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pd = new TH1F("WGamma07pd",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pd = new TH1F("WGamma08pd",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pd = new TH1F("WGamma09pd",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pd = new TH1F("WGamma10pd",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pd = new TH1F("WGamma11pd",histotitle+", Photon mvaID Category",10,-1,3);

  TH1F* hist01pe = new TH1F("WGamma01pe",histotitle+", pt_{#gamma}",60,0,1200); //Fired trigger 165
  TH1F* hist02pe = new TH1F("WGamma02pe",histotitle+", #eta_{#gamma}",50,-5,5);
  TH1F* hist03pe = new TH1F("WGamma03pe",histotitle+", E_{#gamma}",48,0,2400);
  TH1F* hist04pe = new TH1F("WGamma04pe",histotitle+", Et_{#gamma}",48,0,2400);
  TH1F* hist05pe = new TH1F("WGamma05pe",histotitle+", H/E_{#gamma}",50,0,1);
  TH1F *hist06pe = new TH1F("WGamma06pe",histotitle+", isoGamma_{#gamma}",50,0,20);
  TH1F *hist07pe = new TH1F("WGamma07pe",histotitle+", isoCh_{#gamma}",50,0,50);
  TH1F *hist08pe = new TH1F("WGamma08pe",histotitle+", Loose Photon ID",10,-1,3);
  TH1F *hist09pe = new TH1F("WGamma09pe",histotitle+", Pass EleVeto",10,-1,3);
  TH1F *hist10pe = new TH1F("WGamma10pe",histotitle+", Photon mvaID Value",10,-1,3);
  TH1F *hist11pe = new TH1F("WGamma11pe",histotitle+", Photon mvaID Category",10,-1,3);
  
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
  //--Trigger---
  std::map<std::string,bool> *HLT_isFired = new std::map<std::string,bool>();
  
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
    for(UInt_t ifile=0; ifile<10; ifile++) {  
      
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
      eventTree->SetBranchAddress("EVENT_lumiBlock", &lumiBlock);               TBranch *lumiBlockBr = eventTree->GetBranch("EVENT_lumiBlock");
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
      //--Trigger--
      eventTree->SetBranchAddress("HLT_isFired", &HLT_isFired);                 TBranch *HLTisFiredBr = eventTree->GetBranch("HLT_isFired");

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	// Get Events
	photonNBr->GetEntry(ientry);
	ph_pt->clear();                photonPtBr->GetEntry(ientry);
	ph_eta->clear();               photonEtaBr->GetEntry(ientry);
	ph_E->clear();                 photonEBr->GetEntry(ientry);
	ph_Et->clear();                photonEtBr->GetEntry(ientry);
	ph_hOverE->clear();            photonHoverEBr->GetEntry(ientry);
	ph_isoGamma->clear();          photonIsoGammaBr->GetEntry(ientry);
	ph_isoCh->clear();             photonIsoChBr->GetEntry(ientry);
	ph_passEleVeto->clear();       photonPassEleVetoBr->GetEntry(ientry);
	ph_passTightId->clear();       photonPassTightIdBr->GetEntry(ientry);
	ph_mvaVal->clear();            photonMvaValBr->GetEntry(ientry);
	ph_mvaCat->clear();            photonMvaCatBr->GetEntry(ientry);
	HLTisFiredBr->GetEntry(ientry);

	//Only study events contain photons (EM objects here) -- Muon Sample
	if(ph_N < 1) continue;

	//Locate the first photon in EM objects array (highest pt)
	Int_t index_p = -99;
	for(int i=0; i<ph_passTightId->size(); i++){
	  if(ph_passTightId->at(i) == 1 && ph_passEleVeto->at(i) == true){
	    //if(ph_passTightId->at(i) == 1 && ph_passEleVeto->at(i) == true && ph_pt->at(i) > 150){
	    index_p = i;
	    break;
	  }
	}
	if(index_p == -99) continue;
	count1++;
	hist01pa->Fill(ph_pt->at(index_p));
	hist02pa->Fill(ph_eta->at(index_p));
	hist03pa->Fill(ph_E->at(index_p));
	hist04pa->Fill(ph_Et->at(index_p));
	hist05pa->Fill(ph_hOverE->at(index_p));
	hist06pa->Fill(ph_isoGamma->at(index_p));
	hist07pa->Fill(ph_isoCh->at(index_p));
	hist08pa->Fill(ph_passTightId->at(index_p));
	hist09pa->Fill(ph_passEleVeto->at(index_p));
	hist10pa->Fill(ph_mvaVal->at(index_p));
	hist11pa->Fill(ph_mvaCat->at(index_p));

	//Trigger decision
	bool passTrig = false;
	for(map<string,bool>::iterator it = HLT_isFired->begin(); it != HLT_isFired->end(); ++it) {
	  if (it->first.find("HLT_Photon175") != std::string::npos && it->second == 1){
	    hist01pb->Fill(ph_pt->at(index_p));
	    hist02pb->Fill(ph_eta->at(index_p));
	    hist03pb->Fill(ph_E->at(index_p));
	    hist04pb->Fill(ph_Et->at(index_p));
	    hist05pb->Fill(ph_hOverE->at(index_p));
	    hist06pb->Fill(ph_isoGamma->at(index_p));
	    hist07pb->Fill(ph_isoCh->at(index_p));
	    hist08pb->Fill(ph_passTightId->at(index_p));
	    hist09pb->Fill(ph_passEleVeto->at(index_p));
	    hist10pb->Fill(ph_mvaVal->at(index_p));
	    hist11pb->Fill(ph_mvaCat->at(index_p));
	    passTrig = true;
	  }
	  /*
	  if (it->first.find("HLT_Photon200") != std::string::npos && it->second == 1){
	    hist01pc->Fill(ph_pt->at(index_p));
	    hist02pc->Fill(ph_eta->at(index_p));
	    hist03pc->Fill(ph_E->at(index_p));
	    hist04pc->Fill(ph_Et->at(index_p));
	    hist05pc->Fill(ph_hOverE->at(index_p));
	    hist06pc->Fill(ph_isoGamma->at(index_p));
	    hist07pc->Fill(ph_isoCh->at(index_p));
	    hist08pc->Fill(ph_passTightId->at(index_p));
	    hist09pc->Fill(ph_passEleVeto->at(index_p));
	    hist10pc->Fill(ph_mvaVal->at(index_p));
	    hist11pc->Fill(ph_mvaCat->at(index_p));
	    passTrig = true;
	  }
	  */
	  /*
	  if (it->first.find("HLT_Photon150") != std::string::npos && it->second == 1){
	    hist01pd->Fill(ph_pt->at(index_p));
	    hist02pd->Fill(ph_eta->at(index_p));
	    hist03pd->Fill(ph_E->at(index_p));
	    hist04pd->Fill(ph_Et->at(index_p));
	    hist05pd->Fill(ph_hOverE->at(index_p));
	    hist06pd->Fill(ph_isoGamma->at(index_p));
	    hist07pd->Fill(ph_isoCh->at(index_p));
	    hist08pd->Fill(ph_passTightId->at(index_p));
	    hist09pd->Fill(ph_passEleVeto->at(index_p));
	    hist10pd->Fill(ph_mvaVal->at(index_p));
	    hist11pd->Fill(ph_mvaCat->at(index_p));
	    passTrig = true;
	  }
	  if (it->first.find("HLT_Photon165") != std::string::npos && it->second == 1){
	    hist01pe->Fill(ph_pt->at(index_p));
	    hist02pe->Fill(ph_eta->at(index_p));
	    hist03pe->Fill(ph_E->at(index_p));
	    hist04pe->Fill(ph_Et->at(index_p));
	    hist05pe->Fill(ph_hOverE->at(index_p));
	    hist06pe->Fill(ph_isoGamma->at(index_p));
	    hist07pe->Fill(ph_isoCh->at(index_p));
	    hist08pe->Fill(ph_passTightId->at(index_p));
	    hist09pe->Fill(ph_passEleVeto->at(index_p));
	    hist10pe->Fill(ph_mvaVal->at(index_p));
	    hist11pe->Fill(ph_mvaCat->at(index_p));
	    passTrig = true;
	  }
	  */
	}
	if (!passTrig) continue;
	count2++;

	//Trigger pre-scale
	Int_t tri_prescale = -99;
	RunLumiRangeMap::RunLumiPairType rl(runnum, lumiBlock); 
	for(int iprescale=0; iprescale<nprescales; iprescale++){
	  if(Lumi_Photon175[iprescale]->hasRunLumi(rl)){
	    tri_prescale = samp->prescalev[iprescale];
	    cout<<tri_prescale<<endl;
	  }
	}
      }//end of event loop
      cout<<"Number of events in this file: "<<eventTree->GetEntries()<<endl;
      cout<<"Events with more than 1 photons: "<<count1<<endl;
      cout<<"Events passed trigger: "<<count2<<endl;
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
  }//end of sample loop

  //Tirgger Efficiency
  TEfficiency *Eff1 = new TEfficiency(*hist01pb, *hist01pa);
  TEfficiency *Eff2 = new TEfficiency(*hist01pc, *hist01pa);
  TEfficiency *Eff3 = new TEfficiency(*hist01pd, *hist01pa);
  TEfficiency *Eff4 = new TEfficiency(*hist01pe, *hist01pa);

  //Non stacked plots
  TLegend *legend1 = new TLegend(0.6,0.75,0.85,0.85);
  TLegend *legend2 = new TLegend(0.6,0.7,0.85,0.8);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  TCanvas *c01 = new TCanvas("c01","pt_{#gamma}",1200,900);
  xaxis = hist01pa->GetXaxis();
  yaxis = hist01pa->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 20 GeV");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c01->SetLogy();
  c01->cd();
  hist01pa->SetLineWidth(2);
  hist01pa->SetLineColor(1);
  hist01pa->Draw("HIST");
  hist01pb->SetLineWidth(2);
  hist01pb->SetLineColor(2);
  hist01pb->Draw("SAMEHIST");
  hist01pc->SetLineWidth(2);
  hist01pc->SetLineColor(4);
  hist01pc->Draw("SAMEHIST");
  hist01pd->SetLineWidth(2);
  hist01pd->SetLineColor(6);
  hist01pd->Draw("SAMEHIST");
  hist01pe->SetLineWidth(2);
  hist01pe->SetLineColor(8);
  hist01pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist01pa, histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist01pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist01pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist01pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist01pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c01->Print("p_pt.png");

  TCanvas *c02 = new TCanvas("c02","#eta_{#gamma}",1200,900);
  xaxis = hist02pa->GetXaxis();
  yaxis = hist02pa->GetYaxis();
  xaxis->SetTitle("#eta_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 0.2");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c02->SetLogy();
  c02->cd();
  hist02pa->SetLineWidth(2);
  hist02pa->SetLineColor(1);
  hist02pa->Draw("HIST");
  hist02pb->SetLineWidth(2);
  hist02pb->SetLineColor(2);
  hist02pb->Draw("SAMEHIST");
  hist02pc->SetLineWidth(2);
  hist02pc->SetLineColor(4);
  hist02pc->Draw("SAMEHIST");
  hist02pd->SetLineWidth(2);
  hist02pd->SetLineColor(6);
  hist02pd->Draw("SAMEHIST");
  hist02pe->SetLineWidth(2);
  hist02pe->SetLineColor(8);
  hist02pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist02pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist02pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist02pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist02pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist02pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c02->Print("p_eta.png");

  TCanvas *c03 = new TCanvas("c03","E_{#gamma}",1200,900);
  xaxis = hist03pa->GetXaxis();
  yaxis = hist03pa->GetYaxis();
  xaxis->SetTitle("E_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c03->SetLogy();
  c03->cd();
  hist03pa->SetLineWidth(2);
  hist03pa->SetLineColor(1);
  hist03pa->Draw("HIST");
  hist03pb->SetLineWidth(2);
  hist03pb->SetLineColor(2);
  hist03pb->Draw("SAMEHIST");
  hist03pc->SetLineWidth(2);
  hist03pc->SetLineColor(4);
  hist03pc->Draw("SAMEHIST");
  hist03pd->SetLineWidth(2);
  hist03pd->SetLineColor(6);
  hist03pd->Draw("SAMEHIST");
  hist03pe->SetLineWidth(2);
  hist03pe->SetLineColor(8);
  hist03pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist03pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist03pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist03pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist03pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist03pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c03->Print("p_e.png");

  TCanvas *c04 = new TCanvas("c04","Et_{#gamma}",1200,900);
  xaxis = hist04pa->GetXaxis();
  yaxis = hist04pa->GetYaxis();
  xaxis->SetTitle("Et_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c04->SetLogy();
  c04->cd();
  hist04pa->SetLineWidth(2);
  hist04pa->SetLineColor(1);
  hist04pa->Draw("HIST");
  hist04pb->SetLineWidth(2);
  hist04pb->SetLineColor(2);
  hist04pb->Draw("SAMEHIST");
  hist04pc->SetLineWidth(2);
  hist04pc->SetLineColor(4);
  hist04pc->Draw("SAMEHIST");
  hist04pd->SetLineWidth(2);
  hist04pd->SetLineColor(6);
  hist04pd->Draw("SAMEHIST");
  hist04pe->SetLineWidth(2);
  hist04pe->SetLineColor(8);
  hist04pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist04pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist04pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist04pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist04pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist04pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c04->Print("p_et.png");

  TCanvas *c05 = new TCanvas("c05","H/E_{#gamma}",1200,900);
  xaxis = hist05pa->GetXaxis();
  yaxis = hist05pa->GetYaxis();
  xaxis->SetTitle("H/E_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 0.02");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c05->SetLogy();
  c05->cd();
  hist05pa->SetLineWidth(2);
  hist05pa->SetLineColor(1);
  hist05pa->Draw("HIST");
  hist05pb->SetLineWidth(2);
  hist05pb->SetLineColor(2);
  hist05pb->Draw("SAMEHIST");
  hist05pc->SetLineWidth(2);
  hist05pc->SetLineColor(4);
  hist05pc->Draw("SAMEHIST");
  hist05pd->SetLineWidth(2);
  hist05pd->SetLineColor(6);
  hist05pd->Draw("SAMEHIST");
  hist05pe->SetLineWidth(2);
  hist05pe->SetLineColor(8);
  hist05pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist05pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist05pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist05pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist05pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist05pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c05->Print("p_hOverE.png");

  TCanvas *c06 = new TCanvas("c06","isoGamma_{#gamma}",1200,900);
  xaxis = hist06pa->GetXaxis();
  yaxis = hist06pa->GetYaxis();
  xaxis->SetTitle("isoGamma_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 0.4");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c06->SetLogy();
  c06->cd();
  hist06pa->SetLineWidth(2);
  hist06pa->SetLineColor(1);
  hist06pa->Draw("HIST");
  hist06pb->SetLineWidth(2);
  hist06pb->SetLineColor(2);
  hist06pb->Draw("SAMEHIST");
  hist06pc->SetLineWidth(2);
  hist06pc->SetLineColor(4);
  hist06pc->Draw("SAMEHIST");
  hist06pd->SetLineWidth(2);
  hist06pd->SetLineColor(6);
  hist06pd->Draw("SAMEHIST");
  hist06pe->SetLineWidth(2);
  hist06pe->SetLineColor(8);
  hist06pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist06pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist06pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist06pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist06pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist06pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c06->Print("p_isoGamma.png");

  TCanvas *c07 = new TCanvas("c07","isoCh_{#gamma}",1200,900);
  xaxis = hist07pa->GetXaxis();
  yaxis = hist07pa->GetYaxis();
  xaxis->SetTitle("isoCh_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 1");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c07->SetLogy();
  c07->cd();
  hist07pa->SetLineWidth(2);
  hist07pa->SetLineColor(1);
  hist07pa->Draw("HIST");
  hist07pb->SetLineWidth(2);
  hist07pb->SetLineColor(2);
  hist07pb->Draw("SAMEHIST");
  hist07pc->SetLineWidth(2);
  hist07pc->SetLineColor(4);
  hist07pc->Draw("SAMEHIST");
  hist07pd->SetLineWidth(2);
  hist07pd->SetLineColor(6);
  hist07pd->Draw("SAMEHIST");
  hist07pe->SetLineWidth(2);
  hist07pe->SetLineColor(8);
  hist07pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist07pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist07pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist07pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist07pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist07pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c07->Print("p_isoCh.png");

  TCanvas *c08 = new TCanvas("c08","Tight Photon ID",1200,900);
  xaxis = hist08pa->GetXaxis();
  yaxis = hist08pa->GetYaxis();
  xaxis->SetTitle("Tight Photon ID");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,1000000);
  c08->SetLogy();
  c08->cd();
  hist08pa->SetLineWidth(2);
  hist08pa->SetLineColor(1);
  hist08pa->Draw("HIST");
  hist08pb->SetLineWidth(2);
  hist08pb->SetLineColor(2);
  hist08pb->Draw("SAMEHIST");
  hist08pc->SetLineWidth(2);
  hist08pc->SetLineColor(4);
  hist08pc->Draw("SAMEHIST");
  hist08pd->SetLineWidth(2);
  hist08pd->SetLineColor(6);
  hist08pd->Draw("SAMEHIST");
  hist08pe->SetLineWidth(2);
  hist08pe->SetLineColor(8);
  hist08pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist08pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist08pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist08pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist08pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist08pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c08->Print("p_TightID.png");

  TCanvas *c09 = new TCanvas("c09","Electron Veto",1200,900);
  xaxis = hist09pa->GetXaxis();
  yaxis = hist09pa->GetYaxis();
  xaxis->SetTitle("Electron Veto");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,1000000);
  c09->SetLogy();
  c09->cd();
  hist09pa->SetLineWidth(2);
  hist09pa->SetLineColor(1);
  hist09pa->Draw("HIST");
  hist09pb->SetLineWidth(2);
  hist09pb->SetLineColor(2);
  hist09pb->Draw("SAMEHIST");
  hist09pc->SetLineWidth(2);
  hist09pc->SetLineColor(4);
  hist09pc->Draw("SAMEHIST");
  hist09pd->SetLineWidth(2);
  hist09pd->SetLineColor(6);
  hist09pd->Draw("SAMEHIST");
  hist09pe->SetLineWidth(2);
  hist09pe->SetLineColor(8);
  hist09pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist09pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist09pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist09pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist09pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist09pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c09->Print("p_EleVeto.png");

  TCanvas *c10 = new TCanvas("c10","mvaID value",1200,900);
  xaxis = hist10pa->GetXaxis();
  yaxis = hist10pa->GetYaxis();
  xaxis->SetTitle("mvaID Value");
  yaxis->SetTitle("Entries / 0.06");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,1000000);
  c10->SetLogy();
  c10->cd();
  hist10pa->SetLineWidth(2);
  hist10pa->SetLineColor(1);
  hist10pa->Draw("HIST");
  hist10pb->SetLineWidth(2);
  hist10pb->SetLineColor(2);
  hist10pb->Draw("SAMEHIST");
  hist10pc->SetLineWidth(2);
  hist10pc->SetLineColor(4);
  hist10pc->Draw("SAMEHIST");
  hist10pd->SetLineWidth(2);
  hist10pd->SetLineColor(6);
  hist10pd->Draw("SAMEHIST");
  hist10pe->SetLineWidth(2);
  hist10pe->SetLineColor(8);
  hist10pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist10pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist10pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist10pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist10pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist10pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c10->Print("p_mvaVal.png");

  TCanvas *c11 = new TCanvas("c11","mvaID category",1200,900);
  xaxis = hist11pa->GetXaxis();
  yaxis = hist11pa->GetYaxis();
  xaxis->SetTitle("mvaID Category");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,1000000);
  c11->SetLogy();
  c11->cd();
  hist11pa->SetLineWidth(2);
  hist11pa->SetLineColor(1);
  hist11pa->Draw("HIST");
  hist11pb->SetLineWidth(2);
  hist11pb->SetLineColor(2);
  hist11pb->Draw("SAMEHIST");
  hist11pc->SetLineWidth(2);
  hist11pc->SetLineColor(4);
  hist11pc->Draw("SAMEHIST");
  hist11pd->SetLineWidth(2);
  hist11pd->SetLineColor(6);
  hist11pd->Draw("SAMEHIST");
  hist11pe->SetLineWidth(2);
  hist11pe->SetLineColor(8);
  hist11pe->Draw("SAMEHIST");
  legend1->Clear();
  legend1->AddEntry(hist11pa,histotitle+", Photon TightID and EleVeto","f");
  legend1->AddEntry(hist11pb,histotitle+", HLT_Photon175","f");
  legend1->AddEntry(hist11pc,histotitle+", HLT_Photon200","f");
  legend1->AddEntry(hist11pd,histotitle+", HLT_Photon150","f");
  legend1->AddEntry(hist11pe,histotitle+", HLT_Photon165","f");
  legend1->Draw();
  c11->Print("p_mvaCat.png");

  TCanvas *ceff = new TCanvas("ceff","Trigger Efficiency",1200,900);
  ceff->cd();
  Eff1->SetTitle("Trigger Efficiency");
  Eff1->SetLineWidth(2);
  Eff1->SetLineColor(2);
  Eff1->Draw();
  gPad->Update();
  Eff1->SetTitle("Trigger Efficiency; pt_{#gamma} (GeV); Efficiency"); 
  //Eff1->GetPaintedGraph()->GetXaxis()->SetTitle("pt_{#gamma} (GeV)");
  //Eff1->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");
  Eff1->GetPaintedGraph()->GetYaxis()->SetRangeUser(-0.1,1.1);
  Eff1->Draw("SAME");
  Eff2->SetLineWidth(2);
  Eff2->SetLineColor(4);
  Eff2->Draw("SAME");
  Eff3->SetLineWidth(2);
  Eff3->SetLineColor(6);
  Eff3->Draw("SAME");
  Eff4->SetLineWidth(2);
  Eff4->SetLineColor(8);
  Eff4->Draw("SAME");
  legend2->Clear();
  legend2->AddEntry(Eff1,"HLT_Photon175_","f");
  legend2->AddEntry(Eff2,"HLT_Photon200_","f");
  legend2->AddEntry(Eff3,"HLT_Photon150_","f");
  legend2->AddEntry(Eff4,"HLT_Photon165_","f");
  legend2->Draw();
  ceff->Print("eff.png");

  gBenchmark->Show("selectWG");

}
