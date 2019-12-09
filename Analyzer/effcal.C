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
#include <algorithm>
#include <map>
#endif

void effcal(int signal = 700)
{

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);

  int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
  // Open input file
  Float_t s_mass;
  TString signame = std::to_string(signal);
  TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC"+signame+"_WGamma_full_full.root");
  TTree* theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("sys_invmass", &s_mass);
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    count1++;
    if(s_mass < 0.75*signal || s_mass > 1.25*signal) continue;
    count2++;
  }
  cout<<"Entries: "<<count1<<endl;

  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/fullcutdataset/Signal"+signame+"W_Wwindow_full_finalcut.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("m", &s_mass);
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);

    count3++;
    //narrow signals
    if(s_mass < 0.75*signal || s_mass > 1.25*signal) continue;
    count4++;
  }
  cout<<"Entries: "<<count2<<endl;

  cout<<"Mass is: "<<signal<<endl;
  cout<<"Level 1 eff: "<<double(count1) / 20000<<endl;
  cout<<"Level 2 eff: "<<double(count4) / double(count2)<<endl;
  cout<<"Level 3 eff: "<<double(count3) / 20000<<endl;
  cout<<"Rate: "<<(double(count1) / 20000)*(double(count4) / double(count2))*41.527<<" "<<double(count3) / 20000 *41.527<<endl;
  cout<<"TotalEff: "<<(double(count1) / 20000)*(double(count4) / double(count2))<<endl;
  cout<<"=============================="<<endl;
}
