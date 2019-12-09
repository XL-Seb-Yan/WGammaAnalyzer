//make signal shapes
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
#include <algorithm>
#include <map>
#include <RooMsgService.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooCBShape.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooChi2Var.h>
#include <TLatex.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooPlotable.h>
#include <RooWorkspace.h>
#include <RooAddPdf.h>
#endif

double CB_pdf(double x){
  return ROOT::Math::crystalball_pdf(x,1.3194,2.7968,80,3513.59);
}
/*
double com_pdf(double x){
  double frac = 0.88;
  return frac*ROOT::Math::crystalball_pdf(x,1.36,0.68,68.5,1601.4) + (1-frac)*ROOT::Math::gaussian_pdf(x,150,1601.4);
}
*/

void KS_test()
{

  //rlo{600,600,600,800,1000,1100,1350,1500,1700,1800,2000,2200,2400}
  //rhi{1100,1200,1400,1600,1800,2100,2300,2600,2800,3000,3200,3400,3500}

  int rlo = 3500*0.7;
  int rhi = 3500*1.3;
  //gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace std;
  using namespace RooFit;
  //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  int signalmass = 3500;
  TString signalmass_str = std::to_string(signalmass)+"N";

  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/fullcutdataset/Signal"+signalmass_str+"_Wwindow_full_finalcut.root");
  TTree* tree = (TTree*)file.Get("Events");

  // GOF test within ROOT
  float s_mass;
  tree->SetBranchAddress("m", &s_mass);
  TH1F *show = new TH1F("fit test","fit test",145,600,3500);
  int N = 0;
  for(int ievt=0; ievt<tree->GetEntries();ievt++){
    tree->GetEntry(ievt);
    //if(s_mass < rlo || s_mass > rhi) continue;
    N++;
  }
  cout<<N<<endl;
  double *sample = new double[N];
  int i = 0;
  for(int ievt=0; ievt<tree->GetEntries();ievt++){
    tree->GetEntry(ievt);
    //if(s_mass < rlo || s_mass > rhi) continue;
    sample[i] = s_mass;
    i++;
    show->Fill(s_mass);
  }
  ROOT::Math::Functor1D f(&CB_pdf);
  ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(tree->GetEntries(),sample,f,ROOT::Math::GoFTest::kPDF,rlo,rhi);
  double pvalueKS = -99;
  pvalueKS = goftest->KolmogorovSmirnovTest();
  cout<<"p-value of KS test is: "<<pvalueKS<<endl;
}
