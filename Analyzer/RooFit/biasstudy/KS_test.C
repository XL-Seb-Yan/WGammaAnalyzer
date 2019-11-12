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
  return ROOT::Math::crystalball_pdf(x,1.24,49.99,26.94,711.62);
}

void KS_test()
{
  //gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace std;
  using namespace RooFit;
  //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  int signalmass = 700;
  TString signalmass_str = std::to_string(signalmass)+"N";

  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/biasstudy/dataset/Signal"+signalmass_str+"_Wwindow_full_finalcut.root");
  TTree* tree = (TTree*)file.Get("Events");

  // GOF test within ROOT
  float s_mass;
  tree->SetBranchAddress("m", &s_mass);
  TH1F *show = new TH1F("fit test","fit test",145,600,3500);
  double *sample = new double[tree->GetEntries()];
  for(int ievt=0; ievt<tree->GetEntries();ievt++){
    tree->GetEntry(ievt);
    sample[ievt] = s_mass;
    show->Fill(s_mass);
  }
  
  ROOT::Math::Functor1D f(&CB_pdf);
  ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(tree->GetEntries(),sample,f,ROOT::Math::GoFTest::kPDF,600,3500);
  double pvalueKS = -99;
  pvalueKS = goftest->KolmogorovSmirnovTest();
  cout<<"p-value of KS test is: "<<pvalueKS<<endl;
}
