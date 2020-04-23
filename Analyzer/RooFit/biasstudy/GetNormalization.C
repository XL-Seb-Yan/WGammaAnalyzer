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
#include <RooBukinPdf.h>
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

void GetNormalization(int signalmass = 2200, int yhi = 400)
{
  //gErrorIgnoreLevel = kInfo;
  using namespace std;
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);
  
  TString signalmass_str = std::to_string(signalmass);
  
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/rootfit_workspace/anchor_W/"+signalmass_str+"W-shapes-Unbinned-CB2Gaus.root");
  RooWorkspace *w = (RooWorkspace*)file.Get("w");
  RooDataHist *datah = (RooDataHist*)w->data("signal_MC");
  RooAbsPdf *model = w->pdf("CB2Gaus");
  // w->var("CB_alpha")->setConstant();
  // w->var("CB_mean")->setConstant();
  // w->var("CB_n")->setConstant();
  // w->var("CB_sigma")->setConstant();
  // w->var("Gaus_mean")->setConstant();
  // w->var("Gaus_sigma")->setConstant();
  // w->var("frac")->setConstant();
  
  w->var("CB_alpha")->setConstant();
  w->var("CB_mean")->setConstant();
  w->var("CB_n")->setConstant();
  w->var("CB_sigma")->setConstant();
  w->var("Gaus_sigma_1")->setConstant();
  w->var("Gaus_sigma_2")->setConstant();
  w->var("frac1")->setConstant();
  w->var("frac2")->setConstant();
  
  w->var("m")->setRange(600,4000);
  
  // --- Perform extended ML fit ---
  RooRealVar *sig_norm = new RooRealVar("sig_norm","sig_norm",2000,0,5000,"");
  RooAddPdf *ex_model = new RooAddPdf("extended","extended",RooArgList(*model),RooArgList(*sig_norm));
  RooFitResult *ex_r = NULL;
  ex_r = ex_model->fitTo(*datah,Range(signalmass*0.75,signalmass*1.25),RooFit::Minimizer("Minuit2"),Extended(true),SumW2Error(false),Save());
  cout<<"Normalization: "<<sig_norm->getVal()<<" "<<signalmass_str<<endl;
}
