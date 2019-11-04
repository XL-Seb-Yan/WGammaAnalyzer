#include <TMath.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <cstring>
#include <TCanvas.h>
void plottoys()
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("sys_invmass","invmass",600,2500,""); //the name "sys_invmass" will be used by RooDataSet to import data

  // --- Import toy dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.mH120.123456.root");

  TDirectory *d = (TDirectory*)file.Get("toys");

  // --- Histogram to store chi2 test results
  TCanvas *c1 = new TCanvas("data","data",1200,900);
  
  for(int itoy=1; itoy<501; itoy++){

    RooArgList imarglist(*x);
    RooArgSet imargset(imarglist);
    TString toyname = "toy_"+std::to_string(itoy);
    RooDataSet *data = (RooDataSet*)d->Get(toyname);
    cout<<"Processing: ===================                "<<itoy<<"             ==================="<<endl;
    //RooDataSet data = ("Data sideband toys","Data sideband toys",imargset,Import(*import,"toy"));//import branches with names match the "variable name" (not variable) listed in imargset

    // --- Plot ---
    x->setBins(95); //fit is unbinned but chi2 is calculated by binning data with this value
    //gStyle->SetOptStat(111111);
    RooPlot *frame = x->frame(Title("Data Sideband"));
    data->plotOn(frame,RooFit::Name("data"));
    if(itoy%10 == 0){
      c1->Clear();
      c1->cd();
      frame->Draw();
      c1->SetLogy();
      c1->Print(type+" "+toyname+".png");
    }
  }
}
