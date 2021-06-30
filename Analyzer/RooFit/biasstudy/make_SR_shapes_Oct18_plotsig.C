#define fun_type 2
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
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

void make_SR_shapes_Oct18_plotsig(int seed=37)
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37); 

  gROOT->SetBatch(1);
  lumi_13TeV = "137 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Preliminary";
  lumiTextSize = 0.4;
  cmsTextSize = 0.5;
  int iPeriod = 12;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.035);
  gStyle->SetBarWidth(1.03);
  gStyle->SetHistLineWidth(2);

  RooRealVar *x = new RooRealVar("m","m",600,7600,"");
  RooPlot *frame = x->frame();
  
  TString plotmass[6]={"1000N","2600N","4000N","1000W","2600W","4000W"};
  
  for(int j=0; j<6; j++){
    std::vector<float> signalpull;
    float sig_norm = 0;
    float lrange = 0;
    float hrange = 0;
    TString width = "N";
    TString sigfun = "CBGaus";
    int linestyle = 1;
    int linecolor = 2;
    if(plotmass[j].Contains("W")){
        width = "W";
        sigfun="CB2Gaus";
        linestyle = 2;
    }
    if(plotmass[j] == "1000N") {sig_norm=356; lrange=750; hrange=1250;} //20fb
    if(plotmass[j] == "1000W") {sig_norm=343; lrange=750; hrange=1250;} //20fb
    if(plotmass[j] == "2600N") {sig_norm=31.6; linecolor=4; lrange=1950; hrange=3250;}//2fb
    if(plotmass[j] == "2600W") {sig_norm=28.84; linecolor=4; lrange=1950; hrange=3250;}//2fb
    if(plotmass[j] == "4000N") {sig_norm=6.90; linecolor=8; lrange=3000; hrange=5000;}//0.5fb for VVdijet1 and ATLAS1. 0.3fb for others
    if(plotmass[j] == "4000W") {sig_norm=6.0425; linecolor=8; lrange=3000; hrange=5000;}//0.5fb
    TFile *signal = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/RooFitWorkspace/"+width+"/"+plotmass[j]+"-shapes-Unbinned-"+sigfun+".root");
    RooWorkspace *sig_w = (RooWorkspace*)signal->Get("w");
    RooAbsPdf *sig_pdf = sig_w->pdf("CBGaus");
    RooRealVar* m = sig_w->var("m");
    sig_pdf->plotOn(frame,Range(lrange,hrange),LineStyle(linestyle),LineColor(linecolor),Normalization(100,RooAbsReal::NumEvent));
  }

  // --- Visualization ---
  gStyle->SetOptStat(111111);
  TCanvas *c01 = new TCanvas("c01","c01",2100,2000);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{J#gamma} [GeV]");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 40 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.0007,100000);
  xaxis->SetLimits(600,7600);
  c01->SetLogy();
  c01->SetLogx();
  
  std::vector<TLine*> xticks;
  for(int i=700; i<7700; i+=100){
    TLine *line = NULL;
    if(i%1000 != 0)
      line = new TLine(i,0.0007,i,0.0013);
    else
      line = new TLine(i,0.0007,i,0.0018);
    xticks.push_back(line);
  }
  for(int i=0; i<xticks.size(); i++){
    if(xticks[i] != NULL){
      xticks[i]->SetLineWidth(2);
      xticks[i]->Draw();
    }
  }
  
  TLine *line0 = new TLine(700,0.5,6000,0.5);
  line0->Draw();
  
  
  //frame->Draw();


  c01->Print("sig.png");
  c01->Print("sig.pdf");
	c01->Print("sig.svg");



  /*
  // --- Perform extended ML fit of composite PDF to toy data ---
  RooFitResult *ex_r = NULL;
  if(isNorm)
    ex_r = ex_model->fitTo(data_norm,Range(600,3000),RooFit::Minimizer("Minuit2"),Extended(true),SumW2Error(false),Save());
  else
    ex_r = ex_model->fitTo(data,Range(600,3000),RooFit::Minimizer("Minuit2"),Extended(true),Save());
  cout<<"Normalization is: "<<bkg_norm->getVal()<<endl;
  */

}