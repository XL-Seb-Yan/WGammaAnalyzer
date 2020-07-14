#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooBukinPdf.h"
#include "RooIntegralMorph.h"
#include "RooNLLVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TH1.h"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"

using namespace RooFit;

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

void plotPDF(){

  gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 5;
  int iPos = 33;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.015);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);
  
  TFile *file = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/Analyzer/CMSSW_9_4_13/src/WGammaAnalyzer/Analyzer/RooFit/biasstudy/700W-shapes-Unbinned-Bukin.root");
  RooWorkspace *w = (RooWorkspace*)file->Get("w");
  w->Print();
  RooRealVar *m = w->var("m");
  RooRealVar *xi = w->var("xi");
  RooRealVar *rho1 = w->var("rho1");
  RooRealVar *rho2 = w->var("rho2");
  m->setRange(00,4000);
  RooPlot *frame = m->frame(Title("Signal"));
  
  TLegend* leg = new TLegend(0.7,0.2,0.9,0.8);
  
  for(int i=0; i<10; i++){
	rho2->setVal(0.1-i*0.02);
	//rho1->setVal(-0.1);
	RooAbsPdf *model = w->pdf("Bukin");
	TString name = "rho2"+to_str_trim(0.1-i*0.02);
	model->plotOn(frame,Range(0,4000,kTRUE),RooFit::Name(name),LineColor(i));
	leg->AddEntry(frame->findObject(name),name,"l");
  }
  
  TCanvas *c = new TCanvas("","",1200,900);
  c->cd();
  frame->Draw();
  leg->Draw();
  c->Print("pdf.png");
}