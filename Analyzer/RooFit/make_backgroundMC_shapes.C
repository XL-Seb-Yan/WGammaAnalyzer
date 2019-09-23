#include <TMath.h>
#include <TLegend.h>
void make_backgroundMC_shapes(int seed=37)
{
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37); 
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  // --- Create workspace --- 
  RooWorkspace w("w","w");

  //--- background PDF ---
  /*
  RooRealVar *p0 = new RooRealVar("p0","p_0",1,0,3,"");
  RooRealVar *p1 = new RooRealVar("p1","p_1",10,10,40,"");
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-2,10,"");
  RooRealVar *Nbkg = new RooRealVar("pol2_Nbkg","N_{bkg}",5000000,0,10000000,"events");
  RooRealVar *x = new RooRealVar("x","x",600,2000);
  RooGenericPdf *dijet = new RooGenericPdf("dijet","p0*x^(p1+p2*log(x))",RooArgList(*x,*p0,*p1,*p2));
  */
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-5,5);
  RooRealVar *x = new RooRealVar("x","x",600,2000);
  RooRealVar *Nbkg = new RooRealVar("pol2_Nbkg","N_{bkg}",5000000,0,10000000,"events");
  RooExponential *expo = new RooExponential("expo","GJets",*x,*p0);
  RooAddPdf *ex_sum = new RooAddPdf("model","model",RooArgList(*expo),RooArgList(*Nbkg));
  w.import(*ex_sum);
 
  
  // --- Import Binned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/Histogram_GJets.root");
  TH1* histo = (TH1*)file.Get("10");
  RooDataHist data("GJets","GJets",*x, histo);

  // --- Perform extended ML fit of composite PDF to toy data ---
  RooFitResult *r = ex_sum->fitTo(data,Save());

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = x->frame();
  data.plotOn(frame);
  frame->SetTitle("WGamma GJets MC");
  //ex_sum->plotOn(frame,LineColor(kRed));
  ex_sum->plotOn(frame,LineStyle(kDashed),RooFit::Name("bkgfun"));
  ex_sum->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  ex_sum->plotOn(frame,LineStyle(kDashed),RooFit::Name("bkgfun"));
  data.plotOn(frame);
  //w.pdf("model_s")->plotOn(frame,Components("Bern"),LineStyle(kDashed));

  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(1,10000000);
  //frame->SetMaximum(600);
  //frame->SetMinimum(0);
  c1->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.6,0.7,0.8,0.78);
  l->AddEntry(frame->findObject("sigfun"),"Signal Fit","l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err"),"Fit Error 1 #sigma","f");
  //l->Draw("same");
  c1->Print("data.png");
  /*
  
  // --- Output root file ---
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("sys_invmass[400,2000]");
  wUP->import(data,Rename("MC"));
  wUP->import(*gauss1);
  wUP->import(*expo);
  wUP->writeToFile("SignalMC1600-shapes-UnbinnedParam.root");
  */
}
