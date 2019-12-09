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
#include "Math/GoFTest.h"
#include "Math/Functor.h"
#include "Math/DistFunc.h"
#endif

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

double CB_pdf(double x, double mean, double sigma, double alpha, double n){
  return ROOT::Math::crystalball_pdf(x, alpha, n, sigma, mean);
}

void make_signal_shapes(int signalmass = 1800)
{
  //gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace std;
  using namespace RooFit;
  //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);
  TString signalmass_str = std::to_string(signalmass)+"N";

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("m","m",signalmass*0.7,signalmass*1.3,""); //the name "m" will be used by RooDataSet to import data
  x->setBins(int(signalmass*0.06));

  //--- signal PDF ---
  TString fun_name = "CB";
  RooRealVar* sqrtS = new RooRealVar("sqrtS","sqrtS",13000,"");
  sqrtS->setConstant(kTRUE);
  RooFormulaVar* xx = new RooFormulaVar("invmass/sqrtS","m / sqrtS",RooArgSet(*x,*sqrtS));
  RooRealVar* mean = new RooRealVar("mean","mean",signalmass,signalmass-100,signalmass+100,"");
  RooFormulaVar* xmean = new RooFormulaVar("CB_mean","mean / sqrtS",RooArgSet(*mean,*sqrtS));
  RooRealVar* sigma = new RooRealVar("sigma","sigma",30,0,80,"");
  RooFormulaVar* xsigma = new RooFormulaVar("CB_sigma","sigma / sqrtS",RooArgSet(*sigma,*sqrtS));
  RooRealVar* alpha = new RooRealVar("alpha","alpha",1,0,5,"");
  RooRealVar* n = new RooRealVar("n","n",5,0,10,"");
  RooCBShape* model = new RooCBShape("CBShape","Cystal Ball Function",*xx,*xmean,*xsigma,*alpha,*n);
  
  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/fullcutdataset/Signal"+signalmass_str+"_Wwindow_full_finalcut.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(*x);
  RooArgSet imargset(imarglist);
  RooDataSet data("Signal","Signal"+signalmass_str,imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  
  // --- Perform ML fit of composite PDF to data ---
  RooFitResult *r = model->fitTo(data,RooFit::Minimizer("Minuit2"),Save());

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = x->frame(Title("Signal"+signalmass_str));
  data.plotOn(frame,RooFit::Name("data"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  data.plotOn(frame);
  RooDataHist datah("dh","binned data",RooArgSet(*x),data);
  RooChi2Var chi2 ("chi2", "chi2", *model,datah,DataError(RooAbsData::Poisson));//Default: SumW2
  TString chi2txt = "Chi2: "+to_str_trim(chi2.getVal());
  TString NLLtxt = "NLL: "+to_str_trim(model->createNLL(data)->getVal());
  TString CBmean = "#mu_{CB}: "+to_str_trim(mean->getVal())+" #pm "+to_str_trim(mean->getError())+" (GeV)";
  TString CBsigma = "#sigma_{CB}: "+to_str_trim(sigma->getVal())+" #pm "+to_str_trim(sigma->getError())+" (GeV)";
  TString CBalpha = "#alpha_{CB}: "+to_str_trim(alpha->getVal())+" #pm "+to_str_trim(alpha->getError());
  TString CBn = "n_{CB}: "+to_str_trim(n->getVal())+" #pm "+to_str_trim(n->getError());
  TLatex *chi2lax = new TLatex(0.7,0.6,chi2txt);
  TLatex *NLLlax = new TLatex(0.7,0.57,NLLtxt);
  TLatex *CBmeanlax = new TLatex(0.7,0.54,CBmean);
  TLatex *CBsigmalax = new TLatex(0.7,0.51,CBsigma);
  TLatex *CBalphalax = new TLatex(0.7,0.48,CBalpha);
  TLatex *CBnlax = new TLatex(0.7,0.45,CBn);

  // --- Perform extended ML fit ---
  RooRealVar *sig_norm = new RooRealVar("sig_norm","sig_norm",5000,0,100000,"");
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*sig_norm));
  RooFitResult *ex_r = NULL;
  ex_r = ex_model->fitTo(data,RooFit::Minimizer("Minuit2"),Extended(true),Save());
  cout<<"Normalization is: "<<sig_norm->getVal()<<endl;
  TString sigN = "Ns: "+to_str_trim(sig_norm->getVal());
  TLatex *sigNlax = new TLatex(0.7,0.42,sigN);

  chi2lax->SetNDC();
  NLLlax->SetNDC();
  CBmeanlax->SetNDC();
  CBsigmalax->SetNDC();
  CBalphalax->SetNDC();
  CBnlax->SetNDC();
  sigNlax->SetNDC();
  
  chi2lax->SetTextSize(0.022);
  NLLlax->SetTextSize(0.022);
  CBmeanlax->SetTextSize(0.022);
  CBsigmalax->SetTextSize(0.022);
  CBalphalax->SetTextSize(0.022);
  CBnlax->SetTextSize(0.022);
  sigNlax->SetTextSize(0.022);
  
  frame->addObject(chi2lax);
  frame->addObject(NLLlax);
  frame->addObject(CBmeanlax);
  frame->addObject(CBsigmalax);
  frame->addObject(CBalphalax);
  frame->addObject(CBnlax);
  frame->addObject(sigNlax);

  TCanvas *c01 = new TCanvas("c01","c01",1200,900);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 10 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0,300);
  //frame->SetMaximum(600);
  //frame->SetMinimum(0);
  //c01->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.7,0.75,0.85,0.85);
  l->AddEntry(frame->findObject(fun_name),"SignalMC fit "+fun_name,"l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
  c01->Print(fun_name+signalmass_str+".png");

  TCanvas *c02 = new TCanvas("c02","c02",1200,300);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  c02->cd();
  RooHist* hpull = frame->pullHist();
  RooPlot* pull_frame = x->frame();
  pull_frame->addPlotable(hpull,"P") ;
  hpull->SetMarkerStyle(8);
  xaxis = pull_frame->GetXaxis();
  yaxis = pull_frame->GetYaxis();
  xaxis->SetTitle("M_{j#gamma}");
  yaxis->SetTitle("#frac{MC-fit}{#sigma_{stat}}");
  yaxis->SetTitleOffset(0.5);
  yaxis->SetRangeUser(-5,5);
  xaxis->SetLabelSize(0.08);
  xaxis->SetTitleSize(0.08);
  yaxis->SetLabelSize(0.08);
  yaxis->SetTitleSize(0.08);
  yaxis->SetNdivisions(5);
  c02->SetBottomMargin(0.18);
  c02->SetTopMargin(0.18);
  c02->SetGrid();
  pull_frame->Draw();
  c02->Update();
  c02->Print("pull"+fun_name+signalmass_str+".png");
  
  // --- Output root file ---
  
  RooWorkspace *w = new RooWorkspace("w","w");
  w->import(*x);
  w->import(data,Rename("signal_MC"));
  w->import(*model);
  w->writeToFile(signalmass_str+"-shapes-Unbinned-unfitted"+fun_name+".root");
}
