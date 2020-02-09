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
#include <RooGaussian.h>
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

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

void make_signal_wide_gen_shapes(int signalmass = 3450)
{
  //int rlo = 1200;
  //int rhi = 2100;
  double yhi = 0.1;
  //gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace std;
  using namespace RooFit;
  //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);
  TString signalmass_str = std::to_string(signalmass);

  // --- Create obervable --- 
  RooRealVar *m = new RooRealVar("m","m",signalmass*0.75,signalmass*1.25,""); //the name "m" will be used by RooDataSet to import data
  m->setBins((int)signalmass*0.5 / (int)20);
  //RooRealVar *m = new RooRealVar("m","m",600,875,""); //the name "m" will be used by RooDataSet to import data
  //m->setBins(14);

  
  //--- signal PDF ---
  RooRealVar* CB_mean = new RooRealVar("CB_mean","CB_mean",signalmass,signalmass-100,signalmass+100,"");
  RooRealVar* CB_sigma = new RooRealVar("CB_sigma","CB_sigma",50,10,250,"");
  RooRealVar* CB_alpha = new RooRealVar("CB_alpha","CB_alpha",1,0.1,5,"");
  RooRealVar* CB_n = new RooRealVar("CB_n","CB_n",3,0.1,10,"");
  /*
  RooRealVar* CB_sigma = new RooRealVar("CB_sigma","CB_sigma",85,80,86.85,"");
  RooRealVar* CB_alpha = new RooRealVar("CB_alpha","CB_alpha",1,0.75,1.21,"");
  RooRealVar* CB_n = new RooRealVar("CB_n","CB_n",0.8,0.5,1.14,"");
  */
  RooCBShape* CB_model = new RooCBShape("CBShape","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
  
  
  //--- signal PDF ---
  
  //RooRealVar* Gaus_mean = new RooRealVar("Gaus_mean","Gaus_mean",signalmass,signalmass-100,signalmass+100,"");
  RooRealVar* Gaus_sigma = new RooRealVar("Gaus_sigma","Gaus_sigma",150,50,500,"");
  RooGaussian* Gaus_model = new RooGaussian("Gaussian","Gaussian Function",*m,*CB_mean,*Gaus_sigma);
  
  /*
  RooRealVar* exp_p0 = new RooRealVar("exp_p0","exp_p0",0,-20,20,"");
  RooRealVar* exp_p1 = new RooRealVar("exp_p1","exp_p1",0,-20,20,"");
  RooGenericPdf* exp_model = new RooGenericPdf("tail","exp(exp_p0*m)",RooArgList(*m,*exp_p0));
  */

  RooRealVar* frac = new RooRealVar("frac","frac",0.7,0.5,1);
  RooAddPdf* com_model = new RooAddPdf("composite","composite",RooArgList(*CB_model,*Gaus_model),RooArgList(*frac));

  /*
  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/fullcut/Signal"+signalmass_str+"_Wwindow_sigrange_finalcut.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(*m);
  RooArgSet imargset(imarglist);
  RooDataSet data("Signal","Signal"+signalmass_str,imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset
  */
  /*
  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/RooFitWorkspace_Jan12/GenSignalDataset/wide/roodataset_signal-"+signalmass_str+"-wide.root");
  RooDataSet *input = (RooDataSet*)file.Get("lmorphData");
  RooDataSet data("Signal","Signal"+signalmass_str,input,RooArgSet(*m));
  */
  
  // --- Import Binned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/RooFitWorkspace_Jan12/GenSignalDataset/wide/roodataset_signal-"+signalmass_str+"-wide.root");
  TH1D *MChist = (TH1D*)file.Get("distribs_5_10_0__m");
  RooDataHist datah("Signal","Signal"+signalmass_str,RooArgSet(*m),MChist);
  
  
  
  // --- Perform ML fit of composite PDF to data ---
  /*
  RooFitResult *r = NULL;
  if(signalmass*0.7 > 600 && signalmass*1.3 < 3500)
    r = com_model->fitTo(data,Range(signalmass*0.7,signalmass*1.3),RooFit::Minimizer("Minuit2"),Save());
  else if(signalmass*0.7 < 600)
    r = com_model->fitTo(data,Range(600,signalmass*1.3),RooFit::Minimizer("Minuit2"),Save());
  else if(signalmass*1.3 > 3500)
    r = com_model->fitTo(data,Range(signalmass*0.7,3500),RooFit::Minimizer("Minuit2"),Save());
  */
  RooFitResult *r = com_model->fitTo(datah,Range(signalmass*0.75,signalmass*1.25),RooFit::Minimizer("Minuit2"),SumW2Error(false),Save());

  cout<<"OK with fir"<<endl;

  
  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = m->frame(Title("Signal"+signalmass_str));
  datah.plotOn(frame,RooFit::Name("data"));
  com_model->plotOn(frame,LineColor(2),RooFit::Name("fit"));
  //com_model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  //com_model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  com_model->plotOn(frame,Components("CBShape"),LineStyle(kDashed),LineColor(1),RooFit::Name("CB"));
  com_model->plotOn(frame,Components("Gaussian"),LineStyle(kDashed),RooFit::Name("Gaus"));
  com_model->plotOn(frame,LineColor(2),RooFit::Name("fit"));
  datah.plotOn(frame);

  RooChi2Var chi2 ("chi2", "chi2", *com_model,datah,DataError(RooAbsData::Poisson));//Default: SumW2
  TString chi2txt = "Chi2/ndf: "+to_str_trim(chi2.getVal()/datah.numEntries());
  TString NLLtxt = "NLL: "+to_str_trim(com_model->createNLL(datah)->getVal());
  TString CBmean = "#mu_{CB}: "+to_str_trim(CB_mean->getVal())+" #pm "+to_str_trim(CB_mean->getError())+" (GeV)";
  TString CBsigma = "#sigma_{CB}: "+to_str_trim(CB_sigma->getVal())+" #pm "+to_str_trim(CB_sigma->getError())+" (GeV)";
  TString CBalpha = "#alpha_{CB}: "+to_str_trim(CB_alpha->getVal())+" #pm "+to_str_trim(CB_alpha->getError());
  TString CBn = "n_{CB}: "+to_str_trim(CB_n->getVal())+" #pm "+to_str_trim(CB_n->getError());
  TString Gaussigma = "#sigma_{Gaus}: "+to_str_trim(Gaus_sigma->getVal())+" #pm "+to_str_trim(Gaus_sigma->getError())+" (GeV)";
  TString Fraction = "CB frac: "+to_str_trim(frac->getVal())+" #pm "+to_str_trim(frac->getError());
  double position = 0.72;
  TLatex *chi2lax = new TLatex(0.7,0.7-0.03,chi2txt);
  //TLatex *NLLlax = new TLatex(0.7,0.7-0.06,NLLtxt);
  TLatex *CBmeanlax = new TLatex(0.7,0.7-0.06,CBmean);
  TLatex *CBsigmalax = new TLatex(0.7,0.7-0.09,CBsigma);
  TLatex *CBalphalax = new TLatex(0.7,0.7-0.12,CBalpha);
  TLatex *CBnlax = new TLatex(0.7,0.7-0.15,CBn);
  TLatex *Gaussigmalax = new TLatex(0.7,0.7-0.18,Gaussigma);
  TLatex *Fraclax = new TLatex(0.7,0.7-0.21,Fraction);

  /*
  // --- Perform extended ML fit ---
  RooRealVar *sig_norm = new RooRealVar("sig_norm","sig_norm",5000,0,100000,"");
  RooAddPdf *ex_model = new RooAddPdf("extended","extended",RooArgList(*com_model),RooArgList(*sig_norm));
  RooFitResult *ex_r = NULL;
  ex_r = ex_model->fitTo(data,RooFit::Minimizer("Minuit2"),Extended(true),Save());
  cout<<"Normalization is: "<<sig_norm->getVal()<<endl;
  TString sigN = "Ns: "+to_str_trim(sig_norm->getVal());
  TLatex *sigNlax = new TLatex(0.7,0.38,sigN);
  */

  chi2lax->SetNDC();
  //NLLlax->SetNDC();
  CBmeanlax->SetNDC();
  CBsigmalax->SetNDC();
  CBalphalax->SetNDC();
  CBnlax->SetNDC();
  Gaussigmalax->SetNDC();
  Fraclax->SetNDC();
  //sigNlax->SetNDC();
  
  chi2lax->SetTextSize(0.020);
  //NLLlax->SetTextSize(0.020);
  CBmeanlax->SetTextSize(0.020);
  CBsigmalax->SetTextSize(0.020);
  CBalphalax->SetTextSize(0.020);
  CBnlax->SetTextSize(0.020);
  Gaussigmalax->SetTextSize(0.020);
  Fraclax->SetTextSize(0.020);
  //sigNlax->SetTextSize(0.020);
  
  frame->addObject(chi2lax);
  //frame->addObject(NLLlax);
  frame->addObject(CBmeanlax);
  frame->addObject(CBsigmalax);
  frame->addObject(CBalphalax);
  frame->addObject(CBnlax);
  frame->addObject(Gaussigmalax);
  frame->addObject(Fraclax);
  //frame->addObject(sigNlax);
  
  

  TCanvas *c01 = new TCanvas("c01","c01",1200,900);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 20 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0,yhi);
  xaxis->SetRangeUser(signalmass*0.75,signalmass*1.25);
  //xaxis->SetRangeUser(600,875);
  //frame->SetMaximum(600);
  //frame->SetMinimum(0);
  //c01->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.7,0.75,0.85,0.85);
  l->AddEntry(frame->findObject("fit"),"SignalMC fit","l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  //l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  //l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
  c01->Print(signalmass_str+"line.png");

  TCanvas *c02 = new TCanvas("c02","c02",1200,300);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  c02->cd();
  RooHist* hpull = frame->pullHist();
  RooPlot* pull_frame = m->frame();
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
  xaxis->SetRangeUser(signalmass*0.7,signalmass*1.3);
  //xaxis->SetRangeUser(600,875);
  c02->SetBottomMargin(0.18);
  c02->SetTopMargin(0.18);
  c02->SetGrid();
  pull_frame->Draw();
  c02->Update();
  c02->Print("pull"+signalmass_str+".png");

  TCanvas *c03 = new TCanvas("c03","c03",1200,900);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  xaxis = frame->GetXaxis();
  yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 20 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.0001,0.1);
  xaxis->SetRangeUser(signalmass*0.75,signalmass*1.25);
  //xaxis->SetRangeUser(600,875);
  //xaxis->SetRangeUser(600,signalmass*1.3);
  //frame->SetMaximum(600);
  //frame->SetMinimum(0);
  c03->SetLogy();
  frame->Draw();
  l =  new TLegend(0.7,0.75,0.85,0.85);
  l->AddEntry(frame->findObject("fit"),"SignalMC fit","l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  //l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  //l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
  c03->Print(signalmass_str+"log.png");
  
  // --- Output root file ---
  RooWorkspace *w = new RooWorkspace("w","w");
  w->import(*m);
  w->import(datah,Rename("signal_MC"));
  w->import(*com_model);
  w->writeToFile(signalmass_str+"-shapes-Unbinned-"+".root");
  
  
}
