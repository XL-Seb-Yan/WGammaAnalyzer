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
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

void make_signal_wide_shapes_CB2Gaus(int signalmass = 2800, int yhi = 500)
{
  //gErrorIgnoreLevel = kInfo;
  using namespace std;
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);
  
  gROOT->SetBatch(1);
  lumi_13TeV = "";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 5;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.045,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetBarWidth(1.03);
  gStyle->SetHistLineWidth(2);

  TString signalmass_str = std::to_string(signalmass);

  // --- Create obervable --- 
  RooRealVar *m = new RooRealVar("m","m",600,10000,""); //the name "m" will be used by RooDataSet to import data

  
  //--- signal PDF ---
  TString fun_name = "CB2Gaus";
  RooRealVar* CB_mean = new RooRealVar("CB_mean","CB_mean",signalmass,signalmass-100,signalmass+100,"");
  RooRealVar* CB_sigma = new RooRealVar("CB_sigma","CB_sigma",111,50,150,"");
  RooRealVar* CB_alpha = new RooRealVar("CB_alpha","CB_alpha",1.2,1,2,"");
  RooRealVar* CB_n = new RooRealVar("CB_n","CB_n",1,0.5,1.5,"");
  RooCBShape* CB_model = new RooCBShape("CB","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
  
  // RooRealVar* Gaus_mean_1 = new RooRealVar("Gaus_mean_1","Gaus_mean_1",signalmass,signalmass-100,signalmass+100,"");
  RooRealVar* Gaus_sigma_1 = new RooRealVar("Gaus_sigma_1","Gaus_sigma_1",110,50,150,"");
  RooGaussian* Gaus_model_1 = new RooGaussian("Gaussian_1","Gaussian Function",*m,*CB_mean,*Gaus_sigma_1);
  
  // RooRealVar* Gaus_mean_2 = new RooRealVar("Gaus_mean_2","Gaus_mean_2",signalmass,signalmass-100,signalmass+100,"");
  RooRealVar* Gaus_sigma_2 = new RooRealVar("Gaus_sigma_2","Gaus_sigma_2",300,200,500,"");
  RooGaussian* Gaus_model_2 = new RooGaussian("Gaussian_2","Gaussian Function",*m,*CB_mean,*Gaus_sigma_2);
  
  RooRealVar* frac1 = new RooRealVar("frac1","frac1",0.85,0.75,0.9);
  RooAddPdf* com_model_1 = new RooAddPdf("CBGaus","CBGaus",RooArgList(*CB_model,*Gaus_model_1),RooArgList(*frac1));
  RooRealVar* frac2 = new RooRealVar("frac2","frac2",0.80,0.7,0.85);
  RooAddPdf* com_model = new RooAddPdf("CB2Gaus","CB2Gaus",RooArgList(*com_model_1,*Gaus_model_2),RooArgList(*frac2));


  // --- Import Binned dataset ---
  
  int bin = 376; //25GeV
  if(signalmass > 1999)//50GeV
    bin=188;
  if(signalmass > 4999)//100GeV
    bin=94;
  
  float s_mass, xsec_puweight;
  TH1F MChist("MC","MC",bin,600,10000);
									
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/fullcut/SignalMC"+signalmass_str+"W_postproc_WGamma17_SR_sigrange_fullcut_jmcorr_May22.root");
  TTree* tree = (TTree*)file.Get("Events");
  tree->SetBranchAddress("m", &s_mass);
  tree->SetBranchAddress("xsec_puweight", &xsec_puweight);
  for (int ievt = 0; ievt<tree->GetEntries();ievt++) {
    tree->GetEntry(ievt);
    MChist.Fill(s_mass, xsec_puweight);
  }
  RooDataHist datah("Signal MC (W band)","Signal MC (W band)",RooArgSet(*m),&MChist);
  cout<<"number of weighted entries: "<<datah.sum(false)<<endl;
  
  RooFitResult *r = com_model->fitTo(datah,Range(signalmass*0.75,signalmass*1.25),RooFit::Minimizer("Minuit2"),SumW2Error(false),Save()); //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)

  
  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = m->frame(Title("Signal"+signalmass_str));
  datah.plotOn(frame,RooFit::Name("data"));
  frame->getAttMarker()->SetMarkerSize(2);
  com_model->plotOn(frame,LineColor(2),RooFit::Name("fit"),Range(signalmass*0.75,signalmass*1.25));
  //com_model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  com_model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"),Range(signalmass*0.75,signalmass*1.25));
  com_model->plotOn(frame,Components("CB"),LineStyle(kDashed),LineColor(1),RooFit::Name("CB"),Range(signalmass*0.75,signalmass*1.25));
  com_model->plotOn(frame,Components("Gaussian_1"),LineStyle(kDashed),LineColor(30),RooFit::Name("Gaus1"),Range(signalmass*0.75,signalmass*1.25));
  com_model->plotOn(frame,Components("Gaussian_2"),LineStyle(kDashed),LineColor(42),RooFit::Name("Gaus2"),Range(signalmass*0.75,signalmass*1.25));
  com_model->plotOn(frame,LineColor(2),RooFit::Name("fit"),Range(signalmass*0.75,signalmass*1.25));
  datah.plotOn(frame);
  frame->getAttMarker()->SetMarkerSize(2);
  
  TString chi2txt = "#chi^{2}/ndf: "+to_str_trim(frame->chiSquare("fit","data", 7));
  TString NLLtxt = "NLL: "+to_str_trim(com_model->createNLL(datah)->getVal());
  TString CBmean = "#mu_{CB}: "+to_str_trim(CB_mean->getVal())+" #pm "+to_str_trim(CB_mean->getError())+" (GeV)";
  TString CBsigma = "#sigma_{CB}: "+to_str_trim(CB_sigma->getVal())+" #pm "+to_str_trim(CB_sigma->getError())+" (GeV)";
  TString CBalpha = "#alpha_{CB}: "+to_str_trim(CB_alpha->getVal())+" #pm "+to_str_trim(CB_alpha->getError());
  TString CBn = "n_{CB}: "+to_str_trim(CB_n->getVal())+" #pm "+to_str_trim(CB_n->getError());
  TString Gaussigma1 = "#sigma_{Gaus1}: "+to_str_trim(Gaus_sigma_1->getVal())+" #pm "+to_str_trim(Gaus_sigma_1->getError())+" (GeV)";
  TString Gaussigma2 = "#sigma_{Gaus2}: "+to_str_trim(Gaus_sigma_2->getVal())+" #pm "+to_str_trim(Gaus_sigma_2->getError())+" (GeV)";
  TString Fraction1 = "frac1: "+to_str_trim(frac1->getVal())+" #pm "+to_str_trim(frac1->getError());
  TString Fraction2 = "frac2: "+to_str_trim(frac2->getVal())+" #pm "+to_str_trim(frac2->getError());
  TLatex *chi2lax = new TLatex(0.15,0.7-0*0.04,chi2txt);
  TLatex *NLLlax = new TLatex(0.15,0.8-0.06,NLLtxt);
  TLatex *CBmeanlax = new TLatex(0.15,0.7-1*0.04,CBmean);
  TLatex *CBsigmalax = new TLatex(0.15,0.7-2*0.04,CBsigma);
  TLatex *CBalphalax = new TLatex(0.15,0.7-3*0.04,CBalpha);
  TLatex *CBnlax = new TLatex(0.15,0.7-4*0.04,CBn);
  TLatex *Gaussigma1lax = new TLatex(0.15,0.7-5*0.04,Gaussigma1);
  TLatex *Gaussigma2lax = new TLatex(0.15,0.7-6*0.04,Gaussigma2);
  TLatex *Frac1lax = new TLatex(0.15,0.7-7*0.04,Fraction1);
  TLatex *Frac2lax = new TLatex(0.15,0.7-8*0.04,Fraction2);

  CBmeanlax->SetNDC();
  CBsigmalax->SetNDC();
  CBalphalax->SetNDC();
  CBnlax->SetNDC();
  Gaussigma1lax->SetNDC();
  Frac1lax->SetNDC();
  Gaussigma2lax->SetNDC();
  Frac2lax->SetNDC();

  CBmeanlax->SetTextSize(0.026);
  CBsigmalax->SetTextSize(0.026);
  CBalphalax->SetTextSize(0.026);
  CBnlax->SetTextSize(0.026);
  Gaussigma1lax->SetTextSize(0.026);
  Frac1lax->SetTextSize(0.026);
  Gaussigma2lax->SetTextSize(0.026);
  Frac2lax->SetTextSize(0.026);
  
  frame->addObject(CBmeanlax);
  frame->addObject(CBsigmalax);
  frame->addObject(CBalphalax);
  frame->addObject(CBnlax);
  frame->addObject(Gaussigma1lax);
  frame->addObject(Gaussigma2lax);
  frame->addObject(Frac1lax);
  frame->addObject(Frac2lax);
  
  TCanvas *c01 = new TCanvas("c01","c01",2100,2000);
  TPad *p01a = new TPad("p01a","p01a",0.05,0.27,0.95,1.0);
  TPad *p01b = new TPad("p01b","p01b",0.05,0.10,0.95,0.315);
  p01a->Draw();
  p01b->Draw();
  p01a->cd();
  p01a->SetLeftMargin(0.11);
  p01a->SetRightMargin(0.1);
  p01a->SetBottomMargin(0.11);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 50 GeV");
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0,yhi);
  xaxis->SetLimits(signalmass*0.75,signalmass*1.25);
  frame->Draw();
  CMS_lumi(p01a,iPeriod,iPos);
  TLegend *l =  new TLegend(0.60,0.75,0.88,0.88);
  l->AddEntry(frame->findObject("data"),"2017 Signal MC","lep");
  l->AddEntry(frame->findObject("fit"),"Signal fit "+fun_name,"l");
  l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  //l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
  
  p01b->cd();
  p01b->SetLeftMargin(0.11);
  p01b->SetRightMargin(0.1);
  p01b->SetTopMargin(0.2);
  p01b->SetBottomMargin(0.25);
  RooHist* hpull = frame->pullHist();
  RooPlot* pull_frame = m->frame();
  pull_frame->addPlotable(hpull,"B") ;
  hpull->SetMarkerStyle(8);
  hpull->SetMarkerSize(0);
  hpull->SetLineWidth(0);
  hpull->SetFillColor(kAzure+1);
  hpull->SetFillColorAlpha(kAzure+1,0.7);
  xaxis = pull_frame->GetXaxis();
  yaxis = pull_frame->GetYaxis();
  xaxis->SetTitle("M_{j#gamma} (GeV)");
  yaxis->SetTitle("#frac{MC-fit}{#sigma_{stat}}");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetRangeUser(-5,5);
  xaxis->SetLabelSize(0.12);
  xaxis->SetTitleSize(0.12);
  yaxis->SetLabelSize(0.12);
  yaxis->SetTitleSize(0.12);
  yaxis->SetNdivisions(5);
  xaxis->SetLimits(signalmass*0.75,signalmass*1.25);
  p01b->SetGrid();
  pull_frame->Draw();
  p01b->Update();
  c01->Print("CB2Gaus"+signalmass_str+".png");
  c01->Print("CB2Gaus"+signalmass_str+".pdf");
  c01->Print("CB2Gaus"+signalmass_str+".svg");
  
  // --- Output root file ---
  RooWorkspace *w = new RooWorkspace("w","w");
  w->import(*m);
  w->import(datah,Rename("signal_MC"));
  w->import(*com_model);
  w->writeToFile(signalmass_str+"W-shapes-Unbinned-CB2Gaus"+".root");
  
  
}
