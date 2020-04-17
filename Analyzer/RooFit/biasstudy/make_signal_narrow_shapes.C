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
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

void make_signal_narrow_shapes(int signalmass = 900, int yhi = 900)
{
  //gErrorIgnoreLevel = kInfo;
  using namespace std;
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);
  
  gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
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
  RooRealVar *m = new RooRealVar("m","m",600,4000,""); //the name "m" will be used by RooDataSet to import data
  
  //--- signal PDF ---
  TString fun_name = "CB";
  RooRealVar* mean = new RooRealVar("mean","mean",signalmass,signalmass-100,signalmass+100,"");
  RooRealVar* alpha = new RooRealVar("alpha","alpha",1,0.1,2.00,"");
  // RooRealVar* sigma = new RooRealVar("sigma","sigma",50,20,150,"");
  // RooRealVar* n = new RooRealVar("n","n",1,0.1,50,"");
  RooRealVar* sigma = new RooRealVar("sigma","sigma",39,36.97,45,"");
  RooRealVar* n = new RooRealVar("n","n",0.1,0,5,"");
  RooCBShape* model = new RooCBShape("CBShape","Cystal Ball Function",*m,*mean,*sigma,*alpha,*n);

  // --- Import Binned dataset ---
  float s_mass, xsec_puweight;
  TH1F MChist("MC","MC",170,600,4000);
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/fullcut/SignalMC"+signalmass_str+"N_postproc_WGamma17_sigrange_WB_fullcut_jmcorr_Mar17.root");
  TTree* tree = (TTree*)file.Get("Events");
  tree->SetBranchAddress("m", &s_mass);
  tree->SetBranchAddress("xsec_puweight", &xsec_puweight);
  for (int ievt = 0; ievt<tree->GetEntries();ievt++) {
    tree->GetEntry(ievt);
    MChist.Fill(s_mass, xsec_puweight);
  }
  RooDataHist datah("Signal MC (W band)","Signal MC (W band)",RooArgSet(*m),&MChist);
  cout<<"number of weighted entries: "<<datah.sum(false)<<endl;
  
  // --- Perform extended ML fit of composite PDF to toy data ---
  RooFitResult *r = model->fitTo(datah,Range(signalmass*0.75,signalmass*1.25),RooFit::Minimizer("Minuit2"),SumW2Error(false),Save()); //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)
  
  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = m->frame(Title("Signal"+signalmass_str));
  datah.plotOn(frame,RooFit::Name("data"));
  // Change attributes of last added plot elements
  frame->getAttMarker()->SetMarkerSize(2);
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name),Range(signalmass*0.75,signalmass*1.25));
  //model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"),Range(signalmass*0.75,signalmass*1.25));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name),Range(signalmass*0.75,signalmass*1.25));
  datah.plotOn(frame);
  
  TString chi2txt = "#chi^{2}/ndf: "+to_str_trim(frame->chiSquare(fun_name,"data", 4));
  TString NLLtxt = "NLL: "+to_str_trim(model->createNLL(datah)->getVal());
  TString CBmean = "#mu_{CB}: "+to_str_trim(mean->getVal())+" #pm "+to_str_trim(mean->getError())+" (GeV)";
  TString CBsigma = "#sigma_{CB}: "+to_str_trim(sigma->getVal())+" #pm "+to_str_trim(sigma->getError())+" (GeV)";
  TString CBalpha = "#alpha_{CB}: "+to_str_trim(alpha->getVal())+" #pm "+to_str_trim(alpha->getError());
  TString CBn = "n_{CB}: "+to_str_trim(n->getVal())+" #pm "+to_str_trim(n->getError());
  TLatex *chi2lax = new TLatex(0.15,0.7-0*0.04,chi2txt);
  TLatex *CBmeanlax = new TLatex(0.15,0.7-1*0.04,CBmean);
  TLatex *CBsigmalax = new TLatex(0.15,0.7-2*0.04,CBsigma);
  TLatex *CBalphalax = new TLatex(0.15,0.7-3*0.04,CBalpha);
  TLatex *CBnlax = new TLatex(0.15,0.7-4*0.04,CBn);

  // --- Perform extended ML fit ---
  // RooRealVar *sig_norm = new RooRealVar("sig_norm","sig_norm",2000,0,5000,"");
  // RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*sig_norm));
  // RooFitResult *ex_r = NULL;
  // ex_r = ex_model->fitTo(data,RooFit::Minimizer("Minuit2"),Extended(true),SumW2Error(false),Save());
  // cout<<"Normalization is: "<<sig_norm->getVal()<<endl;
  // TString sigN = "Ns: "+to_str_trim(sig_norm->getVal());
  // TLatex *sigNlax = new TLatex(0.15,0.7-5*0.04,sigN);

  //chi2lax->SetNDC();
  CBmeanlax->SetNDC();
  CBsigmalax->SetNDC();
  CBalphalax->SetNDC();
  CBnlax->SetNDC();
  
  //chi2lax->SetTextSize(0.026);
  CBmeanlax->SetTextSize(0.026);
  CBsigmalax->SetTextSize(0.026);
  CBalphalax->SetTextSize(0.026);
  CBnlax->SetTextSize(0.026);
  
  //frame->addObject(chi2lax);
  frame->addObject(CBmeanlax);
  frame->addObject(CBsigmalax);
  frame->addObject(CBalphalax);
  frame->addObject(CBnlax);
  

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
  yaxis->SetTitle("Events / 20 GeV");
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0,yhi);
  xaxis->SetRangeUser(signalmass*0.75,signalmass*1.25);
  frame->Draw();
  CMS_lumi(p01a,iPeriod,iPos);
  TLegend *l =  new TLegend(0.65,0.75,0.88,0.88);
  l->AddEntry(frame->findObject("data"),"2017 Signal MC","lep");
  l->AddEntry(frame->findObject(fun_name),"Signal fit "+fun_name,"l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
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
  xaxis->SetRangeUser(signalmass*0.75,signalmass*1.25);
  p01b->SetGrid();
  pull_frame->Draw();
  p01b->Update();
  c01->Print(fun_name+signalmass_str+".png");
  c01->Print(fun_name+signalmass_str+".pdf");
  c01->Print(fun_name+signalmass_str+".root");
  c01->Print(fun_name+signalmass_str+".svg");
  
  
  // --- Output root file ---
  RooWorkspace *w = new RooWorkspace("w","w");
  m->setRange(600,4000);
  m->setBins(170);
  w->import(*m);
  w->import(datah,Rename("signal_MC"));
  w->import(*model);
  w->writeToFile(signalmass_str+"N-shapes-Unbinned-"+fun_name+".root");
  
}
