#include <TMath.h>
#include <TLegend.h>
void make_signal_shapes(int seed=37)
{
  using namespace RooFit;
  Double_t signalmass = 900;
  RooRandom::randomGenerator()->SetSeed(37); 
  TCanvas *c1 = new TCanvas("c1","c1",900,1200);
  c1->SetLogy();

  // --- Create workspace --- 
  RooWorkspace w("w","w");
  RooRealVar sys_invmass("sys_invmass","invmass",signalmass*0.5,signalmass*1.5,""); //the name "sys_invmass" will be used by RooDataSet to import data
  sys_invmass.setBins(50);

  //--- signal PDF ---
  RooRealVar* mean = new RooRealVar("CB_mean","mean",signalmass,signalmass-50,signalmass+50,"");
  RooRealVar* CBsigma = new RooRealVar("CB_sigma","CBsigma",40,0,100,"");
  RooRealVar* Gaussigma = new RooRealVar("Gaus_sigma","Gaussigma",40,0,100,"");
  RooRealVar* alpha = new RooRealVar("CB_alpha","alpha",2,0,5,"");
  RooRealVar* n = new RooRealVar("CB_n","n",10,0,15,"");
  RooRealVar* frac = new RooRealVar("Gaus_frac","frac",0.7,0,1); 
  RooRealVar* CBNs = new RooRealVar("CB_Ns","CBN_{s}",1000,0,10000,"events");
  RooRealVar* GausNs = new RooRealVar("Gaus_Ns","GausN_{s}",100,0,10000,"events");
  RooGaussian* Gaus = new RooGaussian("Gaus","Gaus",sys_invmass,*mean,*Gaussigma);
  RooCBShape* CBShape = new RooCBShape("CBShape","Cystal Ball Function",sys_invmass,*mean,*CBsigma,*alpha,*n);
  RooAddPdf* ex_Sum = new RooAddPdf("Sum_ex","",RooArgList(*CBShape, *Gaus),RooArgList(*CBNs,*GausNs));
  w.import(*ex_Sum);

  //--- background PDF ---
  /*
  RooRealVar *pC = new RooRealVar("pol2_pC","C",1,0,100,"");
  pC->setConstant(kFALSE);
  RooRealVar *p0 = new RooRealVar("pol2_p0","p_0",-0.1,-1,0,"");
  p0->setConstant(kFALSE);
  RooRealVar *Nbkg   = new RooRealVar("pol2_Nbkg","N_{bkg}",500,0,2000,"events");
  Nbkg->setConstant(kFALSE);
  RooExponential* expo = new RooExponential("expo","",sys_invmass,*p0);
  RooAddPdf* ex_expo = new RooAddPdf("expo_ex","",RooArgList(*expo),RooArgList(*Nbkg));
  w.import(*ex_expo);
  */

  //--- combine PDF ---
  /*
  RooAddPdf* ex_sum = new RooAddPdf("SUM_ex","",RooArgList(*expo,*gauss1),RooArgList(*Nbkg,*Ns));
  w.import(*ex_sum);
  */
  
  // --- Import unbinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/Signal900N_Wwindow.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(sys_invmass);
  RooArgSet imargset(imarglist);
  RooDataSet data("Signal MC","Signal MC",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  // --- Perform extended ML fit of composite PDF to toy data ---
  //w.pdf("model_s")->fitTo(data);
  //sys_invmass.setRange("SB1",1.47682,1.67682);
  //sysinvmass.setRange("SB2",1.87682,2.07682);
  //RooFitResult *r = ex_bern->fitTo(data,Range("SB1,SB2"),Save());
  RooFitResult *r = ex_Sum->fitTo(data,Save());

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = sys_invmass.frame();
  data.plotOn(frame, RooFit::Name("data"));
  frame->SetTitle("WGamma M-900 MC");
  ex_Sum->plotOn(frame,LineColor(kBlack),RooFit::Name("sigfun"));
  ex_Sum->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  ex_Sum->plotOn(frame,LineColor(kBlack),RooFit::Name("sigfun"));
  ex_Sum->plotOn(frame,Components("CBShape"),LineColor(kRed),RooFit::Name("CBfun"));
  ex_Sum->plotOn(frame,Components("Gaus"),LineColor(kBlue),RooFit::Name("Gausfun"));
  //ex_bern->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //ex_CB->plotOn(frame,Components("expo"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  data.plotOn(frame);
  //w.pdf("model_s")->plotOn(frame,Components("Bern"),LineStyle(kDashed));

  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{W#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events");
  yaxis->SetTitleOffset(1.2);
  frame->SetMaximum(1000);
  frame->SetMinimum(0.1);
  frame->Draw();
  TLegend *l =  new TLegend(0.65,0.75,0.85,0.85);
  l->AddEntry(frame->findObject("sigfun"),"Signal Fit","l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err"),"Fit Error 1 #sigma","f");
  l->Draw("same");
  TText *text = new TText(0.5,0.5,"test");
  text->Draw("same");
  // TString chi2 = Form("%f",frame->chiSquare("sigfun","data"));
  //TText *chi2show = new TText(.65,.35,"12344");
  //chi2show->Draw("SAME");
  c1->Print("data.png");
  cout<<"chisquare is: "<<frame->chiSquare("sigfun","data")<<endl;
  

  // --- Output root file ---
  /*
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("sys_invmass[400,2000]");
  wUP->import(data,Rename("MC"));
  wUP->import(*gauss1);
  wUP->import(*expo);
  wUP->writeToFile("SignalMC1600-shapes-UnbinnedParam.root");
  */
}
