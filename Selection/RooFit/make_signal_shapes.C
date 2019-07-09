#include <TMath.h>
#include <TLegend.h>
void make_signal_shapes(int seed=37)
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37); 
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  // --- Create workspace --- 
  RooWorkspace w("w","w");
  RooRealVar sys_invmass("sys_invmass","invmass",400,2000,""); //the name "sysinvmass" will be used by RooDataSet to import data
  sys_invmass.setBins(25);

  //--- signal PDF ---
  RooRealVar* mu1 = new RooRealVar("Gaus_mu1","#mu1",1600,1500,1700,"");
  //RooRealVar* mu2 = new RooRealVar("DG_mu2","#mu2",1.96847,1.9,2.1,"");
  RooRealVar* sigma1 = new RooRealVar("Gaus_sigma1","#sigma_{1}",80,0,150,"");
  //RooRealVar* sigma2 = new RooRealVar("DG_sigma2","#sigma_{2}",0.03,0.02,0.08,"");
  mu1->setConstant(kFALSE);
  //mu2->setConstant(kFALSE);
  sigma1->setConstant(kFALSE);
  //sigma2->setConstant(kFALSE);
  RooRealVar* frac = new RooRealVar("Gaus_frac","frac",0.7,0,1); 
  RooRealVar* Ns = new RooRealVar("Gaus_Ns","N_{s}",1000,0,10000,"events");
  Ns->setConstant(kFALSE);
  RooGaussian* gauss1 = new RooGaussian("Gaus1","",sys_invmass,*mu1,*sigma1);
  //RooGaussian* gauss2 = new RooGaussian("G2","",sysinvmass,*mu2,*sigma2);
  //RooAddPdf* doublegauss = new RooAddPdf("SummedG1G2","",RooArgList(*gauss1,*gauss2),*frac);
  //RooAddPdf* ex_doublegauss = new RooAddPdf("DG_ex","extDgauss",RooArgList(*doublegauss),RooArgList(*Ns));
  RooAddPdf* ex_gauss = new RooAddPdf("Gaus_ex","",RooArgList(*gauss1),RooArgList(*Ns));
  w.import(*ex_gauss);

  //--- background PDF ---
  RooRealVar *pC = new RooRealVar("pol2_pC","C",1,0,100,"");
  pC->setConstant(kFALSE);
  RooRealVar *p0 = new RooRealVar("pol2_p0","p_0",-0.1,-1,0,"");
  p0->setConstant(kFALSE);
  RooRealVar *Nbkg   = new RooRealVar("pol2_Nbkg","N_{bkg}",500,0,2000,"events");
  Nbkg->setConstant(kFALSE);
  RooExponential* expo = new RooExponential("expo","",sys_invmass,*p0);
  RooAddPdf* ex_expo = new RooAddPdf("expo_ex","",RooArgList(*expo),RooArgList(*Nbkg));
  w.import(*ex_expo);

  //--- combine PDF ---
  RooAddPdf* ex_sum = new RooAddPdf("SUM_ex","",RooArgList(*expo,*gauss1),RooArgList(*Nbkg,*Ns));
  w.import(*ex_sum);
  
  // --- Import unbinned dataset ---
  TFile file("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/MC_WGamma_select.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(sys_invmass);
  RooArgSet imargset(imarglist);
  RooDataSet data("Signal MC","Signal MC",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  // --- Perform extended ML fit of composite PDF to toy data ---
  //w.pdf("model_s")->fitTo(data);
  //sys_invmass.setRange("SB1",1.47682,1.67682);
  //sysinvmass.setRange("SB2",1.87682,2.07682);
  //RooFitResult *r = ex_bern->fitTo(data,Range("SB1,SB2"),Save());
  //RooFitResult *r = ex_doublegauss->fitTo(data,Save());
  //RooFitResult *r = gauss1->fitTo(data,Save());
  RooFitResult *r = ex_sum->fitTo(data,Save());

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = sys_invmass.frame();
  data.plotOn(frame);
  frame->SetTitle("WGamma M1600 MC");
  //ex_sum->plotOn(frame,LineColor(kRed));
  ex_sum->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //gauss1->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  ex_sum->plotOn(frame,Components("expo"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  ex_sum->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  //gauss1->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  ex_sum->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //ex_bern->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  ex_sum->plotOn(frame,Components("expo"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  data.plotOn(frame);
  //w.pdf("model_s")->plotOn(frame,Components("Bern"),LineStyle(kDashed));

  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{W#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events");
  yaxis->SetTitleOffset(1.2);
  //frame->SetMaximum(600);
  //frame->SetMinimum(0);
  //c1->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.6,0.7,0.8,0.78);
  l->AddEntry(frame->findObject("sigfun"),"Signal Fit","l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err"),"Fit Error 1 #sigma","f");
  //l->Draw("same");
  c1->Print("data.png");
  
  // --- Output root file ---
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("sys_invmass[400,2000]");
  wUP->import(data,Rename("MC"));
  wUP->import(*gauss1);
  wUP->import(*expo);
  wUP->writeToFile("SignalMC1600-shapes-UnbinnedParam.root");
}
