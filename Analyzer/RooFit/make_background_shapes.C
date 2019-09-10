#include <TMath.h>
#include <TLegend.h>
void make_background_shapes(int seed=37)
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37); 
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  // --- Create workspace --- 
  RooWorkspace w("w","w");
  RooRealVar sys_invmass("sys_invmass","invmass",700,2500,""); //the name "sysinvmass" will be used by RooDataSet to import data
  sys_invmass.setBins(90);

  //--- background PDF ---
  RooRealVar *p0 = new RooRealVar("p0","p_0",1,0,3,"");
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1,1,"");
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1,1,"");
  RooRealVar *Nbkg = new RooRealVar("pol2_Nbkg","N_{bkg}",500,0,5000,"events");
  RooRealVar *x = new RooRealVar("x","x",700,2500);
  RooGenericPdf *dijet = new RooGenericPdf("dijet","p0*x^(p1+p2*log(x))",RooArgList(*x,*p0,*p1,*p2));
  //RooExponential* expo = new RooExponential("expo","",sys_invmass,*p0);
  RooAddPdf* ex_sum = new RooAddPdf("ex_sum","",RooArgList(*dijet),RooArgList(*Nbkg));
  w.import(*ex_sum);
 
  
  // --- Import unbinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/SinglePhoton2017_sideband.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(sys_invmass);
  RooArgSet imargset(imarglist);
  RooDataSet data("Data sideband","Data sideband",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

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
  ex_sum->plotOn(frame,Components("dijet"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  ex_sum->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  ex_sum->plotOn(frame,Components("dijet"),LineStyle(kDashed),RooFit::Name("bkgfun"));
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
