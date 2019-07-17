#include <TMath.h>
#include <TLegend.h>
void make_sig_shapes(int seed=37)
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37); 
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  // --- Create workspace --- 
  RooWorkspace w("w","w");
  RooRealVar sys_invmass("sys_invmass","invmass",1400,2600,""); //the name "sysinvmass" will be used by RooDataSet to import data
  sys_invmass.setBins(48);

  /*
  //create roofit gaussian
  RooRealVar mean("mean","mean",1600,1360,1840) ;
  RooRealVar sigma("sigma","sigma",50,0,120) ;
  RooRealVar width("width", "width",5.,0.,120.) ;
  RooRealVar sig_yield("sig_yield","signal yield",1000,0,10000);
  RooVoigtian gauss("gauss","signal p.d.f.",sys_invmass,mean,width,sigma) ;
  */

  //create roofit gaussian
  RooRealVar mean("mean","mean",2000,1900,2100) ;
  RooRealVar sigma("sigma","sigma",100,0,200) ;
  RooRealVar sig_yield("sig_yield","signal yield",1000,0,10000);
  RooGaussian gauss("gauss","signal",sys_invmass,mean,sigma) ;

  //create roofit quadratic
  RooRealVar a1("a1","a1 coefficient of polynomial",-0.004,-0.01,-0.0001) ;
  RooRealVar bkg_yield("bkg_yield","background yield",1000,0,10000);
  RooExponential bkg("bkg","background", sys_invmass, a1) ;

  /*
  RooRealVar a1("a1","a1 coefficient of polynomial",-0.7, -2., 2.) ;
  RooRealVar a2("a2","a2 coefficient of polynomial",0.3, -2, 2.) ;
  RooRealVar a3("a3","a3 coefficient of polynomial",-0.03, -2, 2.) ;
  RooRealVar bkg_yield("bkg_yield","background yield",1000,0,10000);
  RooChebychev bkg("bkg","background", sys_invmass, RooArgList(a1,a2,a3)) ;
  */

  //create composite model model(x) = sig_yield*gauss(x) + bkg_yield*bkg(x)
  RooAddPdf model("model","model",RooArgList(gauss,bkg),RooArgList(sig_yield,bkg_yield));
  
  // --- Import unbinned dataset ---
  TFile file("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC2000_WGamma_select.root");
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
  RooFitResult *r = model.fitTo(data,Save());
  sys_invmass.setRange("window",1700,2300);
  RooAbsReal* bkg_N = gauss.createIntegral(sys_invmass,RooFit::NormSet(sys_invmass),Range("window"));

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = sys_invmass.frame();
  data.plotOn(frame);
  frame->SetTitle("WGamma M2000 MC");
  //ex_sum->plotOn(frame,LineColor(kRed));
  model.plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //gauss1->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  model.plotOn(frame,Components("bkg"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  model.plotOn(frame,VisualizeError(*r,3),FillColor(kOrange),RooFit::Name("err"));
  //gauss1->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  model.plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //ex_bern->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  model.plotOn(frame,Components("bkg"),LineStyle(kDashed),RooFit::Name("bkgfun"));
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
  cout<<bkg_N->getValV()<<" "<<bkg_N->getValV()*sig_yield.getValV()<<endl;
  
  // --- Output root file ---
  /*
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("sys_invmass[1120,2080]");
  wUP->import(data,Rename("MC"));
  wUP->import(gauss);
  wUP->import(bkg);
  wUP->writeToFile("SignalMC1600-shapes-UnbinnedParam.root");
  */
}
