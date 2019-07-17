#include <TMath.h>
#include <TLegend.h>
void make_bkg_shapes(int seed=37)
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37); 
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  // --- Create workspace --- 
  RooWorkspace w("w","w");
  RooRealVar sys_invmass("sys_invmass","invmass",1120,2080,""); //the name "sysinvmass" will be used by RooDataSet to import data
  sys_invmass.setBins(80);

  //create bkg function
  RooRealVar a1("a1","a1 coefficient of polynomial",-0.004,-0.005,-0.003) ;
  RooRealVar bkg_yield("bkg_yield","background yield",10000,0,100000);
  RooExponential bkg("bkg","background", sys_invmass, a1) ;
  RooAddPdf model("model","model",RooArgList(bkg),RooArgList(bkg_yield)) ;
   
  // --- Import unbinned dataset ---
  TFile file("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017_WGamma_select.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(sys_invmass);
  RooArgSet imargset(imarglist);
  RooDataSet data("Signal MC","Signal MC",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  // --- Perform fit of composite PDF to toy data ---
  sys_invmass.setRange("SB1",1120,1360);
  sys_invmass.setRange("SB2",1840,2080);
  //RooFitResult *r = ex_bern->fitTo(data,Range("SB1,SB2"),Save());
  RooFitResult *r = model.fitTo(data,Range("SB1,SB2"),Save());
  sys_invmass.setRange("window",1360,1840);
  RooAbsReal* bkg_N = model.createIntegral(sys_invmass,RooFit::NormSet(sys_invmass),Range("window"));

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = sys_invmass.frame();
  data.plotOn(frame);
  frame->SetTitle("2017 SinglePhoton");
  model.plotOn(frame,LineColor(kRed),RooFit::Name("bkgfun"));
  model.plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  model.plotOn(frame,LineColor(kRed),RooFit::Name("bkgfun"));
  data.plotOn(frame);
  
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
  l->AddEntry(frame->findObject("bkgfun"),"Data Fit","l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err"),"Fit Error 1 #sigma","f");
  //l->Draw("same");
  c1->Print("data.png");
  
  // --- Output root file ---
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("sys_invmass[1120,2080]");
  wUP->import(data,Rename("MC"));
  //wUP->import(gauss);
  //wUP->import(bkg);
  wUP->writeToFile("SignalMC1600-shapes-UnbinnedParam.root");
  std::cout<<"==============================================="<<std::endl;
  std::cout<<"==============================================="<<std::endl;
  std::cout<<"==============================================="<<std::endl;
  std::cout<<"signal: "<<bkg_N->getValV()<<std::endl;
  std::cout<<"==============================================="<<std::endl;
  std::cout<<"==============================================="<<std::endl;
  std::cout<<"==============================================="<<std::endl;
}
