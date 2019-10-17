#include <TMath.h>
#include <TLegend.h>
void make_signal_shapes(int seed=37)
{
  //gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37); 

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("sys_invmass","invmass",0,1500,""); //the name "sys_invmass" will be used by RooDataSet to import data
  x->setBins(150);

  //--- signal PDF ---
  double signalmass = 800;
  TString fun_name = "CB";
  RooRealVar* mean = new RooRealVar("CB_mean","mean",0,0,0.1,"");
  RooRealVar* CBsigma = new RooRealVar("CB_sigma","CBsigma",0,0,1,"");
  //RooRealVar* Gaussigma = new RooRealVar("Gaus_sigma","Gaussigma",0,0,1,"");
  RooRealVar* alpha = new RooRealVar("CB_alpha","alpha",5,1,10,"");
  RooRealVar* n = new RooRealVar("CB_n","n",1,1,2,"");
  RooRealVar* sqrtS = new RooRealVar("sqrtS","sqrtS",13000,"");
  sqrtS->setConstant(kTRUE);
  RooFormulaVar* xx = new RooFormulaVar("invmass/sqrtS","sys_invmass / sqrtS",RooArgSet(*x,*sqrtS));
  //RooRealVar* frac = new RooRealVar("Gaus_frac","frac",0.01,0,1); 
  //RooGaussian* Gaus = new RooGaussian("Gaus","Gaus",*xx,*mean,*Gaussigma);
  RooCBShape* model = new RooCBShape("CBShape","Cystal Ball Function",*xx,*mean,*CBsigma,*alpha,*n);
  //RooAddPdf* model = new RooAddPdf(fun_name,fun_name,RooArgList(*Gaus,*CB),*frac);
  
  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/Signal800N_Wwindow_full_finalcut.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(*x);
  RooArgSet imargset(imarglist);
  RooDataSet data("Signal800N","Signal800N",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  /*
  // --- Perform ML fit of composite PDF to data ---
  RooFitResult *r = model->fitTo(data,Range(600,1000),RooFit::Minimizer("Minuit2"),Save());

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = x->frame(Title("Signal800N"));
  data.plotOn(frame,RooFit::Name("data"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  data.plotOn(frame);
  TString chi2txt = "Chi2/Ndof: "+std::to_string(frame->chiSquare(fun_name,"data"));
  TString NLLtxt = "NLL: "+std::to_string(r->minNll());
  TLatex *chi2lax = new TLatex(1830,180,chi2txt);
  TLatex *NLLlax = new TLatex(1830,130,NLLtxt);
  chi2lax->SetTextSize(2.5);
  chi2lax->SetTextColor(kBlack);
  NLLlax->SetTextSize(2.5);
  NLLlax->SetTextColor(kBlack);
  frame->addObject(chi2lax);
  frame->addObject(NLLlax);

  TCanvas *c01 = new TCanvas("c01","c01",1200,900);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 20 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.1,10000);
  //frame->SetMaximum(600);
  //frame->SetMinimum(0);
  c01->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.6,0.7,0.8,0.78);
  l->AddEntry(frame->findObject(fun_name),"SignalMC fit "+fun_name,"l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
  c01->Print(fun_name+".png");

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
  c02->Print("pull"+fun_name+".png");
  */
  
  // --- Output root file ---
  
  RooWorkspace *w = new RooWorkspace("w","w");
  w->import(*x);
  w->import(data,Rename("data_sideband"));
  w->import(*model);
  w->writeToFile("siganl800N-shapes-Unbinned-unfitted"+fun_name+".root");
  
}
