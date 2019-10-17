#include <TMath.h>
#include <TLegend.h>
void make_sideband_shapes(int seed=37)
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37); 

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("sys_invmass","invmass",600,3000,""); //the name "sys_invmass" will be used by RooDataSet to import data

  //--- background PDF ---
  
  //-----------------------------dijet 1-----------------------------------
  /*
  TString fun_name = "dijet1";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0))",RooArgList(*x,*p0));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------dijet 2-----------------------------------
  /*
  TString fun_name = "dijet2";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------dijet3-----------------------------------
  /*
  TString fun_name = "dijet3";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)+p2*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------ATLAS1-----------------------------------
  
  TString fun_name = "ATLAS1";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1))",RooArgList(*x,*p0,*p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  
  //-----------------------------ATLAS2-----------------------------------
  /*
  TString fun_name = "ATLAS2";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1,*p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------ATLAS3-----------------------------------
  /*
  TString fun_name = "ATLAS3";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)+p3*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2,*p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------VVdijet1-----------------------------------
  /*
  TString fun_name = "VVdijet1";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1))",RooArgList(*x,*p0,*p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------VVdijet2-----------------------------------
  /*
  TString fun_name = "VVdijet2";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1,*p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------VVdijet3-----------------------------------
  /*
  TString fun_name = "VVdijet3";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)+p3*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2,*p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------expow1-----------------------------------
  /*
  TString fun_name = "expow1";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1)",RooArgList(*x,*p0,*p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------expow2-----------------------------------
  /*
  TString fun_name = "expow2";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1+p2*(sys_invmass/13000))",RooArgList(*x,*p0,*p1,*p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
  //-----------------------------expow3-----------------------------------
  /*
  TString fun_name = "expow3";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1+p2*(sys_invmass/13000)+p3*pow(sys_invmass/13000,2))",RooArgList(*x,*p0,*p1,*p2,*p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  */
 
  
  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/Singlephoton2017_Wsideband_full_finalcut.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(*x);
  RooArgSet imargset(imarglist);
  RooDataSet data("Data sideband","Data sideband",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  // --- Perform extended ML fit of composite PDF to toy data ---
  RooFitResult *r = model->fitTo(data,Range(600,3000),RooFit::Minimizer("Minuit2"),Save());

  
  // --- Plot ---
  x->setBins(120); //fit is unbinned but chi2 is calculated by binning data with this value
  gStyle->SetOptStat(111111);
  RooPlot *frame = x->frame(Title("Data Sideband"));
  data.plotOn(frame,RooFit::Name("data"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  data.plotOn(frame);
  RooDataHist datah("dh","binned data",RooArgSet(*x),data);
  RooChi2Var chi2 ("chi2", "chi2", *model,datah,DataError(RooAbsData::Poisson));
  TString chi2txt = "Chi2/Ndof: "+std::to_string(frame->chiSquare(fun_name,"data",3));
  cout<<"chi2 is: "<<chi2.getVal()<<endl;
  TString NLLtxt = "NLL: "+std::to_string(r->minNll());
  TLatex *chi2lax = new TLatex(1830,180,chi2txt);
  TLatex *NLLlax = new TLatex(1830,130,NLLtxt);
  chi2lax->SetTextSize(0.025);
  chi2lax->SetTextColor(kBlack);
  NLLlax->SetTextSize(0.025);
  NLLlax->SetTextColor(kBlack);
  frame->addObject(chi2lax);
  frame->addObject(NLLlax);

  // access residuals
  double SSR = 0;
  RooCurve* curve = (RooCurve*)frame->findObject(fun_name,RooCurve::Class());
  RooHist* hist = (RooHist*)frame->findObject("data",RooHist::Class());
  Double_t xstart,xstop,y;
  curve->GetPoint(0,xstart,y);
  curve->GetPoint(curve->GetN()-1,xstop,y);
  for(Int_t i=0 ; i<hist->GetN() ; i++) {
    Double_t x,point;
    hist->GetPoint(i,x,point);
    // Only calculate pull for bins inside curve range
    if (x<xstart || x>xstop) continue;
    Double_t yy;
    yy = point - curve->interpolate(x);
    SSR += pow(yy,2);
  }
  cout<<"SSR is: "<<SSR<<" NdoF is: "<<hist->GetN()<<endl;
  TString SSRtxt = "SSR: "+std::to_string(SSR);
  TString NdoFtxt = "NLL: "+std::to_string(r->minNll());
  TLatex *SSRlax = new TLatex(1830,80,SSRtxt);
  SSRlax->SetTextSize(0.025);
  SSRlax->SetTextColor(kBlack);
  frame->addObject(SSRlax);
  
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
  l->AddEntry(frame->findObject(fun_name),"Data SB fit "+fun_name,"l");
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
  yaxis->SetTitle("#frac{data-fit}{#sigma_{stat}}");
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
  
  /*
  // --- Output root file ---
  RooWorkspace *w = new RooWorkspace("w","w");
  w->import(*x);
  w->import(data,Rename("data_sideband"));
  w->import(*model);
  w->writeToFile("sideband-shapes-Unbinned-"+fun_name+".root");

  // --- Perform extended ML fit of composite PDF to toy data ---
  RooFitResult *ex_r = ex_model->fitTo(data,Range(600,3000),RooFit::Minimizer("Minuit2"),Extended(true),Save());
  cout<<"Normalization is: "<<bkg_norm->getVal()<<endl;
  */
}
