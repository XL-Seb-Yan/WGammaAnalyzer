#include <TMath.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <cstring>
#include <TCanvas.h>
void NLL()
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);
  TString funname;
  int nParam = 0;

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("sys_invmass","invmass",600,3000,""); //the name "sys_invmass" will be used by RooDataSet to import data

  // --- Import toy dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/higgsCombineTest.GenerateOnly.dijet2.123456.root");
  TDirectory *d = (TDirectory*)file.Get("toys");

  // --- Histogram to store chi2 test results
  TH1F *nllhist = new TH1F("NLL","",100,4000,9000);
  TCanvas *c1 = new TCanvas("data","data",1200,900);
  
  for(int itoy=1; itoy<501; itoy++){

    RooArgList imarglist(*x);
    RooArgSet imargset(imarglist);
    TString toyname = "toy_"+std::to_string(itoy);
    RooDataSet *data = (RooDataSet*)d->Get(toyname);
    cout<<"Processing: "<<itoy<<"             ==================="<<endl;
    //RooDataSet data = ("Data sideband toys","Data sideband toys",imargset,Import(*import,"toy"));//import branches with names match the "variable name" (not variable) listed in imargset

    //--- background PDF ---
    //-----------------------------dijet 1-----------------------------------
    /*
      TString fun_name = "dijet1";
      nParam = 1;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0))",RooArgList(*x,*p0));
      //RooRealVar *nbkg = new RooRealVar("nbkg","nbkg",5000,0,100000,""); 
      //RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*nbkg));
      */
    //-----------------------------dijet 2-----------------------------------
    
      TString fun_name = "dijet2";
      nParam = 2;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1));
      //RooRealVar *nbkg = new RooRealVar("nbkg","nbkg",5000,0,100000,""); 
      //RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*nbkg));
      
    //-----------------------------dijet3-----------------------------------
    /*
      TString fun_name = "dijet3";
      nParam = 3;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)+p2*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2));
    */
    //-----------------------------ATLAS1-----------------------------------
    /*
      TString fun_name = "ATLAS1";
      nParam = 2;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1))",RooArgList(*x,*p0,*p1));
    */
    //-----------------------------ATLAS2-----------------------------------
    /*
      TString fun_name = "ATLAS2";
      nParam = 3;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1,*p2));
    */
    //-----------------------------ATLAS3-----------------------------------
    /*
      TString fun_name = "ATLAS3";
      nParam = 4;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)+p3*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2,*p3));
    */
    //-----------------------------VVdijet1-----------------------------------
    /*
      TString fun_name = "VVdijet1";
      nParam = 2;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1))",RooArgList(*x,*p0,*p1));
    */
    //-----------------------------VVdijet2-----------------------------------
    /*
      TString fun_name = "VVdijet2";
      nParam = 3;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1,*p2));
    */
    //-----------------------------VVdijet3-----------------------------------
    /*
      TString fun_name = "VVdijet3";
      nParam = 4;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)+p3*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2,*p3));
    */
    //-----------------------------expow1-----------------------------------
    /*
      TString fun_name = "expow1";
      nParam = 2;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1)",RooArgList(*x,*p0,*p1));
    */
    //-----------------------------expow2-----------------------------------
    /*
    TString fun_name = "expow2";
    nParam = 3;
    RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
    RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1+p2*(sys_invmass/13000))",RooArgList(*x,*p0,*p1,*p2));
    */
    //-----------------------------expow3-----------------------------------
    /*
      TString fun_name = "expow3";
      nParam = 4;
      RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
      RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
      RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1+p2*(sys_invmass/13000)+p3*pow(sys_invmass/13000,2))",RooArgList(*x,*p0,*p1,*p2,*p3));
    */

    // --- Perform extended ML fit of composite PDF to toy data ---
    RooFitResult *r = model->fitTo(*data,Range(600,3000),RooFit::Minimizer("Minuit2"),Save());
  
    // --- Plot ---
    x->setBins(120); //fit is unbinned but chi2 is calculated by binning data with this value
    //gStyle->SetOptStat(111111);
    RooPlot *frame = x->frame(Title("Data Sideband"));
    data->plotOn(frame,RooFit::Name("data"));
    model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
    if(r->status() == 1 || r->status() == 0 || r->status() == 4 || r->status() == 3){
      //nllhist->Fill(frame->chiSquare(fun_name,"data",nParam));
      RooAbsReal* nll = model->createNLL(*data);
      nllhist->Fill(nll->getVal());
    }
    nllhist->SetTitle("NLL"+fun_name);
    /*
      if(itoy%10 == 0){
      c1->Clear();
      c1->cd();
      frame->Draw();
      c1->SetLogy();
      c1->Print("data"+toyname+".png");
      }
    */
    funname = fun_name;
  }
  gStyle->SetOptStat(0);
  TCanvas *c0 = new TCanvas("01","01",1200,900);
  c0->cd();
  nllhist->SetLineWidth(2);
  TAxis *xaxis = nllhist->GetXaxis();
  TAxis *yaxis = nllhist->GetYaxis();
  xaxis->SetTitle("nll");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.2);
  TLatex latex;  
  TString Ntoystr = "# of valid toys: "+std::to_string((int)nllhist->GetEntries());
  nllhist->Draw();
  latex.SetTextSize(0.04);
  latex.DrawLatexNDC(.65,.75,Ntoystr);
  c0->Print("nll"+funname+".png");
  cout<<"mean of chi2 is: "<<nllhist->GetMean(1)<<" std is: "<<nllhist->GetStdDev(1)<<endl;
}
