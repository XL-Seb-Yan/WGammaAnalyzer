#define function 7 //
#define data_type 3 //1:SB, 2:SBS(caled to)WB, 3:MCWB
#include <TMath.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <cstring>
#include <TCanvas.h>
void CHI2NLL()
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);
  TString funname;
  int nParam = 0;

#if data_type == 1
  TString type = "SB";
#elif data_type == 2
  TString type = "SBSWB";
#elif data_type == 3
  TString type = "MCWB";
#endif

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("sys_invmass","invmass",600,2500,""); //the name "sys_invmass" will be used by RooDataSet to import data

  // --- Import toy dataset ---
#if function == 1 && data_type == 1
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet1.SB.123456.root");
#elif function == 1 && data_type == 2
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet1.SBSWB.123456.root");
#elif function == 1 && data_type == 3
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet1.MCWB.123456.root");
#elif function == 2 && data_type == 1
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet2.SB.123456.root");
#elif function == 2 && data_type == 2
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet2.SBSWB.123456.root");
#elif function == 2 && data_type == 3
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet2.MCWB.123456.root");
#elif function == 3 && data_type == 1
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet3.SB.123456.root");
#elif function == 3 && data_type == 2
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet3.SBSWB.123456.root");
#elif function == 3 && data_type == 3
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.dijet3.MCWB.123456.root");
#elif function == 4 && data_type == 1
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.VVdijet1.SB.123456.root");
#elif function == 4 && data_type == 2
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.VVdijet1.SBSWB.123456.root");
#elif function == 4 && data_type == 3
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.VVdijet1.MCWB.123456.root");
#elif function == 5 && data_type == 1
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.VVdijet2.SB.123456.root");
#elif function == 5 && data_type == 2
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.VVdijet2.SBSWB.123456.root");
#elif function == 5 && data_type == 3
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.VVdijet2.MCWB.123456.root");
#elif function == 6 && data_type == 1
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.ATLAS1.SB.123456.root");
#elif function == 6 && data_type == 2
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.ATLAS1.SBSWB.123456.root");
#elif function == 6 && data_type == 3
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.ATLAS1.MCWB.123456.root");
#elif function == 7 && data_type == 1
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.ATLAS2.SB.123456.root");
#elif function == 7 && data_type == 2
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.ATLAS2.SBSWB.123456.root");
#elif function == 7 && data_type == 3
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/toys/higgsCombineTest.GenerateOnly.ATLAS2.MCWB.123456.root");
#endif

  TDirectory *d = (TDirectory*)file.Get("toys");

  // --- Histogram to store chi2 test results
  TH1F *chi2hist = new TH1F("chi2","",100,0,50);
  TH1F *NLLhist = new TH1F("chi2","",500,0,100000);
  TCanvas *c1 = new TCanvas("data","data",1200,900);
  
  for(int itoy=1; itoy<501; itoy++){

    RooArgList imarglist(*x);
    RooArgSet imargset(imarglist);
    TString toyname = "toy_"+std::to_string(itoy);
    RooDataSet *data = (RooDataSet*)d->Get(toyname);
    cout<<"Processing: ===================                "<<itoy<<"             ==================="<<endl;
    //RooDataSet data = ("Data sideband toys","Data sideband toys",imargset,Import(*import,"toy"));//import branches with names match the "variable name" (not variable) listed in imargset

    //--- background PDF ---
#if function == 1
    //-----------------------------dijet 1-----------------------------------
    TString fun_name = "dijet1";
    nParam = 1;
    RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
    RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0))",RooArgList(*x,*p0));
    //RooRealVar *nbkg = new RooRealVar("nbkg","nbkg",5000,0,100000,""); 
    //RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*nbkg));
#elif function == 2
    //-----------------------------dijet 2-----------------------------------
    TString fun_name = "dijet2";
    nParam = 2;
    RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
    RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1));
    //RooRealVar *nbkg = new RooRealVar("nbkg","nbkg",5000,0,100000,""); 
    //RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*nbkg));
#elif function == 3
    //-----------------------------dijet3-----------------------------------
    TString fun_name = "dijet3";
    nParam = 3;
    RooRealVar *p0 = new RooRealVar("p0","p_0",1,-1000,1000,""); //-10.5471 +- 3.06516
    RooRealVar *p1 = new RooRealVar("p1","p_1",-1,-1000,1000,""); //-0.799413 +- 0.560435
    RooRealVar *p2 = new RooRealVar("p2","p_2",0.5,-30,30,""); //-0.799413 +- 0.560435
    RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)+p2*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2));
#elif function == 4
    //-----------------------------VVdijet1-----------------------------------
    
    TString fun_name = "VVdijet1";
    nParam = 2;
    RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
    RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1))",RooArgList(*x,*p0,*p1));
#elif function == 5
    //-----------------------------VVdijet2-----------------------------------
    
    TString fun_name = "VVdijet2";
    nParam = 3;
    RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
    RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1,*p2));
#elif function == 6
    //-----------------------------ATLAS1-----------------------------------
    TString fun_name = "ATLAS1";
    nParam = 2;
    RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
    RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1))",RooArgList(*x,*p0,*p1));
#elif function == 7
    //-----------------------------ATLAS2-----------------------------------
    TString fun_name = "ATLAS2";
    nParam = 3;
    RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
    RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
    RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1,*p2));
#endif

    // --- Perform extended ML fit of composite PDF to toy data ---
    RooFitResult *r = model->fitTo(*data,Range(600,2500),RooFit::Minimizer("Minuit2"),Save());
  
    // --- Plot ---
    x->setBins(95); //fit is unbinned but chi2 is calculated by binning data with this value
    //gStyle->SetOptStat(111111);
    RooPlot *frame = x->frame(Title("Data Sideband"));
    data->plotOn(frame,RooFit::Name("data"));
    model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
    RooDataHist datah("dh","binned data",RooArgSet(*x),*data);
    RooChi2Var chi2 ("chi2", "chi2", *model,datah,DataError(RooAbsData::Poisson));
    RooAbsReal* nll = model->createNLL(*data);
    if(r->status() == 1 || r->status() == 0 || r->status() == 4 || r->status() == 3){
      //chi2hist->Fill(frame->chiSquare(fun_name,"data",nParam));
      chi2hist->Fill(chi2.getVal());
      NLLhist->Fill(nll->getVal());
    }
    chi2hist->SetTitle("#chi^{2} "+fun_name+" "+type);
    /*
      if(itoy%10 == 0){
      c1->Clear();
      c1->cd();
      frame->Draw();
      c1->SetLogy();
      c1->Print(type+" "+toyname+".png");
      }
    */
    funname = fun_name;
  }
  gStyle->SetOptStat(0);
  TCanvas *c0 = new TCanvas("01","01",1200,900);
  c0->cd();
  chi2hist->SetLineWidth(2);
  TAxis *xaxis = chi2hist->GetXaxis();
  TAxis *yaxis = chi2hist->GetYaxis();
  xaxis->SetTitle("#chi^{2}");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.2);
  TLatex latex;  
  TString Ntoystr = "# of valid toys: "+std::to_string((int)chi2hist->GetEntries());
  chi2hist->Draw();
  latex.SetTextSize(0.03);
  latex.DrawLatexNDC(.65,.75,Ntoystr);
  c0->Print("chi2"+funname+" "+type+".png");
  cout<<"mean of chi2 is: "<<chi2hist->GetMean(1)<<" std is: "<<chi2hist->GetStdDev(1)<<endl;
  cout<<"mean of NLL is: "<<NLLhist->GetMean(1)<<" std is: "<<NLLhist->GetStdDev(1)<<endl;
  cout<<funname<<" "<<type<<endl;
}
