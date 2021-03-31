#define fun_type 3
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
#include <RooBukinPdf.h>
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

void make_SR_binned_shapes(int seed=37)
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37); 
  
  gROOT->SetBatch(1);
  lumi_13TeV = "137.19 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Preliminary";
  lumiTextSize = 0.4;
  cmsTextSize = 0.5;
  int iPeriod = 12;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.035);
  gStyle->SetBarWidth(1.03);
  gStyle->SetHistLineWidth(2);

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("m","m",600,5000,""); //the name "m" will be used by RooDataSet to import data, normalization range is 600-3500 but plot range can be defined to like 600-3000

  //--- background PDF ---
#if fun_type == 1
  //-----------------------------dijet 1-----------------------------------
  TString fun_name = "dijet-1";
  RooRealVar *dijet1_p0 = new RooRealVar("dijet1_p0","dijet1_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(m/13000,dijet1_p0))",RooArgList(*x,*dijet1_p0));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 2
//-----------------------------dijet 2-----------------------------------
  TString fun_name = "dijet-2";
  RooRealVar *dijet2_p0 = new RooRealVar("dijet2_p0","dijet2_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *dijet2_p1 = new RooRealVar("dijet2_p1","dijet2_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(m/13000,dijet2_p0+dijet2_p1*log(m/13000)))",RooArgList(*x,*dijet2_p0,*dijet2_p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 3
//-----------------------------dijet 3-----------------------------------
  TString fun_name = "dijet-3";
  RooRealVar *dijet3_p0 = new RooRealVar("dijet3_p0","dijet3_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *dijet3_p1 = new RooRealVar("dijet3_p1","dijet3_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *dijet3_p2 = new RooRealVar("dijet3_p2","dijet3_p_2",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(m/13000,dijet3_p0+dijet3_p1*log(m/13000)+dijet3_p2*pow(log(m/13000),2)))",RooArgList(*x,*dijet3_p0,*dijet3_p1,*dijet3_p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 4
  //-----------------------------dijet 4 -----------------------------------
  TString fun_name = "dijet-4";
  RooRealVar *dijet4_p0 = new RooRealVar("dijet4_p0","dijet4_p_0",0,-20,20,""); //-10.5471 +- 3.06516
  RooRealVar *dijet4_p1 = new RooRealVar("dijet4_p1","dijet4_p_1",0,-20,20,""); //-0.799413 +- 0.560435
  RooRealVar *dijet4_p2 = new RooRealVar("dijet4_p2","dijet4_p_2",0,-20,20,""); //-0.799413 +- 0.560435
  RooRealVar *dijet4_p3 = new RooRealVar("dijet4_p3","dijet4_p_3",0,-20,20,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(m/13000,dijet4_p0+dijet4_p1*log(m/13000)+dijet4_p2*pow(log(m/13000),2)+dijet4_p3*pow(log(m/13000),3)))",RooArgList(*x,*dijet4_p0,*dijet4_p1,*dijet4_p2,*dijet4_p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 5
//-----------------------------VVdijet1-----------------------------------
  TString fun_name = "VVdijet-1";
  RooRealVar *VVdijet1_p0 = new RooRealVar("VVdijet1_p0","VVdijet1_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *VVdijet1_p1 = new RooRealVar("VVdijet1_p1","VVdijet1_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-m/13000,VVdijet1_p0)/pow(m/13000,VVdijet1_p1))",RooArgList(*x,*VVdijet1_p0,*VVdijet1_p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 6
  //-----------------------------VVdijet2-----------------------------------
  TString fun_name = "VVdijet-2";
  RooRealVar *VVdijet2_p0 = new RooRealVar("VVdijet2_p0","VVdijet2_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *VVdijet2_p1 = new RooRealVar("VVdijet2_p1","VVdijet2_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet2_p2 = new RooRealVar("VVdijet2_p2","VVdijet2_p_2",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-m/13000,VVdijet2_p0)/pow(m/13000,VVdijet2_p1+VVdijet2_p2*log(m/13000)))",RooArgList(*x,*VVdijet2_p0,*VVdijet2_p1,*VVdijet2_p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 7
//-----------------------------VVdijet3-----------------------------------
  TString fun_name = "VVdijet-3";
  RooRealVar *VVdijet3_p0 = new RooRealVar("VVdijet3_p0","VVdijet3_p_0",-17,-30,0,""); //-10.5471 +- 3.06516
  RooRealVar *VVdijet3_p1 = new RooRealVar("VVdijet3_p1","VVdijet3_p_1",16,0,25,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet3_p2 = new RooRealVar("VVdijet3_p2","VVdijet3_p_2",1,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet3_p3 = new RooRealVar("VVdijet3_p3","VVdijet3_p_3",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-m/13000,VVdijet3_p0)/pow(m/13000,VVdijet3_p1+VVdijet3_p2*log(m/13000)+VVdijet3_p3*pow(log(m/13000),2)))",RooArgList(*x,*VVdijet3_p0,*VVdijet3_p1,*VVdijet3_p2,*VVdijet3_p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 8
//-----------------------------ATLAS1-----------------------------------
  TString fun_name = "ATLAS-1";
  RooRealVar *ATLAS1_p0 = new RooRealVar("ATLAS1_p0","ATLAS1_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *ATLAS1_p1 = new RooRealVar("ATLAS1_p1","ATLAS1_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(m/13000,1/3),ATLAS1_p0)/pow(m/13000,ATLAS1_p1))",RooArgList(*x,*ATLAS1_p0,*ATLAS1_p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 9
  //-----------------------------ATLAS2-----------------------------------
  TString fun_name = "ATLAS-2";
  RooRealVar *ATLAS2_p0 = new RooRealVar("ATLAS2_p0","ATLAS2_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *ATLAS2_p1 = new RooRealVar("ATLAS2_p1","ATLAS2_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *ATLAS2_p2 = new RooRealVar("ATLAS2_p2","ATLAS2_p_2",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(m/13000,1/3),ATLAS2_p0)/pow(m/13000,ATLAS2_p1+ATLAS2_p2*log(m/13000)))",RooArgList(*x,*ATLAS2_p0,*ATLAS2_p1,*ATLAS2_p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 10
  //-----------------------------ATLAS3-----------------------------------
  TString fun_name = "ATLAS-3";
  RooRealVar *ATLAS3_p0 = new RooRealVar("ATLAS3_p0","ATLAS3_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *ATLAS3_p1 = new RooRealVar("ATLAS3_p1","ATLAS3_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *ATLAS3_p2 = new RooRealVar("ATLAS3_p2","ATLAS3_p_2",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *ATLAS3_p3 = new RooRealVar("ATLAS3_p3","ATLAS3_p_3",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(m/13000,1/3),ATLAS3_p0)/pow(m/13000,ATLAS3_p1+ATLAS3_p2*log(m/13000)+ATLAS3_p3*pow(log(m/13000),2)))",RooArgList(*x,*ATLAS3_p0,*ATLAS3_p1,*ATLAS3_p2,*ATLAS3_p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 11
  //-----------------------------dijet5-----------------------------------
  TString fun_name = "dijet-5";
  RooRealVar *dijet5_p0 = new RooRealVar("dijet5_p0","dijet5_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *dijet5_p1 = new RooRealVar("dijet5_p1","dijet5_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *dijet5_p2 = new RooRealVar("dijet5_p2","dijet5_p_2",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *dijet5_p3 = new RooRealVar("dijet5_p3","dijet5_p_3",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *dijet5_p4 = new RooRealVar("dijet5_p4","dijet5_p_4",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(m/13000,dijet5_p0+dijet5_p1*log(m/13000)+dijet5_p2*pow(log(m/13000),2)+dijet5_p3*pow(log(m/13000),3)+dijet5_p4*pow(log(m/13000),4)))",RooArgList(*x,*dijet5_p0,*dijet5_p1,*dijet5_p2,*dijet5_p3,*dijet5_p4));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 12
  //-----------------------------VVdijet4-----------------------------------
  TString fun_name = "VVdijet-4";
  RooRealVar *VVdijet4_p0 = new RooRealVar("VVdijet4_p0","VVdijet4_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *VVdijet4_p1 = new RooRealVar("VVdijet4_p1","VVdijet4_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet4_p2 = new RooRealVar("VVdijet4_p2","VVdijet4_p_2",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet4_p3 = new RooRealVar("VVdijet4_p3","VVdijet4_p_3",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet4_p4 = new RooRealVar("VVdijet4_p4","VVdijet4_p_4",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-m/13000,VVdijet4_p0)/pow(m/13000,VVdijet4_p1+VVdijet4_p2*log(m/13000)+VVdijet4_p3*pow(log(m/13000),2)+VVdijet4_p4*pow(log(m/13000),3)))",RooArgList(*x,*VVdijet4_p0,*VVdijet4_p1,*VVdijet4_p2,*VVdijet4_p3,*VVdijet4_p4));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 13
  //-----------------------------ATLAS4-----------------------------------
  TString fun_name = "ATLAS-4";
  RooRealVar *ATLAS4_p0 = new RooRealVar("ATLAS4_p0","ATLAS4_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *ATLAS4_p1 = new RooRealVar("ATLAS4_p1","ATLAS4_p_1",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *ATLAS4_p2 = new RooRealVar("ATLAS4_p2","ATLAS4_p_2",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *ATLAS4_p3 = new RooRealVar("ATLAS4_p3","ATLAS4_p_3",0,-100,100,""); //-0.799413 +- 0.560435
  RooRealVar *ATLAS4_p4 = new RooRealVar("ATLAS4_p4","ATLAS4_p_4",0,-100,100,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(m/13000,1/3),ATLAS4_p0)/pow(m/13000,ATLAS4_p1+ATLAS4_p2*log(m/13000)+ATLAS4_p3*pow(log(m/13000),2)+ATLAS4_p4*pow(log(m/13000),3)))",RooArgList(*x,*ATLAS4_p0,*ATLAS4_p1,*ATLAS4_p2,*ATLAS4_p3,*ATLAS4_p4));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#endif

  // --- Import Binned dataset ---
  
  //float xbinlow[53] = {600,650,700,750,800,850,900,950,1000,1050,1100,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3200,3400,3600,3800,400,4200,4400,4600,4800,5000,5400,5800,6200,6600,7000,7400,7800};
  // ------------------------Photons------------------------
  //TH1F* DATAh = new TH1F("DATAh","DATAh",52,&xbinlow[0]);
  TH1F* DATAh = new TH1F("DATAh","DATAh",110,600,5000);
  float s_mass;
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/fullcut/Run2Data_postproc_WGammaRun2_SR_sigrange_fullcut_jmcorr_May22.root");
  TTree* tree = (TTree*)file.Get("Events");
  tree->SetBranchAddress("m", &s_mass);
  for (int ievt = 0; ievt<tree->GetEntries();ievt++) {
    tree->GetEntry(ievt);
    DATAh->Fill(s_mass, 1);
    //DATAh->Fill(s_mass, 2.4826);
  }
  RooDataHist datah("Data SR binned","Data SR binned",RooArgSet(*x),DATAh);
  cout<<"number of weighted entries: "<<datah.sum(false)<<endl;
  
  // --- Import unBinned dataset ---
  RooDataSet data("Data SR","Data SR",RooArgSet(*x),Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset
  
  // --- Perform extended ML fit of composite PDF to toy data ---
  RooFitResult *r = model->fitTo(datah,RooFit::Minimizer("Minuit2"),Save()); //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)
    
  // --- plot for chi2 calculation and visualization ---
  RooPlot *frame = x->frame();

  frame->SetTitle("Data Sideband");
  datah.plotOn(frame,RooFit::Name("datah"),DataError(RooAbsData::Poisson),MarkerStyle(20)); //for unweighted data
  // Change attributes of last added plot elements
  frame->getAttMarker()->SetMarkerSize(2);
  model->plotOn(frame,LineWidth(2),RooFit::Name(fun_name));
  model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  model->plotOn(frame,LineWidth(2),RooFit::Name(fun_name));

    datah.plotOn(frame,RooFit::Name("datah"),DataError(RooAbsData::Poisson)); //for unweighted data

  frame->Print("V");
  RooAbsReal* nll = NULL;
    nll = model->createNLL(datah);
  
  cout<<"NLL is: "<<nll->getVal()<<endl;

  cout<<datah<<endl;

  cout<<"data bins: "<<datah.numEntries()<<endl;
  int n_0 = 0;
  for(int i=0; i<datah.numEntries(); i++){
    datah.get(i) ;
    if(datah.weight() == 0)
      n_0++;
  }
  cout<<"Number of empty bins: "<<n_0<<endl;
  
  RooArgSet *flparams = model->getParameters(*x);
  int nfloparam = (flparams->selectByAttrib("Constant",kFALSE))->getSize();
  cout<<"# of floating params: "<<nfloparam<<endl;

  // Chi2 problem solved:
  // Unweighted data:
  // RooChi2Var uses RooAbsData::Possion to obtain Possion interval (which can treat 0 entry bins), only skip bins with both PDF and data = 0
  // chi2square() obtains chi2 by skipping 0 entry bins, divide by numbr of non-0 entry bins
  // Weighted data:
  // RooChi2Var uses RooAbsData::SumW2, fail when arrive 0-entry bins
  // chi2square() obtains chi2 by skipping 0 entry bins, divide by numbr of non-0 entry bins
  RooChi2Var chi2 ("chi2", "chi2", *model,datah,DataError(RooAbsData::Poisson));//Default: SumW2
  cout<<"Chi2 from RooRealVar with possion is: "<<chi2.getVal()<<endl;
  cout<<"Reduced chi2(from RooFrame)/ndof is: "<<frame->chiSquare(fun_name,"datah", nfloparam)<<endl;
  cout<<"Chi2(from RooFrame)/ndof is: "<<frame->chiSquare(fun_name,"datah")<<endl;
  
  /*
  TString NLLtxt = "NLL: "+to_str_trim(nll->getVal());
  TLatex *NLLlax = new TLatex(0.5,0.6,NLLtxt);
  NLLlax->SetNDC();
  NLLlax->SetTextSize(0.028);
  NLLlax->SetTextColor(kBlack);
  frame->addObject(NLLlax);

  // access residuals
  double SSR = 0;
  RooCurve* curve = (RooCurve*)frame->findObject(fun_name,RooCurve::Class());
  RooHist* hist = (RooHist*)frame->findObject("datah",RooHist::Class());
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
    cout<<x<<" "<<point<<" "<<curve->interpolate(x)<<" "<<yy<<endl;
  }
  cout<<"SSR is: "<<SSR<<" NdoF is: "<<hist->GetN()<<endl;
  TString SSRtxt = "SSR: "+to_str_trim(SSR);
  TLatex *SSRlax = new TLatex(0.5,0.57,SSRtxt);
  SSRlax->SetNDC();
  SSRlax->SetTextSize(0.028);
  SSRlax->SetTextColor(kBlack);
  frame->addObject(SSRlax);
  */
  
  // write Chi2
  TString chi2txt = "#chi^{2}/ndof: " + to_str_trim(frame->chiSquare(fun_name,"datah") * (datah.numEntries() - n_0)) + "/" + std::to_string(datah.numEntries() - n_0 - nfloparam) + "=" + to_str_trim(frame->chiSquare(fun_name,"datah", nfloparam));
  TLatex *chi2lax = new TLatex(0.5,0.55,chi2txt);
  chi2lax->SetNDC();
  chi2lax->SetTextSize(0.028);
  chi2lax->SetTextColor(kBlack);
  frame->addObject(chi2lax);
  
   // --- Perform extended ML fit of composite PDF to toy data ---
  RooFitResult *ex_r = NULL;
  // if(isNorm)
    // ex_r = ex_model->fitTo(data_norm,Range(600,5000),RooFit::Minimizer("Minuit2"),Extended(true),SumW2Error(false),Save());
  // else
    ex_r = ex_model->fitTo(datah,Range(600,5000),RooFit::Minimizer("Minuit2"),Extended(true),Save());
  cout<<"Normalization is: "<<bkg_norm->getVal()<<endl;
  
  std::vector<float> massbin;
  std::vector<TGraph*> signalpullplot_vector;
  
  TString plotmass[6]={"1000N","2600N","4000N","1000W","2600W","4000W"};
  
  for(int j=0; j<6; j++){
    std::vector<float> signalpull;
    float sig_norm = 0;
    TString width = "N";
    TString sigfun = "CBGaus";
    int linestyle = 1;
    int linecolor = 2;
    if(plotmass[j].Contains("W")){
        width = "W";
        sigfun="CB2Gaus";
        linestyle = 2;
    }
    if(plotmass[j] == "1000N") sig_norm=267; //15fb
    if(plotmass[j] == "1000W") sig_norm=257; //15fb
    if(plotmass[j] == "2600N") {sig_norm=15.80; linecolor=4;}//1fb
    if(plotmass[j] == "2600W") {sig_norm=14.42; linecolor=4;}//1fb
    if(plotmass[j] == "4000N") {sig_norm=6.894; linecolor=8;}//0.5fb
    if(plotmass[j] == "4000W") {sig_norm=6.041; linecolor=8;}//0.5fb
    TFile *signal = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/RooFitWorkspace/"+width+"/"+plotmass[j]+"-shapes-Unbinned-"+sigfun+".root");
    RooWorkspace *sig_w = (RooWorkspace*)signal->Get("w");
    RooAbsPdf *sig_pdf = sig_w->pdf("CBGaus");
    for(int i=600; i<8001; i+=20){
      float masspoint = i;
      massbin.push_back(i+10);
      x->setRange("sig",i,i+20);
      RooAbsReal* fracInt_bkg = model->createIntegral(*x,NormSet(*x),Range("sig"));
      RooRealVar* m = sig_w->var("m");
      m->setRange("sig",i,i+20);
      RooAbsReal* fracInt_sig = sig_pdf->createIntegral(*m,NormSet(*m),Range("sig"));
      signalpull.push_back(fracInt_sig->getValV()*sig_norm/sqrt(fracInt_bkg->getValV()*bkg_norm->getVal()));
    }
    TGraph *gr = new TGraph(185,&massbin[0],&signalpull[0]);
    gr->SetLineWidth(2);
    gr->SetLineColor(linecolor);
    gr->SetLineStyle(linestyle);
    signalpullplot_vector.push_back(gr);
  }

  // --- Visualization ---
  gStyle->SetOptStat(111111);
  TCanvas *c01 = new TCanvas("c01","c01",2100,2000);
  TPad *p01a = new TPad("p01a","p01a",0,0.178,1,1.0);
  TPad *p01b = new TPad("p01b","p01b",0,0,1,0.26);
  p01a->Draw();
  p01b->Draw();
  p01a->cd();
  p01a->SetLeftMargin(0.13);
  p01a->SetRightMargin(0.1);
  p01a->SetBottomMargin(0.11);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 40 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.002,100000);
  xaxis->SetLimits(600,5000);
  p01a->SetLogy();
  frame->Draw();
  CMS_lumi(p01a,iPeriod,iPos);
  TLegend *l =  new TLegend(0.48,0.65,0.85,0.85);
  l->AddEntry(frame->findObject("datah"),"Data SR","lep");
  l->AddEntry(frame->findObject(fun_name),"Data SR Fit "+fun_name,"l");
  l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
  
  p01b->cd();
  p01b->SetLeftMargin(0.13);
  p01b->SetRightMargin(0.1);
  p01b->SetTopMargin(0);
  p01b->SetBottomMargin(0.36);
  RooHist* hpull = frame->pullHist();
  RooPlot* pull_frame = x->frame();
  // for(int i=0; i<6; i++){
    // pull_frame->addObject(signalpullplot_vector.at(i),"l");
  // }
  pull_frame->addPlotable(hpull,"P");
  hpull->SetMarkerStyle(8);
  hpull->SetMarkerSize(2);
  xaxis = pull_frame->GetXaxis();
  yaxis = pull_frame->GetYaxis();
  xaxis->SetTitle("M_{j#gamma} (GeV)");
  yaxis->SetTitle("#frac{data-fit}{#sigma_{stat.}}");
  yaxis->SetTitleOffset(0.37);
  yaxis->SetRangeUser(-5,5);
  xaxis->SetLabelSize(0.15);
  xaxis->SetTitleSize(0.15);
  xaxis->SetLimits(600,5000);
  yaxis->SetLabelSize(0.15);
  yaxis->SetTitleSize(0.15);
  yaxis->SetNdivisions(5);
  p01b->SetGrid();
  pull_frame->Draw();
  p01b->Update();
  
  c01->Print(fun_name+".png");
  c01->Print(fun_name+".pdf");
  c01->Print(fun_name+".svg");
  
  frame->GetXaxis()->SetLimits(600,5000);
  pull_frame->GetXaxis()->SetLimits(600,5000);
  c01->Print(fun_name+"_5000plotrange.png");
  c01->Print(fun_name+"_5000plotrange.pdf");
  c01->Print(fun_name+"_5000plotrange.svg");

// --- Output root file ---
  RooWorkspace *w = new RooWorkspace("w","w");
  w->import(*x);
  w->import(data,Rename("data_SR"));
  w->import(datah,Rename("data_SR_binned"));
  w->import(*model);
  w->writeToFile("SR-shapes-Unbinned-"+fun_name+"_5k.root");
  
}
