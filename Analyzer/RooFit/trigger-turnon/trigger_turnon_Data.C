#define fun_type 13
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
void trigger_turnon_Data(int seed=37)
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37); 
  
  lumi_13TeV = "137.19 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Preliminary";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 12;
  int iPos = 0;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.027);
  gStyle->SetBarWidth(1.03);
  gStyle->SetHistLineWidth(2);

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("m","m",400,2000,""); //the name "m" will be used by RooDataSet to import data, normalization range is 600-3500 but plot range can be defined to like 600-3000

  //--- background PDF ---
  
  TString fun_name = "Erf-dijet3";
  RooRealVar *dijet3_p0 = new RooRealVar("dijet3_p0","dijet3_p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *dijet3_p1 = new RooRealVar("dijet3_p1","dijet3_p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *dijet3_p2 = new RooRealVar("dijet3_p2","dijet3_p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar* offset = new RooRealVar("offset","offset",600,400,1000);
  RooRealVar* width  = new RooRealVar("width","width",60,20,100);
  RooGenericPdf * model = new RooGenericPdf(fun_name,fun_name,"pow(m/13000,dijet3_p0+dijet3_p1*log(m/13000)+dijet3_p2*pow(log(m/13000),2))*(1.+TMath::Erf((m-offset)/width))/2.",RooArgList(*x, *dijet3_p0, *dijet3_p1, *dijet3_p2, *offset, *width));

  // --- Import Binned dataset ---
  Float_t p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, evt_pdfunc;
  TH1F Datahist("Data","Data",180,400,4000);
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/ntuple/cutid/Run2Data_nominal_pileup_WGamma_full_full_May22.root");
  TTree *tree = (TTree*)file.Get("Events");
  tree->SetBranchAddress("photon_pt", &p_pt);
  tree->SetBranchAddress("photon_eta", &p_eta);
  tree->SetBranchAddress("photon_phi", &p_phi);
  tree->SetBranchAddress("photon_e", &p_e);
  tree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  tree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  tree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  tree->SetBranchAddress("ak8puppijet_e", &j_e);
  tree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  tree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  tree->SetBranchAddress("sys_costhetastar", &s_cos);
  tree->SetBranchAddress("sys_ptoverm", &s_ptm);
  tree->SetBranchAddress("sys_invmass", &s_mass);
  for (int ievt = 0; ievt<tree->GetEntries();ievt++) {
    tree->GetEntry(ievt);
    if(j_mass < 68 || j_mass > 94) continue;
	if(p_pt < 225) continue;
	if(j_pt < 225) continue;
    if(abs(p_eta) > 1.44) continue;
    if(abs(j_eta) > 2) continue;
    if(j_tau21 > 0.35) continue;
    if(s_ptm < 0.37) continue;
    if(s_cos > 0.6) continue;
    Datahist.Fill(s_mass);
  }
  RooDataHist datah("Datahist","Datahist",RooArgSet(*x),&Datahist);
  cout<<"number of weighted entries: "<<datah.sum(false)<<endl;
  
  // --- Perform extended ML fit of composite PDF to toy data ---
    RooFitResult *r = model->fitTo(datah,Range(400,4000),RooFit::Minimizer("Minuit2"),Save());  //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)
    cout<<0;
  // --- plot for chi2 calculation and visualization ---
  //x->setBins(100); //fit is unbinned but chi2 is calculated by binning data with this value
  //x->setBins(95); //fit is unbinned but chi2 is calculated by binning data with this value
  RooPlot *frame = x->frame();
  cout<<1;
  frame->SetTitle("Data SR");
  datah.plotOn(frame,RooFit::Name("datah"),Binning(180,400,4000),DataError(RooAbsData::Poisson)); //for weighted data
  cout<<2;
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  cout<<3;
  model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  cout<<4;
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  cout<<5;
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  cout<<6;
  datah.plotOn(frame,RooFit::Name("datah"),Binning(180,400,4000),DataError(RooAbsData::Poisson)); //for 
  cout<<7;
  
  // --- Visualization ---
  gStyle->SetOptStat(111111);
  TCanvas *c01 = new TCanvas("c01","c01",2400,1800);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 20 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0,4500);
  xaxis->SetLimits(400,2000);
  //c01->SetLogy();
  frame->Draw();
  CMS_lumi(c01,iPeriod,iPos);
  TLegend *l =  new TLegend(0.5,0.68,0.85,0.8);
  l->AddEntry(frame->findObject(fun_name),"Data SR fit "+fun_name,"l");
  l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
  c01->Print(fun_name+".png");
  c01->Print(fun_name+".pdf");
  c01->Print(fun_name+".svg");
  
  // TCanvas *c02 = new TCanvas("c02","c02",1200,300);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  // c02->cd();
  // RooHist* hpull = frame->pullHist();
  // RooPlot* pull_frame = x->frame();
  // pull_frame->addPlotable(hpull,"P") ;
  // hpull->SetMarkerStyle(8);
  // xaxis = pull_frame->GetXaxis();
  // yaxis = pull_frame->GetYaxis();
  // xaxis->SetTitle("M_{j#gamma}");
  // yaxis->SetTitle("#frac{MC-fit}{#sigma_{stat}}");
  // yaxis->SetTitleOffset(0.5);
  // yaxis->SetRangeUser(-5,5);
  // xaxis->SetLabelSize(0.08);
  // xaxis->SetTitleSize(0.08);
  // xaxis->SetRangeUser(0,2000);
  // yaxis->SetLabelSize(0.08);
  // yaxis->SetTitleSize(0.08);
  // yaxis->SetNdivisions(5);
  // c02->SetBottomMargin(0.18);
  // c02->SetTopMargin(0.18);
  // c02->SetGrid();
  // pull_frame->Draw();
  // c02->Update();
  // c02->Print("pull"+fun_name+".png");
  // c02->Print("pull"+fun_name+".pdf"); 
  // c02->Print("pull"+fun_name+".svg"); 

  TCanvas *c03 = new TCanvas("c03","c03",1200,900);
  c03->cd();
  TF1 *f1 = new TF1("erf","3900*(1.+TMath::Erf((x-[0])/[1]))/2.",400,2000);
  f1->SetParameters(offset->getValV(), width->getValV());
  TF1 *f2 = new TF1("exp","5*pow(10,-7)*pow(x/13000,[0]+[1]*log(x/13000)+[2]*pow(log(x/13000),2))",400,2000);
  f2->SetParameters(dijet3_p0->getValV(), dijet3_p1->getValV(), dijet3_p2->getValV());
  cout<<"Turn-on end point is: "<<f1->GetX(0.99)<<endl;
  f1->GetYaxis()->SetRangeUser(0,4500);
  f1->Draw();
  f2->Draw("SAME");
  c03->Print("pdf.png");
  c03->Print("pdf.svg");
}
