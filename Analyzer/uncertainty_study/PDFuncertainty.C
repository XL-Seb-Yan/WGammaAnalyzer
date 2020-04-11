//MC selection
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
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

void PDFuncertainty(int sigm = 2800)
{

  gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.4;
  cmsTextSize = 0.7;
  int iPeriod = 5;
  int iPos = 0;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetLabelSize(0.06,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(1);
  // Plots
  
  // 17
  double totalweight1N[5] = {19747.2,20653.6,21540.7,22130.7,22551.4};
  double totalweight1W[5] = {19731.2,20612.1,21467.2,22008.1,22317.9};
  double PDF1_N[5] = {3179.05,3347.91,3120.7,3013.15,2878.2};
  double PDF1_W[5] = {3127.35,3152.14,2951.7,2623.12,2379.62};
  double nominal_N[5] = {3229,3250,2896,2700,2523};
  double nominal_W[5] = {3177,3062,2739,2354,2090};
  
  double massN[5] = {700,1200,2000,2800,3500};
  double massW[5] = {700,1200,2000,2800,3500};
  
  TGraph PDF1N;
  TGraph PDF1W;
  
  for(int i=0; i<5; i++){
	  PDF1N.SetPoint(i,massN[i],(PDF1_N[i] / totalweight1N[i])/(nominal_N[i]/20000) - 1);
	  PDF1W.SetPoint(i,massW[i],(PDF1_W[i] / totalweight1W[i])/(nominal_W[i]/20000) - 1);
  }
  
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;
  gStyle->SetBarWidth(0.5);
  TCanvas *c0 = new TCanvas("c0","",2400,1800);
  c0->Divide(1,2);
  c0->cd(1);
  xaxis = PDF1N.GetXaxis();
  yaxis = PDF1N.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitleOffset(1.2);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.05,0.05);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.7);
  yaxis->SetNdivisions(5);
  PDF1N.SetFillColor(2);
  //PDF1N.SetFillColorAlpha(kViolet,0.5);
  PDF1N.Draw("AB");
  TLegend *l01 = new TLegend(0.65,0.8,0.9,0.9);
  l01->AddEntry(&PDF1N,"NNPDF31_nnlo_hessian replicas (N)","f");
  l01->Draw();
  c0->cd(2);
  xaxis = PDF1W.GetXaxis();
  yaxis = PDF1W.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitleOffset(1.2);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.05,0.05);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.7);
  yaxis->SetNdivisions(5);
  PDF1W.SetFillColor(8);
  //PDF1W.SetFillColorAlpha(kGreen+3,0.5);
  PDF1W.Draw("AB");
  TLegend *l03 = new TLegend(0.65,0.8,0.9,0.9);
  l03->AddEntry(&PDF1W,"NNPDF31_nnlo_hessian replicas (W)","f");
  l03->Draw();
  ((TPad*)(c0->GetPrimitive("c0_1")))->SetGrid();
  ((TPad*)(c0->GetPrimitive("c0_2")))->SetGrid();
  ((TPad*)(c0->GetPrimitive("c0_1")))->SetTopMargin(0.1);
  ((TPad*)(c0->GetPrimitive("c0_2")))->SetTopMargin(0.1);
  ((TPad*)(c0->GetPrimitive("c0_1")))->SetBottomMargin(0.2);
  ((TPad*)(c0->GetPrimitive("c0_2")))->SetBottomMargin(0.2);
  CMS_lumi(((TPad*)(c0->GetPrimitive("c0_1"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c0->GetPrimitive("c0_2"))),iPeriod,iPos);
  c0->Print("PDFUnc.png");
  c0->Print("PDFUnc.pdf");
  c0->Print("PDFUnc.svg");
  c0->Print("PDFUnc.root");
}