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

void Pileupuncertainty(int sigm = 2800)
{

  gROOT->SetBatch(1);
  lumi_13TeV = "35.92 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.8;
  cmsTextSize = 1;
  int iPeriod = 4;
  int iPos = 0;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.09,"XYZ");
  gStyle->SetLabelSize(0.09,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.06);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(1);
  // Plots
  
  //16
  double nominal_N[15] = {3141,3058,3130,3046,3090,3028,2936,2793,2643,2724,2553,2532,2449,2403,2452};
  double up_N[15] = {3132,3035,3109,3048,3076,3011,2924,2780,2621,2713,2533,2526,2443,2385,2447};
  double down_N[15] = {3147,3084,3149,3047,3103,3044,2946,2805,2664,2736,2572,2541,2455,2420,2454};
  
  double nominal_W[14] = {2911,2858,2938,3004,2919,2865,2740,2765,2557,2432,2452,2194,2182,2078};
  double up_W[14] = {2904,2842,2917,2988,2927,2854,2728,2755,2555,2422,2436,2187,2174,2069};
  double down_W[14] = {2918,2873,2959,3018,2911,2875,2749,2772,2556,2445,2471,2201,2189,2086};
  
  double massN[15] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500};
  double massW[14] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2800,3000,3500};
  
  // 17
  // double nominal_N[14] = {3225,3270,3316,3397,3271,3232,3170,2933,2796,2740,2662,2699,2536,2510};
  // double up_N[14] = {3185,3241,3306,3355,3245,3188,3129,2905,2767,2714,2636,2665,2520,2497};
  // double down_N[14] = {3255,3302,3328,3435,3295,3266,3207,2956,2819,2761,2682,2727,2552,2527};
  
  // double nominal_W[14] = {3174,3157,3224,3272,3048,3059,2878,2913,2726,2671,2506,2472,2364,2108};
  // double up_W[14] = {3153,3141,3204,3250,2999,3019,2846,2882,2701,2636,2450,2469,2355,2087};
  // double down_W[14] = {3194,3166,3250,3294,3090,3099,2914,2940,2752,2707,2554,2473,2374,2121};
  
  // double massN[14] = {700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3500};
  // double massW[14] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3500};
  
  //18
  // double nominal_N[3] = {3247,3040,3088};
  // double up_N[3] = {3246,3028,3081};
  // double down_N[3] = {3257,3055,3090};
  
  // double nominal_W[4] = {3268,3166,3122,2995};
  // double up_W[4] = {3255,3154,3116,2980};
  // double down_W[4] = {3277,3173,3131,3008};
  
  // double massN[3] = {1200,1800,2000};
  // double massW[4] = {900,1000,1200,1600};
  
  TGraph Nup;
  TGraph Ndown;
  TGraph Wup;
  TGraph Wdown;
  
  for(int i=0; i<14; i++){
	  Nup.SetPoint(i,massN[i],up_N[i]/nominal_N[i] - 1);
	  Ndown.SetPoint(i,massN[i],down_N[i]/nominal_N[i] - 1);
  }
  
    for(int i=0; i<14; i++){
	  Wup.SetPoint(i,massW[i],up_W[i]/nominal_W[i] - 1);
	  Wdown.SetPoint(i,massW[i],down_W[i]/nominal_W[i] - 1);
  }
  
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;
  gStyle->SetBarWidth(0.5);
  TCanvas *c0 = new TCanvas("c0","",2400,2000);
  c0->Divide(1,4);
  c0->cd(1);
  xaxis = Nup.GetXaxis();
  yaxis = Nup.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  Nup.SetFillColor(2);
  //Nup.SetFillColorAlpha(kViolet,0.5);
  Nup.Draw("AB");
  TLegend *l1 = new TLegend(0.70,0.75,0.9,0.9);
  l1->AddEntry(&Nup,"Pileup reweighting UP (N)","f");
  l1->Draw();
  c0->cd(2);
  xaxis = Ndown.GetXaxis();
  yaxis = Ndown.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  Ndown.SetFillColor(2);
  //Ndown.SetFillColorAlpha(kViolet,0.5);
  Ndown.Draw("AB");
  TLegend *l2 = new TLegend(0.70,0.75,0.9,0.9);
  l2->AddEntry(&Ndown,"Pileup reweighting DOWN (N)","f");
  l2->Draw();
  c0->cd(3);
  xaxis = Wup.GetXaxis();
  yaxis = Wup.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  Wup.SetFillColor(8);
  //Wup.SetFillColorAlpha(kGreen+3,0.5);
  Wup.Draw("AB");
  TLegend *l3 = new TLegend(0.70,0.75,0.9,0.9);
  l3->AddEntry(&Wup,"Pileup reweighting UP (W)","f");
  l3->Draw();
  c0->cd(4);
  xaxis = Wdown.GetXaxis();
  yaxis = Wdown.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  Wdown.SetFillColor(8);
  //Wdown.SetFillColorAlpha(8,0.5);
  Wdown.Draw("AB");
  TLegend *l4 = new TLegend(0.70,0.75,0.9,0.9);
  l4->AddEntry(&Wdown,"Pileup reweighting DOWN (W)","f");
  l4->Draw();
  ((TPad*)(c0->GetPrimitive("c0_1")))->SetGrid();
  ((TPad*)(c0->GetPrimitive("c0_2")))->SetGrid();
  ((TPad*)(c0->GetPrimitive("c0_3")))->SetGrid();
  ((TPad*)(c0->GetPrimitive("c0_4")))->SetGrid();
  ((TPad*)(c0->GetPrimitive("c0_1")))->SetTopMargin(0.1);
  ((TPad*)(c0->GetPrimitive("c0_2")))->SetTopMargin(0.1);
  ((TPad*)(c0->GetPrimitive("c0_3")))->SetTopMargin(0.1);
  ((TPad*)(c0->GetPrimitive("c0_4")))->SetTopMargin(0.1);
  ((TPad*)(c0->GetPrimitive("c0_1")))->SetBottomMargin(0.2);
  ((TPad*)(c0->GetPrimitive("c0_2")))->SetBottomMargin(0.2);
  ((TPad*)(c0->GetPrimitive("c0_3")))->SetBottomMargin(0.2);
  ((TPad*)(c0->GetPrimitive("c0_4")))->SetBottomMargin(0.2);
  CMS_lumi(((TPad*)(c0->GetPrimitive("c0_1"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c0->GetPrimitive("c0_2"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c0->GetPrimitive("c0_3"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c0->GetPrimitive("c0_4"))),iPeriod,iPos);
  c0->Print("pileup.png");
  c0->Print("pileup.pdf");
  c0->Print("pileup.svg");
  c0->Print("pileup.root");
}