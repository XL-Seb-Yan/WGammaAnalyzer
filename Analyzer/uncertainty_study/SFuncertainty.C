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

void SFuncertainty(int sigm = 2800)
{

  gROOT->SetBatch(1);
  lumi_13TeV = "59.74 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.8;
  cmsTextSize = 1;
  int iPeriod = 6;
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
  // double ph_energy_up_N[15] = {3122,3114,3146,3072,3096,3056,2956,2821,2679,2762,2594,2573,2451,2424,2455};
  // double ph_energy_up_W[14] = {2914,2874,2979,3015,2911,2881,2755,2775,2537,2466,2519,2211,2197,2090};
  // double ph_energy_down_N[15] = {3118,3106,3145,3064,3085,3039,2940,2805,2673,2751,2588,2565,2445,2413,2448};
  // double ph_energy_down_W[14] = {2913,2872,2975,3018,2896,2879,2739,2760,2522,2458,2506,2199,2181,2080};
  // double jet_energy_up_N[15] = {3103,3091,3135,3056,3089,3037,2942,2807,2684,2757,2592,2569,2463,2420,2451};
  // double jet_energy_up_W[14] = {2895,2861,2958,3006,2896,2874,2743,2775,2527,2468,2511,2204,2191,2087};
 
  // double jet_energy_down_N[15] = {3141,3138,3154,3075,3098,3049,2947,2820,2668,2746,2578,2572,2444,2418,2456};
  // double jet_energy_down_W[14] = {2931,2884,2989,3030,2912,2879,2750,2771,2521,2450,2509,2208,2192,2081};
  // double jet_res_up_N[15] = {3095,3088,3137,3049,3074,3029,2939,2811,2666,2747,2578,2564,2446,2418,2446};
  // double jet_res_up_W[14] = {2886,2856,2955,3009,2889,2866,2733,2756,2518,2446,2508,2191,2177,2079};
  // double jet_res_down_N[15] = {3145,3131,3158,3083,3102,3061,2953,2817,2687,2761,2595,2578,2456,2432,2455};
  // double jet_res_down_W[14] = {2937,2888,2997,3037,2916,2889,2754,2774,2534,2468,2521,2216,2209,2084};
  // double nominal_N[15] = {3119,3109,3148,3069,3093,3049,2948,2816,2677,2756,2592,2568,2448,2420,2450};
  // double nominal_W[14] = {2915,2871,2978,3016,2902,2881,2748,2766,2530,2459,2513,2204,2187,2086};
  
  // double massN[15] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500};
  // double massW[14] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2800,3000,3500};
  
  // 17
  // double ph_energy_up_N[14] = {3212,3275,3302,3391,3280,3221,3193,2926,2821,2759,2648,2710,2543,2523};
  // double ph_energy_up_W[14] = {3171,3153,3244,3311,3072,3087,2929,2941,2742,2690,2558,2465,2368,2093};
  // double ph_energy_down_N[14] = {3192,3259,3285,3374,3263,3210,3162,2905,2808,2742,2640,2700,2529,2518};
  // double ph_energy_down_W[14] = {3166,3146,3232,3300,3058,3071,2912,2913,2720,2678,2521,2450,2354,2079};
  // double jet_energy_up_N[14] = {3197,3262,3300,3373,3278,3221,3171,2927,2817,2755,2649,2711,2545,2524};
  // double jet_energy_up_W[14] = {3162,3139,3234,3307,3064,3078,2923,2930,2735,2681,2546,2463,2362,2094};
  // double jet_energy_down_N[14] = {3212,3277,3306,3383,3268,3214,3174,2913,2813,2736,2639,2705,2536,2514};
  // double jet_energy_down_W[14] = {3173,3157,3242,3307,3068,3090,2918,2929,2725,2676,2532,2452,2353,2084};
  // double jet_res_up_N[14] = {3194,3251,3280,3365,3261,3198,3168,2916,2803,2735,2640,2692,2536,2515};
  // double jet_res_up_W[14] = {3145,3121,3212,3286,3052,3072,2910,2916,2724,2673,2534,2447,2355,2084};
  // double jet_res_down_N[14] = {3238,3291,3320,3406,3286,3237,3204,2927,2823,2763,2652,2720,2550,2531};
  // double jet_res_down_W[14] = {3191,3172,3260,3342,3073,3105,2932,2950,2738,2681,2544,2478,2372,2101};
  // double nominal_N[14] = {3208,3270,3294,3383,3272,3217,3178,2922,2813,2751,2643,2705,2535,2523};
  // double nominal_W[14] = {3167,3143,3233,3304,3063,3081,2921,2930,2732,2683,2537,2458,2361,2088};
  
  // double massN[14] = {700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3500};
  // double massW[14] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3500};
  
  //18
  double ph_energy_up_N[3] = {3307,3085,3095};
  double ph_energy_up_W[4] = {3272,3189,3152,3028};
  double ph_energy_down_N[3] = {3273,3061,3069};
  double ph_energy_down_W[4] = {3262,3169,3130,3019};
  double jet_energy_up_N[3] = {3296,3078,3087};
  double jet_energy_up_W[4] = {3260,3181,3133,3030};
 
  double jet_energy_down_N[3] = {3303,3076,3082};
  double jet_energy_down_W[4] = {3266,3190,3146,3027};
  double jet_res_up_N[3] = {3293,3074,3078};
  double jet_res_up_W[4] = {3258,3175,3126,3022};
  double jet_res_down_N[3] = {3304,3079,3091};
  double jet_res_down_W[4] = {3279,3192,3152,3036};
  double nominal_N[3] = {3298,3077,3081};
  double nominal_W[4] = {3269,3182,3143,3024};
  
  double massN[3] = {1200,1800,2000};
  double massW[4] = {900,1000,1200,1600};
  
  TGraph ph_corrNup;
  TGraph ph_corrNdown;
  TGraph ph_corrWup;
  TGraph ph_corrWdown;
  TGraph jet_corrNup;
  TGraph jet_corrNdown;
  TGraph jet_corrWup;
  TGraph jet_corrWdown;
  TGraph jet_resNup;
  TGraph jet_resWup;
  TGraph jet_resNdown;
  TGraph jet_resWdown;
  
  for(int i=0; i<3; i++){
	  ph_corrNup.SetPoint(i,massN[i],ph_energy_up_N[i]/nominal_N[i] - 1);
	  ph_corrNdown.SetPoint(i,massN[i],ph_energy_down_N[i]/nominal_N[i] - 1);
	  jet_corrNup.SetPoint(i,massN[i],jet_energy_up_N[i]/nominal_N[i] - 1);
	  jet_corrNdown.SetPoint(i,massN[i],jet_energy_down_N[i]/nominal_N[i] - 1);
	  jet_resNup.SetPoint(i,massN[i],jet_res_up_N[i]/nominal_N[i] - 1);
	  jet_resNdown.SetPoint(i,massN[i],jet_res_down_N[i]/nominal_N[i] - 1);
	  
  }
  
    for(int i=0; i<4; i++){
	  ph_corrWup.SetPoint(i,massW[i],ph_energy_up_W[i]/nominal_W[i] - 1);
	  ph_corrWdown.SetPoint(i,massW[i],ph_energy_down_W[i]/nominal_W[i] - 1);
	  jet_corrWup.SetPoint(i,massW[i],jet_energy_up_W[i]/nominal_W[i] - 1);
	  jet_corrWdown.SetPoint(i,massW[i],jet_energy_down_W[i]/nominal_W[i] - 1);
	  jet_resWup.SetPoint(i,massW[i],jet_res_up_W[i]/nominal_W[i] - 1);
	  jet_resWdown.SetPoint(i,massW[i],jet_res_down_W[i]/nominal_W[i] - 1);
	  
  }
  
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;
  gStyle->SetBarWidth(0.5);
  TCanvas *c0 = new TCanvas("c0","",2400,2000);
  c0->Divide(1,4);
  c0->cd(1);
  xaxis = ph_corrNup.GetXaxis();
  yaxis = ph_corrNup.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  ph_corrNup.SetFillColor(2);
  //ph_corrNup.SetFillColorAlpha(kViolet,0.5);
  ph_corrNup.Draw("AB");
  TLegend *l01 = new TLegend(0.7,0.75,0.9,0.9);
  l01->AddEntry(&ph_corrNup,"Photon energy corr. UP (N)","f");
  l01->Draw();
  c0->cd(2);
  xaxis = ph_corrNdown.GetXaxis();
  yaxis = ph_corrNdown.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  ph_corrNdown.SetFillColor(2);
  //ph_corrNdown.SetFillColorAlpha(kViolet,0.5);
  ph_corrNdown.Draw("AB");
  TLegend *l02 = new TLegend(0.7,0.75,0.9,0.9);
  l02->AddEntry(&ph_corrNdown,"Photon energy corr. DOWN (N)","f");
  l02->Draw();
  c0->cd(3);
  xaxis = ph_corrWup.GetXaxis();
  yaxis = ph_corrWup.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  ph_corrWup.SetFillColor(8);
  //ph_corrWup.SetFillColorAlpha(kGreen+3,0.5);
  ph_corrWup.Draw("AB");
  TLegend *l03 = new TLegend(0.7,0.75,0.9,0.9);
  l03->AddEntry(&ph_corrWup,"Photon energy corr. UP (W)","f");
  l03->Draw();
  c0->cd(4);
  xaxis = ph_corrWdown.GetXaxis();
  yaxis = ph_corrWdown.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  ph_corrWdown.SetFillColor(8);
  //ph_corrWdown.SetFillColorAlpha(kGreen+3,0.5);
  ph_corrWdown.Draw("AB");
  TLegend *l04 = new TLegend(0.7,0.75,0.9,0.9);
  l04->AddEntry(&ph_corrWdown,"Photon energy corr. DOWN (W)","f");
  l04->Draw();
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
  c0->Print("ph_corr.png");
  c0->Print("ph_corr.pdf");
  c0->Print("ph_corr.svg");
  c0->Print("ph_corr.root");
  
  TCanvas *c1 = new TCanvas("c1","",2400,2000);
  c1->Divide(1,4);
  c1->cd(1);
  xaxis = jet_corrNup.GetXaxis();
  yaxis = jet_corrNup.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  jet_corrNup.SetFillColor(2);
  //jet_corrNup.SetFillColorAlpha(kViolet,0.5);
  jet_corrNup.Draw("AB");
  TLegend *l11 = new TLegend(0.7,0.75,0.9,0.9);
  l11->AddEntry(&jet_corrNup,"Jet energy corr. UP (N)","f");
  l11->Draw();
  c1->cd(2);
  xaxis = jet_corrNdown.GetXaxis();
  yaxis = jet_corrNdown.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  jet_corrNdown.SetFillColor(2);
  //jet_corrNdown.SetFillColorAlpha(2,0.9);
  jet_corrNdown.Draw("AB");
  TLegend *l12 = new TLegend(0.7,0.75,0.9,0.9);
  l12->AddEntry(&jet_corrNdown,"Jet energy corr. DOWN (N)","f");
  l12->Draw();
  c1->cd(3);
  xaxis = jet_corrWup.GetXaxis();
  yaxis = jet_corrWup.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  jet_corrWup.SetFillColor(8);
  //jet_corrWup.SetFillColorAlpha(8,0.9);
  jet_corrWup.Draw("AB");
  TLegend *l13 = new TLegend(0.7,0.75,0.9,0.9);
  l13->AddEntry(&jet_corrWup,"Jet energy corr. UP (W)","f");
  l13->Draw();
  c1->cd(4);
  xaxis = jet_corrWdown.GetXaxis();
  yaxis = jet_corrWdown.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  jet_corrWdown.SetFillColor(8);
  //jet_corrWdown.SetFillColorAlpha(8,0.9);
  jet_corrWdown.Draw("AB");
  TLegend *l14 = new TLegend(0.7,0.75,0.9,0.9);
  l14->AddEntry(&jet_corrWdown,"Jet energy corr. DOWN (W)","f");
  l14->Draw();
  ((TPad*)(c1->GetPrimitive("c1_1")))->SetGrid();
  ((TPad*)(c1->GetPrimitive("c1_2")))->SetGrid();
  ((TPad*)(c1->GetPrimitive("c1_3")))->SetGrid();
  ((TPad*)(c1->GetPrimitive("c1_4")))->SetGrid();
  ((TPad*)(c1->GetPrimitive("c1_1")))->SetTopMargin(0.1);
  ((TPad*)(c1->GetPrimitive("c1_2")))->SetTopMargin(0.1);
  ((TPad*)(c1->GetPrimitive("c1_3")))->SetTopMargin(0.1);
  ((TPad*)(c1->GetPrimitive("c1_4")))->SetTopMargin(0.1);
  ((TPad*)(c1->GetPrimitive("c1_1")))->SetBottomMargin(0.2);
  ((TPad*)(c1->GetPrimitive("c1_2")))->SetBottomMargin(0.2);
  ((TPad*)(c1->GetPrimitive("c1_3")))->SetBottomMargin(0.2);
  ((TPad*)(c1->GetPrimitive("c1_4")))->SetBottomMargin(0.2);
  CMS_lumi(((TPad*)(c1->GetPrimitive("c1_1"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c1->GetPrimitive("c1_2"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c1->GetPrimitive("c1_3"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c1->GetPrimitive("c1_4"))),iPeriod,iPos);
  c1->Print("jet_corr.png");
  c1->Print("jet_corr.pdf");
  c1->Print("jet_corr.svg");
  c1->Print("jet_corr.root");
  
  TCanvas *c2 = new TCanvas("c2","",2400,2000);
  c2->Divide(1,4);
  c2->cd(1);
  xaxis = jet_resNup.GetXaxis();
  yaxis = jet_resNup.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  jet_resNup.SetFillColor(2);
  //jet_resNup.SetFillColorAlpha(2,0.9);
  jet_resNup.Draw("AB");
  TLegend *l21 = new TLegend(0.7,0.75,0.9,0.9);
  l21->AddEntry(&jet_resNup,"Jet energy res. UP (N)","f");
  l21->Draw();
  c2->cd(2);
  xaxis = jet_resNdown.GetXaxis();
  yaxis = jet_resNdown.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  jet_resNdown.SetFillColor(2);
  //jet_resNdown.SetFillColorAlpha(2,0.9);
  jet_resNdown.Draw("AB");
  TLegend *l22 = new TLegend(0.7,0.75,0.9,0.9);
  l22->AddEntry(&jet_resNdown,"Jet energy res. DOWN (N)","f");
  l22->Draw();
  c2->cd(3);
  xaxis = jet_resWup.GetXaxis();
  yaxis = jet_resWup.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  jet_resWup.SetFillColor(8);
 // jet_resWup.SetFillColorAlpha(8,0.9);
  jet_resWup.Draw("AB");
  TLegend *l23 = new TLegend(0.7,0.75,0.9,0.9);
  l23->AddEntry(&jet_resWup,"Jet energy res. UP (W)","f");
  l23->Draw();
  c2->cd(4);
  xaxis = jet_resWdown.GetXaxis();
  yaxis = jet_resWdown.GetYaxis();
  xaxis->SetLimits(500,4000);
  xaxis->SetLabelOffset(0.03);
  xaxis->SetTitle("M_{x} (GeV)");
  yaxis->SetRangeUser(-0.02,0.02);
  yaxis->SetTitle("Variation");
  yaxis->SetTitleOffset(0.4);
  yaxis->SetNdivisions(5);
  jet_resWdown.SetFillColor(8);
  //jet_resWdown.SetFillColorAlpha(8,0.9);
  jet_resWdown.Draw("AB");
  TLegend *l24 = new TLegend(0.7,0.75,0.9,0.9);
  l24->AddEntry(&jet_resWdown,"Jet energy res. DOWN (W)","f");
  l24->Draw();
  ((TPad*)(c2->GetPrimitive("c2_1")))->SetGrid();
  ((TPad*)(c2->GetPrimitive("c2_2")))->SetGrid();
  ((TPad*)(c2->GetPrimitive("c2_3")))->SetGrid();
  ((TPad*)(c2->GetPrimitive("c2_4")))->SetGrid();
  ((TPad*)(c2->GetPrimitive("c2_1")))->SetTopMargin(0.1);
  ((TPad*)(c2->GetPrimitive("c2_2")))->SetTopMargin(0.1);
  ((TPad*)(c2->GetPrimitive("c2_3")))->SetTopMargin(0.1);
  ((TPad*)(c2->GetPrimitive("c2_4")))->SetTopMargin(0.1);
  ((TPad*)(c2->GetPrimitive("c2_1")))->SetBottomMargin(0.2);
  ((TPad*)(c2->GetPrimitive("c2_2")))->SetBottomMargin(0.2);
  ((TPad*)(c2->GetPrimitive("c2_3")))->SetBottomMargin(0.2);
  ((TPad*)(c2->GetPrimitive("c2_4")))->SetBottomMargin(0.2);
  CMS_lumi(((TPad*)(c2->GetPrimitive("c2_1"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c2->GetPrimitive("c2_2"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c2->GetPrimitive("c2_3"))),iPeriod,iPos);
  CMS_lumi(((TPad*)(c2->GetPrimitive("c2_4"))),iPeriod,iPos);
  c2->Print("jet_res.png");
  c2->Print("jet_res.pdf");
  c2->Print("jet_res.svg");
  c2->Print("jet_res.root");
}