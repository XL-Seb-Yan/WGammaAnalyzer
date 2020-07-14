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

void compare_plot_2_mc()
{
  gROOT->SetBatch(1);
  lumi_13TeV = "";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  relPosX = 0.11;
  bool plot_CMS = true;
  extraText = "Preliminary";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 12;
  int iPos = 0;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.035,"XYZ");
  gStyle->SetLabelSize(0.035,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);

  TH1 *hist11 = new TH1F("11","pt_{#gamma}",100,0,3000);
  TH1 *hist12 = new TH1F("12","eta_{#gamma}",50,-2,2);
  TH1 *hist13 = new TH1F("13","pt_{j}",100,0,3000);
  TH1 *hist14 = new TH1F("14","eta_{j}",50,-2,2);
  TH1 *hist15 = new TH1F("15","E_{j}",100,0,3000);
  TH1 *hist16 = new TH1F("16","masssoftdrop_{j}",40,40,120);
  TH1 *hist17 = new TH1F("17","tau21_{j}",50,0,1);
  TH1 *hist18 = new TH1F("18","cos(#theta*)_{p}",50,0,1);
  TH1 *hist19 = new TH1F("19","pt/M",50,0,2);
  TH1 *hist110 = new TH1F("110","invariant mass",100,0,4000);
  TH1 *hist111 = new TH1F("111","#gamma Data SR",60,-0.2,1);
  TH1 *hist21 = new TH1F("21","pt_{#gamma}",100,0,3000);
  TH1 *hist22 = new TH1F("22","eta_{#gamma}",50,-2,2);
  TH1 *hist23 = new TH1F("23","pt_{j}",100,0,3000);
  TH1 *hist24 = new TH1F("24","eta_{j}",50,-2,2);
  TH1 *hist25 = new TH1F("25","E_{j}",100,0,3000);
  TH1 *hist26 = new TH1F("26","masssoftdrop_{j}",40,40,120);
  TH1 *hist27 = new TH1F("27","tau21_{j}",50,0,1);
  TH1 *hist28 = new TH1F("28","cos(#theta*)_{p}",50,0,1);
  TH1 *hist29 = new TH1F("29","pt/M",50,0,2);
  TH1 *hist210 = new TH1F("210","invariant mass",100,0,4000);
  TH1 *hist211 = new TH1F("211","#gamma Data SR",60,-0.2,1);
  // Local variables to store to outfile
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21;
  // System
  float sys_costhetastar, sys_ptoverm, m;
  
  // Open input file
  Float_t p_pt, p_eta, p_phi, p_e, p_mva, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight, x_puweight, x_kfactor, x_sf;;
  
  TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/fullcut/Run2Data_postproc_WGammaRun2_SR_sigrange_fullcut_jmcorr_May22.root");
  TTree* theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  //theTree->SetBranchAddress("xsec_weight", &x_weight);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    
    hist11->Fill(p_pt);
    hist12->Fill(p_eta);
    hist13->Fill(j_pt);
    hist14->Fill(j_eta);
    hist15->Fill(j_e);
    hist16->Fill(j_mass);
    hist17->Fill(j_tau21);
    hist18->Fill(s_cos);
    hist19->Fill(s_ptm);
    hist110->Fill(s_mass);
    hist111->Fill(p_mva);  
  }

  input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/fullcut/BkgMC_postproc_WGamma17_SR_full_fullcut_jmcorr_May22.root");
  theTree = (TTree*)input->Get("Events");
  // Improt variables for cutting
  theTree->SetBranchAddress("photon_pt", &p_pt);
  theTree->SetBranchAddress("photon_eta", &p_eta);
  theTree->SetBranchAddress("photon_phi", &p_phi);
  theTree->SetBranchAddress("photon_e", &p_e);
  theTree->SetBranchAddress("photon_mvaval", &p_mva);
  theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
  theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
  theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
  theTree->SetBranchAddress("ak8puppijet_e", &j_e);
  theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
  theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
  theTree->SetBranchAddress("sys_costhetastar", &s_cos);
  theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
  theTree->SetBranchAddress("m", &s_mass);
  theTree->SetBranchAddress("xsec_weight", &x_weight);
  theTree->SetBranchAddress("xsec_kfactor", &x_kfactor);
  theTree->SetBranchAddress("xsec_puweight", &x_puweight);
  theTree->SetBranchAddress("xsec_sf", &x_sf);
  
  double sumW = 0;
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    
    float weight = x_weight*x_puweight*x_kfactor*x_sf;
    
    hist21->Fill(p_pt, weight);
    hist22->Fill(p_eta, weight);
    hist23->Fill(j_pt, weight);
    hist24->Fill(j_eta, weight);
    hist25->Fill(j_e, weight);
    hist26->Fill(j_mass, weight);
    hist27->Fill(j_tau21, weight);
    hist28->Fill(s_cos, weight);
    hist29->Fill(s_ptm, weight);
    hist210->Fill(s_mass, weight);
    hist211->Fill(p_mva, weight);  
    
    if(s_mass > 2000)
        sumW+=weight;
  }

  //=================================================================================
  int color1 = 4;
  int color2 = 2;
  hist11->SetLineColor(color1);
  hist12->SetLineColor(color1);
  hist13->SetLineColor(color1);
  hist14->SetLineColor(color1);
  hist15->SetLineColor(color1);
  hist16->SetLineColor(color1);
  hist17->SetLineColor(color1);
  hist18->SetLineColor(color1);
  hist19->SetLineColor(color1);
  hist110->SetLineColor(color1);
  hist111->SetLineColor(color1);
  hist21->SetLineColor(color2);
  hist22->SetLineColor(color2);
  hist23->SetLineColor(color2);
  hist24->SetLineColor(color2);
  hist25->SetLineColor(color2);
  hist26->SetLineColor(color2);
  hist27->SetLineColor(color2);
  hist28->SetLineColor(color2);
  hist29->SetLineColor(color2);
  hist210->SetLineColor(color2);
  hist211->SetLineColor(color2);

  // hist11->SetLineStyle(2);
  // hist12->SetLineStyle(2);
  // hist13->SetLineStyle(2);
  // hist14->SetLineStyle(2);
  // hist15->SetLineStyle(2);
  // hist16->SetLineStyle(2);
  // hist17->SetLineStyle(2);
  // hist18->SetLineStyle(2);
  // hist19->SetLineStyle(2);
  // hist110->SetLineStyle(2);
  // hist111->SetLineStyle(2);
  // hist21->SetLineStyle(3);
  // hist22->SetLineStyle(3);
  // hist23->SetLineStyle(3);
  // hist24->SetLineStyle(3);
  // hist25->SetLineStyle(3);
  // hist26->SetLineStyle(3);
  // hist27->SetLineStyle(3);
  // hist28->SetLineStyle(3);
  // hist29->SetLineStyle(3);
  // hist210->SetLineStyle(3);
  // hist211->SetLineStyle(3);
  // hist31->SetLineStyle(2);
  // hist32->SetLineStyle(2);
  // hist33->SetLineStyle(2);
  // hist34->SetLineStyle(2);
  // hist35->SetLineStyle(2);
  // hist36->SetLineStyle(2);
  // hist37->SetLineStyle(2);
  // hist38->SetLineStyle(2);
  // hist39->SetLineStyle(2);
  // hist310->SetLineStyle(2);
  // hist311->SetLineStyle(2);
  
  hist11->SetLineWidth(3);
  hist12->SetLineWidth(3);
  hist13->SetLineWidth(3);
  hist14->SetLineWidth(3);
  hist15->SetLineWidth(3);
  hist16->SetLineWidth(3);
  hist17->SetLineWidth(3);
  hist18->SetLineWidth(3);
  hist19->SetLineWidth(3);
  hist110->SetLineWidth(3);
  hist111->SetLineWidth(3);
  hist21->SetLineWidth(3);
  hist22->SetLineWidth(3);
  hist23->SetLineWidth(3);
  hist24->SetLineWidth(3);
  hist25->SetLineWidth(3);
  hist26->SetLineWidth(3);
  hist27->SetLineWidth(3);
  hist28->SetLineWidth(3);
  hist29->SetLineWidth(3);
  hist210->SetLineWidth(3);
  hist211->SetLineWidth(3);
  
  //fix strange error bars...
  // hist11->Sumw2();
  // hist12->Sumw2();
  // hist13->Sumw2();
  // hist14->Sumw2();
  // hist15->Sumw2();
  // hist16->Sumw2();
  // hist17->Sumw2();
  // hist18->Sumw2();
  // hist19->Sumw2();
  // hist110->Sumw2();
  // hist111->Sumw2();
  // hist21->Sumw2();
  // hist22->Sumw2();
  // hist23->Sumw2();
  // hist24->Sumw2();
  // hist25->Sumw2();
  // hist26->Sumw2();
  // hist27->Sumw2();
  // hist28->Sumw2();
  // hist29->Sumw2();
  // hist210->Sumw2();
  // hist211->Sumw2();
  
  double norm1 = 1;
  double norm2 = 1/3.3033;

  // double norm1 = (double) hist11->GetEntries();
  // double norm2 = (double) hist21->GetEntries();
  hist11->Scale(1/norm1);
  hist12->Scale(1/norm1);
  hist13->Scale(1/norm1);
  hist14->Scale(1/norm1);
  hist15->Scale(1/norm1);
  hist16->Scale(1/norm1);
  hist17->Scale(1/norm1);
  hist18->Scale(1/norm1);
  hist19->Scale(1/norm1);
  hist110->Scale(1/norm1);
  hist111->Scale(1/norm1);
  hist21->Scale(1/norm2);
  hist22->Scale(1/norm2);
  hist23->Scale(1/norm2);
  hist24->Scale(1/norm2);
  hist25->Scale(1/norm2);
  hist26->Scale(1/norm2);
  hist27->Scale(1/norm2);
  hist28->Scale(1/norm2);
  hist29->Scale(1/norm2);
  hist210->Scale(1/norm2);
  hist211->Scale(1/norm2);

  TLegend *legend = new TLegend(0.7,0.82,0.9,0.9);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;
  
  // Residual plot
  TH1 *pull1 = NULL;
  TH1 *pull2 = NULL;
  TH1 *main= NULL;
  TAxis *xaxisadd1 = NULL;
  TAxis *yaxisadd1 = NULL;
  TAxis *xaxisadd2 = NULL;
  TAxis *yaxisadd2 = NULL;
  
  //gStyle->SetHistMinimumZero();
  TCanvas *c01 = new TCanvas("c01","",2400,2100);
  c01->cd();
  TPad *p01a = new TPad("p01a","p01a",0,0.10,1,1.0);
  TPad *p01b = new TPad("p01b","p01b",0,0,1,0.19);
  p01a->Draw();
  p01b->Draw();
  p01a->cd();
  p01a->SetBottomMargin(0.2);
  p01a->SetLogy();
  xaxis = hist11->GetXaxis();
  yaxis = hist11->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} (GeV)");
  yaxis->SetTitle("Entries / 30 GeV");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist11->Draw("HIST");
  hist21->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p01a,iPeriod,iPos);
  legend->Draw();
  
  p01b->cd();
  p01b->SetTopMargin(0.1);
  p01b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist11->Clone();
  main = (TH1*)hist21->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("pt_{#gamma} (GeV)");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p01b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");
  c01->Print("p_pt.png");
  c01->Print("p_pt.pdf");
  c01->Print("p_pt.svg");
  c01->Print("p_pt.root");
  //==========================================================
  
    //gStyle->SetHistMinimumZero();
  TCanvas *c02 = new TCanvas("c02","",2400,2100);
  c02->cd();
  TPad *p02a = new TPad("p02a","p02a",0,0.10,1,1.0);
  TPad *p02b = new TPad("p02b","p02b",0,0,1,0.19);
  p02a->Draw();
  p02b->Draw();
  p02a->cd();
  p02a->SetBottomMargin(0.2);
  p02a->SetLogy();
  xaxis = hist12->GetXaxis();
  yaxis = hist12->GetYaxis();
  xaxis->SetTitle("#eta_{#gamma}");
  yaxis->SetTitle("Entries / 0.08");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist12->Draw("HIST");
  hist22->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p02a,iPeriod,iPos);
  legend->Draw();
  
  p02b->cd();
  p02b->SetTopMargin(0.1);
  p02b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist12->Clone();
  main = (TH1*)hist22->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("#eta_{#gamma}");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p02b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");
  c02->Print("p_eta.png");
  c02->Print("p_eta.pdf");
  c02->Print("p_eta.svg");
  c02->Print("p_eta.root");
  //==========================================================
  
     //gStyle->SetHistMinimumZero();
  TCanvas *c03 = new TCanvas("c03","",2400,2100);
  c03->cd();
  TPad *p03a = new TPad("p03a","p03a",0,0.10,1,1.0);
  TPad *p03b = new TPad("p03b","p03b",0,0,1,0.19);
  p03a->Draw();
  p03b->Draw();
  p03a->cd();
  p03a->SetBottomMargin(0.2);
  p03a->SetLogy();
  xaxis = hist13->GetXaxis();
  yaxis = hist13->GetYaxis();
  xaxis->SetTitle("pt_{j} (GeV)");
  yaxis->SetTitle("Entries / 30 GeV");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist13->Draw("HIST");
  hist23->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p03a,iPeriod,iPos);
  legend->Draw();
  
  p03b->cd();
  p03b->SetTopMargin(0.1);
  p03b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist13->Clone();
  main = (TH1*)hist23->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("pt_{j} (GeV)");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p03b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");

  c03->Print("j_pt.png");
  c03->Print("j_pt.pdf");
  c03->Print("j_pt.svg");
  c03->Print("j_pt.root");
  
  //==========================================================
  
      //gStyle->SetHistMinimumZero();
  TCanvas *c04 = new TCanvas("c04","",2400,2100);
  c04->cd();
  TPad *p04a = new TPad("p04a","p04a",0,0.10,1,1.0);
  TPad *p04b = new TPad("p04b","p04b",0,0,1,0.19);
  p04a->Draw();
  p04b->Draw();
  p04a->cd();
  p04a->SetBottomMargin(0.2);
  p04a->SetLogy();
  xaxis = hist14->GetXaxis();
  yaxis = hist14->GetYaxis();
  xaxis->SetTitle("#eta_{j}");
  yaxis->SetTitle("Entries / 0.08");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist14->Draw("HIST");
  hist24->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p04a,iPeriod,iPos);
  legend->Draw();
  
  p04b->cd();
  p04b->SetTopMargin(0.1);
  p04b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist14->Clone();
  main = (TH1*)hist24->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("#eta_{j}");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p04b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");

  c04->Print("j_eta.png");
  c04->Print("j_eta.pdf");
  c04->Print("j_eta.svg");
  c04->Print("j_eta.root");
  
  //==========================================================
  
       //gStyle->SetHistMinimumZero();
  TCanvas *c05 = new TCanvas("c05","",2400,2100);
  c05->cd();
  TPad *p05a = new TPad("p05a","p05a",0,0.10,1,1.0);
  TPad *p05b = new TPad("p05b","p05b",0,0,1,0.19);
  p05a->Draw();
  p05b->Draw();
  p05a->cd();
  p05a->SetBottomMargin(0.2);
  p05a->SetLogy();
  xaxis = hist15->GetXaxis();
  yaxis = hist15->GetYaxis();
  xaxis->SetTitle("E_{j} (GeV)");
  yaxis->SetTitle("Entries / 30 GeV");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist15->Draw("HIST");
  hist25->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p05a,iPeriod,iPos);
  legend->Draw();
  
  p05b->cd();
  p05b->SetTopMargin(0.1);
  p05b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist15->Clone();
  main = (TH1*)hist25->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("E_{j} (GeV)");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p05b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");

  c05->Print("j_e.png");
  c05->Print("j_e.pdf");
  c05->Print("j_e.svg");
  c05->Print("j_e.root");
  
  //==========================================================
  
        //gStyle->SetHistMinimumZero();
  TCanvas *c06 = new TCanvas("c06","",2400,2100);
  c06->cd();
  TPad *p06a = new TPad("p06a","p06a",0,0.10,1,1.0);
  TPad *p06b = new TPad("p06b","p06b",0,0,1,0.19);
  p06a->Draw();
  p06b->Draw();
  p06a->cd();
  p06a->SetBottomMargin(0.2);
  p06a->SetLogy();
  xaxis = hist16->GetXaxis();
  yaxis = hist16->GetYaxis();
  xaxis->SetTitle("SD m_{j} (GeV)");
  yaxis->SetTitle("Entries / 2 GeV");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist16->Draw("HIST");
  hist26->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p06a,iPeriod,iPos);
  legend->Draw();
  
  p06b->cd();
  p06b->SetTopMargin(0.1);
  p06b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist16->Clone();
  main = (TH1*)hist26->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("SD m_{j} (GeV)");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p06b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");
  c06->Print("j_sdm.png");
  c06->Print("j_sdm.pdf");
  c06->Print("j_sdm.svg");
  c06->Print("j_sdm.root");
  
  //==========================================================
  
          //gStyle->SetHistMinimumZero();
  TCanvas *c07 = new TCanvas("c07","",2400,2100);
  c07->cd();
  TPad *p07a = new TPad("p07a","p07a",0,0.10,1,1.0);
  TPad *p07b = new TPad("p07b","p07b",0,0,1,0.19);
  p07a->Draw();
  p07b->Draw();
  p07a->cd();
  p07a->SetBottomMargin(0.2);
  p07a->SetLogy();
  xaxis = hist17->GetXaxis();
  yaxis = hist17->GetYaxis();
  xaxis->SetTitle("#tau21_{j}");
  yaxis->SetTitle("Entries / 0.02");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist17->Draw("HIST");
  hist27->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p07a,iPeriod,iPos);
  legend->Draw();
  
  p07b->cd();
  p07b->SetTopMargin(0.1);
  p07b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist17->Clone();
  main = (TH1*)hist27->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("#tau21_{j}");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p07b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");
  
  c07->Print("j_tau21.png");
  c07->Print("j_tau21.pdf");
  c07->Print("j_tau21.svg");
  c07->Print("j_tau21.root");
  
  //==========================================================
  
            //gStyle->SetHistMinimumZero();
  TCanvas *c08 = new TCanvas("c08","",2400,2100);
  c08->cd();
  TPad *p08a = new TPad("p08a","p08a",0,0.10,1,1.0);
  TPad *p08b = new TPad("p08b","p08b",0,0,1,0.19);
  p08a->Draw();
  p08b->Draw();
  p08a->cd();
  p08a->SetBottomMargin(0.2);
  p08a->SetLogy();
  xaxis = hist18->GetXaxis();
  yaxis = hist18->GetYaxis();
  xaxis->SetTitle("cos(#theta*)");
  yaxis->SetTitle("Entries / 0.02");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist18->Draw("HIST");
  hist28->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p08a,iPeriod,iPos);
  legend->Draw();
  
  p08b->cd();
  p08b->SetTopMargin(0.1);
  p08b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist18->Clone();
  main = (TH1*)hist28->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("cos(#theta*)");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p08b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");

  c08->Print("s_cos.png");
  c08->Print("s_cos.pdf");
  c08->Print("s_cos.svg");
  c08->Print("s_cos.root");
  
  //==========================================================
  
    //==========================================================
  
            //gStyle->SetHistMinimumZero();
  TCanvas *c09 = new TCanvas("c09","",2400,2100);
  c09->cd();
  TPad *p09a = new TPad("p09a","p09a",0,0.10,1,1.0);
  TPad *p09b = new TPad("p09b","p09b",0,0,1,0.19);
  p09a->Draw();
  p09b->Draw();
  p09a->cd();
  p09a->SetBottomMargin(0.2);
  p09a->SetLogy();
  xaxis = hist19->GetXaxis();
  yaxis = hist19->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} / M_{j#gamma}");
  yaxis->SetTitle("Entries / 0.04");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist19->Draw("HIST");
  hist29->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p09a,iPeriod,iPos);
  legend->Draw();
  
  p09b->cd();
  p09b->SetTopMargin(0.1);
  p09b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist19->Clone();
  main = (TH1*)hist29->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("pt_{#gamma} / M_{j#gamma}");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p09b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");
  
  c09->Print("s_ptm.png");
  c09->Print("s_ptm.pdf");
  c09->Print("s_ptm.svg");
  c09->Print("s_ptm.root");
  
  //==========================================================
  
  
            //gStyle->SetHistMinimumZero();
  TCanvas *c10 = new TCanvas("c10","",2400,2100);
  c10->cd();
  TPad *p10a = new TPad("p10a","p10a",0,0.10,1,1.0);
  TPad *p10b = new TPad("p10b","p10b",0,0,1,0.19);
  p10a->Draw();
  p10b->Draw();
  p10a->cd();
  p10a->SetBottomMargin(0.2);
  p10a->SetLogy();
  xaxis = hist110->GetXaxis();
  yaxis = hist110->GetYaxis();
  xaxis->SetTitle("M_{j#gamma} (GeV)");
  yaxis->SetTitle("Entries / 40 GeV");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist110->Draw("HIST");
  hist210->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p10a,iPeriod,iPos);
  legend->Draw();
  
  p10b->cd();
  p10b->SetTopMargin(0.1);
  p10b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist110->Clone();
  main = (TH1*)hist210->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("M_{j#gamma} (GeV)");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p10b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");
 
  c10->Print("s_M.png");
  c10->Print("s_M.pdf");
  c10->Print("s_M.svg");
  c10->Print("s_M.root");
  
  //==========================================================
  
  
            //gStyle->SetHistMinimumZero();
  TCanvas *c11 = new TCanvas("c11","",2400,2100);
  c11->cd();
  TPad *p11a = new TPad("p11a","p11a",0,0.10,1,1.0);
  TPad *p11b = new TPad("p11b","p11b",0,0,1,0.19);
  p11a->Draw();
  p11b->Draw();
  p11a->cd();
  p11a->SetBottomMargin(0.2);
  p11a->SetLogy();
  xaxis = hist111->GetXaxis();
  yaxis = hist111->GetYaxis();
  xaxis->SetTitle("#gamma Data SR");
  yaxis->SetTitle("Entries / 0.04");
  xaxis->SetTitleOffset(1.1);
  yaxis->SetTitleOffset(1.35);
  yaxis->SetRangeUser(0.1,10000);
  hist111->Draw("HIST");
  hist211->Draw("HISTSAME");
  legend->Clear();
  legend->AddEntry(hist11,"Data SR","f");
  legend->AddEntry(hist21,"MC SR","f");
  legend->Draw();
  CMS_lumi(p11a,iPeriod,iPos);
  legend->Draw();
  
  p11b->cd();
  p11b->SetTopMargin(0.1);
  p11b->SetBottomMargin(0.3);
  pull1 = (TH1*)hist111->Clone();
  main = (TH1*)hist211->Clone();
  pull1->Divide(main);
  xaxisadd1 = pull1->GetXaxis();
  yaxisadd1 = pull1->GetYaxis();
  xaxisadd1->SetTitle("#gamma Data SR");
  xaxisadd1->SetTitleOffset(1.25);
  yaxisadd1->SetTitle("Data SR / MC SR");
  yaxisadd1->SetTitleOffset(0.4);
  yaxisadd1->SetRangeUser(0,2);
  xaxisadd1->SetLabelSize(0.1);
  xaxisadd1->SetTitleSize(0.1);
  yaxisadd1->SetLabelSize(0.1);
  yaxisadd1->SetNdivisions(5);
  yaxisadd1->SetTitleSize(0.1);
  p11b->SetGrid();
  //pull1->SetFillColor(kViolet);
  //pull1->SetLineColor(kViolet);
  //pull1->SetFillColorAlpha(kViolet, 0.35);
  pull1->Draw("PE1");
  //pull1->Draw("BAR HIST");
  c11->Print("p_mva.png");
  c11->Print("p_mva.pdf");
  c11->Print("p_mva.svg");
  c11->Print("p_mva.root");
  
  //==========================================================

  cout<<"Sum weights: "<<sumW<<endl;
}
