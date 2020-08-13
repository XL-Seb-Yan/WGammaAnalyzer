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

void compare_plot_N()
{
  gROOT->SetBatch(1);
  lumi_13TeV = "35.92 fb^{-1} + 41.53 fb^{-1} + 59.74 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.2;
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
  
  std::vector<TH1*> p_pt_hist;
  std::vector<TH1*> p_eta_hist;
  std::vector<TH1*> j_pt_hist;
  std::vector<TH1*> j_eta_hist;
  std::vector<TH1*> j_E_hist;
  std::vector<TH1*> j_tau21_hist;
  std::vector<TH1*> j_sdmass_hist;
  std::vector<TH1*> s_ptoverm_hist;
  std::vector<TH1*> s_costheta_hist;
  std::vector<TH1*> s_M_hist;
  
  std::vector<TString> variable_names, binsize;
  variable_names.push_back("pt_{#gamma}"); binsize.push_back("30 GeV");
  variable_names.push_back("eta_{#gamma}"); binsize.push_back("0.08");
  variable_names.push_back("pt_{jet}"); binsize.push_back("30 GeV");
  variable_names.push_back("eta_{jet}"); binsize.push_back("0.08");
  variable_names.push_back("E_{jet}"); binsize.push_back("50 GeV");
  variable_names.push_back("mass_{jet}"); binsize.push_back("2 GeV");
  variable_names.push_back("tau21_{jet}"); binsize.push_back("0.02");
  variable_names.push_back("pt/M"); binsize.push_back("0.02");
  variable_names.push_back("cos(#theta*)"); binsize.push_back("0.04");
  variable_names.push_back("invariant mass"); binsize.push_back("70 GeV");
  
  std::vector<TString> sample_names;
  sample_names.push_back("../selection/SignalMC7000_N_0-4");
  sample_names.push_back("../selection/SignalMC7000_N_5-9");
  sample_names.push_back("../selection/SignalMC7000_N_10-14");
  sample_names.push_back("../selection/SignalMC7000_N_15-19");
  
  for(int ifile=0; ifile<sample_names.size();ifile++){
    TH1 *hist1 = new TH1F("11","pt_{#gamma}",100,0,6000);
    TH1 *hist2 = new TH1F("12","eta_{#gamma}",50,-2,2);
    TH1 *hist3 = new TH1F("13","pt_{j}",100,0,10000);
    TH1 *hist4 = new TH1F("14","eta_{j}",50,-2,2);
    TH1 *hist5 = new TH1F("15","E_{j}",100,0,6000);
    TH1 *hist6 = new TH1F("16","masssoftdrop_{j}",40,40,120);
    TH1 *hist7 = new TH1F("17","tau21_{j}",50,0,1);
    TH1 *hist8 = new TH1F("18","cos(#theta*)_{p}",50,0,1);
    TH1 *hist9 = new TH1F("19","pt/M",50,0,2);
    TH1 *hist10 = new TH1F("110","invariant mass",100,0,10000);
    
    // Open input file
    Float_t p_pt, p_eta, p_phi, p_e, p_mva, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight, x_puweight, x_kfactor;
    
    TFile *input = TFile::Open(sample_names.at(ifile)+"_postproc_WGamma17_full_full_presel_jmcorr_May22.root");
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
    theTree->SetBranchAddress("xsec_weight", &x_weight);
    theTree->SetBranchAddress("xsec_puweight", &x_puweight);
    theTree->SetBranchAddress("xsec_kfactor", &x_kfactor);
    
    for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
      theTree->GetEntry(ievt);
      x_puweight = 1;
      x_kfactor = 1;
      hist1->Fill(p_pt,x_puweight*x_weight*x_kfactor);
      hist2->Fill(p_eta,x_puweight*x_weight*x_kfactor);
      hist3->Fill(j_pt,x_puweight*x_weight*x_kfactor);
      hist4->Fill(j_eta,x_puweight*x_weight*x_kfactor);
      hist5->Fill(j_e,x_puweight*x_weight*x_kfactor);
      hist6->Fill(j_mass,x_puweight*x_weight*x_kfactor);
      hist7->Fill(j_tau21,x_puweight*x_weight*x_kfactor);
      hist8->Fill(s_cos,x_puweight*x_weight*x_kfactor);
      hist9->Fill(s_ptm,x_puweight*x_weight*x_kfactor);
      hist10->Fill(s_mass,x_puweight*x_weight*x_kfactor);
    }
    
    int color = ifile+1;
    hist1->SetLineColor(color);
    hist2->SetLineColor(color);
    hist3->SetLineColor(color);
    hist4->SetLineColor(color);
    hist5->SetLineColor(color);
    hist6->SetLineColor(color);
    hist7->SetLineColor(color);
    hist8->SetLineColor(color);
    hist9->SetLineColor(color);
    hist10->SetLineColor(color);
    
    hist1->SetLineWidth(3);
    hist2->SetLineWidth(3);
    hist3->SetLineWidth(3);
    hist4->SetLineWidth(3);
    hist5->SetLineWidth(3);
    hist6->SetLineWidth(3);
    hist7->SetLineWidth(3);
    hist8->SetLineWidth(3);
    hist9->SetLineWidth(3);
    hist10->SetLineWidth(3);
    
    double norm = (double) hist1->GetEntries();
    hist1->Scale(1/norm);
    hist2->Scale(1/norm);
    hist3->Scale(1/norm);
    hist4->Scale(1/norm);
    hist5->Scale(1/norm);
    hist6->Scale(1/norm);
    hist7->Scale(1/norm);
    hist8->Scale(1/norm);
    hist9->Scale(1/norm);
    hist10->Scale(1/norm);
    
    p_pt_hist.push_back(hist1);
    p_eta_hist.push_back(hist2);
    j_pt_hist.push_back(hist3);
    j_eta_hist.push_back(hist4);
    j_E_hist.push_back(hist5);
    j_tau21_hist.push_back(hist6);
    j_sdmass_hist.push_back(hist7);
    s_ptoverm_hist.push_back(hist8);
    s_costheta_hist.push_back(hist9);
    s_M_hist.push_back(hist10);
}
  std::vector<std::vector<TH1*>> all_hist;
  all_hist.push_back(p_pt_hist);
  all_hist.push_back(p_eta_hist);
  all_hist.push_back(j_pt_hist);
  all_hist.push_back(j_eta_hist);
  all_hist.push_back(j_E_hist);
  all_hist.push_back(j_tau21_hist);
  all_hist.push_back(j_sdmass_hist);
  all_hist.push_back(s_ptoverm_hist);
  all_hist.push_back(s_costheta_hist);
  all_hist.push_back(s_M_hist);
  
  TLegend *legend = new TLegend(0.60,0.78,0.9,0.9);
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
  
  for(int ivar = 0; ivar<variable_names.size(); ivar++){
    TCanvas *c01 = new TCanvas("c01","",1200,900);
    c01->cd();
    c01->cd();
    //c01>SetBottomMargin(0.2);
    c01->SetLogy();
    xaxis = all_hist.at(ivar).at(0)->GetXaxis();
    yaxis = all_hist.at(ivar).at(0)->GetYaxis();
    xaxis->SetTitle(variable_names.at(ivar));
    yaxis->SetTitle("Entries / "+binsize.at(ivar));
    xaxis->SetTitleOffset(1.1);
    yaxis->SetTitleOffset(1.35);
    yaxis->SetRangeUser(0.0001,1);
    legend->Clear();
    all_hist.at(ivar).at(0)->Draw("HIST");
    for(int isample=0; isample<all_hist.at(ivar).size();isample++){
        all_hist.at(ivar).at(isample)->Draw("HIST SAME");
        legend->AddEntry(all_hist.at(ivar).at(isample),sample_names.at(isample),"f");
    }
    legend->Draw();
    CMS_lumi(c01,iPeriod,iPos);
    legend->Draw();
    c01->Print(variable_names.at(ivar)+".png");
    c01->Print(variable_names.at(ivar)+".pdf");
    c01->Print(variable_names.at(ivar)+".svg");
    //==========================================================
  }

}
