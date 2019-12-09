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
#include "../Utils/interface/ConfParse.hh"             // input conf file parser
#include "../Utils/interface/CSample.hh"      // helper class to handle samples
#include <algorithm>
#include <map>
#endif
#define mode 2 //1 for MC W band, 2 for sideband

double cut_peta(TTree* sigTree, TTree* bkgTree, int sigmass){
  // Data structures to store info from produced flat ntuples
  float photon_pt;
  float photon_eta;
  float photon_phi;
  float photon_e;
  float photon_mvaval;
  float photon_mvacat;
  float ak8puppijet_pt;
  float ak8puppijet_eta;
  float ak8puppijet_phi;
  float ak8puppijet_e;
  float ak8puppijet_masssoftdropcorr;
  float ak8puppijet_tau21;
  float sys_seperation;
  float sys_costhetastar;
  float sys_ptoverm;
  float sys_invmass;
  float sys_BDT;
  float xsec_weight;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<int> sig_N, bkg_N;
  std::vector<float> recordarr;

  // Access Event Tree
  sigTree->SetBranchAddress("photon_pt", &photon_pt);                       
  sigTree->SetBranchAddress("photon_eta", &photon_eta);                       
  sigTree->SetBranchAddress("photon_phi", &photon_phi);                         
  sigTree->SetBranchAddress("photon_e", &photon_e);                                             
  sigTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  sigTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  sigTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  sigTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  sigTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  sigTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  sigTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  sigTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  sigTree->SetBranchAddress("sys_invmass", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 2; cut+=0.05){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) < cut){
	N++;
      }
    }
    sig_N.push_back(N);
    recordarr.push_back(record);
  }

  // Access Event Tree
  bkgTree->SetBranchAddress("photon_pt", &photon_pt);                       
  bkgTree->SetBranchAddress("photon_eta", &photon_eta);                       
  bkgTree->SetBranchAddress("photon_phi", &photon_phi);                         
  bkgTree->SetBranchAddress("photon_e", &photon_e);                                             
  bkgTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  bkgTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  bkgTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  bkgTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  bkgTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  bkgTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  bkgTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  bkgTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  bkgTree->SetBranchAddress("sys_invmass", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  bkgTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 2; cut+=0.05){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) < cut){
#if mode == 1
	N = N+1*xsec_weight;
#else
	N++;
#endif
      }
    }
    bkg_N.push_back(N);
    recordarr.push_back(record);
  }

  TString Graphname ="|#eta| cut on photon";

  //Non stacked plots
  TLegend *legend1 = new TLegend(0.15,0.75,0.3,0.85);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
    if(bkg_N[i] == 0)
      SBratio.push_back(0);
    else
      SBratio.push_back(float(sig_N[i])/sqrt(float(bkg_N[i])));
  }
  int dim = SBratio.size();

  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(2);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,1200,900);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  xaxis->SetRangeUser(0,2);
  yaxis->SetRangeUser(0,250);
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  c01->Print(name+"_peta.png");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

double cut_petajeta(TTree* sigTree, TTree* bkgTree, int sigmass, double peta_cut){
  // Data structures to store info from produced flat ntuples
  float photon_pt;
  float photon_eta;
  float photon_phi;
  float photon_e;
  float photon_mvaval;
  float photon_mvacat;
  float ak8puppijet_pt;
  float ak8puppijet_eta;
  float ak8puppijet_phi;
  float ak8puppijet_e;
  float ak8puppijet_masssoftdropcorr;
  float ak8puppijet_tau21;
  float sys_seperation;
  float sys_costhetastar;
  float sys_ptoverm;
  float sys_invmass;
  float sys_BDT;
  float xsec_weight;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<int> sig_N, bkg_N;
  std::vector<float> recordarr;

  // Access Event Tree
  sigTree->SetBranchAddress("photon_pt", &photon_pt);                       
  sigTree->SetBranchAddress("photon_eta", &photon_eta);                       
  sigTree->SetBranchAddress("photon_phi", &photon_phi);                         
  sigTree->SetBranchAddress("photon_e", &photon_e);                                             
  sigTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  sigTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  sigTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  sigTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  sigTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  sigTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  sigTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  sigTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  sigTree->SetBranchAddress("sys_invmass", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 2; cut+=0.05){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) < cut){
	N++;
      }
    }
    sig_N.push_back(N);
    recordarr.push_back(record);
  }

  // Access Event Tree
  bkgTree->SetBranchAddress("photon_pt", &photon_pt);                       
  bkgTree->SetBranchAddress("photon_eta", &photon_eta);                       
  bkgTree->SetBranchAddress("photon_phi", &photon_phi);                         
  bkgTree->SetBranchAddress("photon_e", &photon_e);                                             
  bkgTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  bkgTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  bkgTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  bkgTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  bkgTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  bkgTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  bkgTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  bkgTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  bkgTree->SetBranchAddress("sys_invmass", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  bkgTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 2; cut+=0.05){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) < cut){
#if mode == 1
	N = N+1*xsec_weight;
#else
	N++;
#endif
      }
    }
    bkg_N.push_back(N);
    recordarr.push_back(record);
  }

  TString Graphname ="|#eta| cut on jet";

  //Non stacked plots
  TLegend *legend1 = new TLegend(0.15,0.75,0.3,0.85);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
    if(bkg_N[i] == 0)
      SBratio.push_back(0);
    else
      SBratio.push_back(float(sig_N[i])/sqrt(float(bkg_N[i])));
  }
  int dim = SBratio.size();
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(2);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,1200,900);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  xaxis->SetRangeUser(0,2);
  yaxis->SetRangeUser(0,250);
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  c01->Print(name+"_petajeta.png");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

double cut_petajetaptm(TTree* sigTree, TTree* bkgTree, int sigmass, double peta_cut, double jeta_cut){
  // Data structures to store info from produced flat ntuples
  float photon_pt;
  float photon_eta;
  float photon_phi;
  float photon_e;
  float photon_mvaval;
  float photon_mvacat;
  float ak8puppijet_pt;
  float ak8puppijet_eta;
  float ak8puppijet_phi;
  float ak8puppijet_e;
  float ak8puppijet_masssoftdropcorr;
  float ak8puppijet_tau21;
  float sys_seperation;
  float sys_costhetastar;
  float sys_ptoverm;
  float sys_invmass;
  float sys_BDT;
  float xsec_weight;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<int> sig_N, bkg_N;
  std::vector<float> recordarr;

  // Access Event Tree
  sigTree->SetBranchAddress("photon_pt", &photon_pt);                       
  sigTree->SetBranchAddress("photon_eta", &photon_eta);                       
  sigTree->SetBranchAddress("photon_phi", &photon_phi);                         
  sigTree->SetBranchAddress("photon_e", &photon_e);                                             
  sigTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  sigTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  sigTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  sigTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  sigTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  sigTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  sigTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  sigTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  sigTree->SetBranchAddress("sys_invmass", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 1; cut+=0.01){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm > cut){
	N++;
      }
    }
    sig_N.push_back(N);
    recordarr.push_back(record);
  }

  // Access Event Tree
  bkgTree->SetBranchAddress("photon_pt", &photon_pt);                       
  bkgTree->SetBranchAddress("photon_eta", &photon_eta);                       
  bkgTree->SetBranchAddress("photon_phi", &photon_phi);                         
  bkgTree->SetBranchAddress("photon_e", &photon_e);                                             
  bkgTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  bkgTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  bkgTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  bkgTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  bkgTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  bkgTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  bkgTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  bkgTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  bkgTree->SetBranchAddress("sys_invmass", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  bkgTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 1; cut+=0.01){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm > cut){
#if mode == 1
	N = N+1*xsec_weight;
#else
	N++;
#endif
      }
    }
    bkg_N.push_back(N);
    recordarr.push_back(record);
  }

  TString Graphname ="pt/M cut on photon";

  //Non stacked plots
  TLegend *legend1 = new TLegend(0.15,0.75,0.3,0.85);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
    if(bkg_N[i] == 0)
      SBratio.push_back(0);
    else
      SBratio.push_back(float(sig_N[i])/sqrt(float(bkg_N[i])));
  }
  int dim = SBratio.size();
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(2);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,1200,900);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  xaxis->SetRangeUser(0,1);
  yaxis->SetRangeUser(0,250);
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  c01->Print(name+"_petajetaptm.png");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

double cut_petajetaptmcos(TTree* sigTree, TTree* bkgTree, int sigmass, double peta_cut, double jeta_cut, double ptm_cut){
  // Data structures to store info from produced flat ntuples
  float photon_pt;
  float photon_eta;
  float photon_phi;
  float photon_e;
  float photon_mvaval;
  float photon_mvacat;
  float ak8puppijet_pt;
  float ak8puppijet_eta;
  float ak8puppijet_phi;
  float ak8puppijet_e;
  float ak8puppijet_masssoftdropcorr;
  float ak8puppijet_tau21;
  float sys_seperation;
  float sys_costhetastar;
  float sys_ptoverm;
  float sys_invmass;
  float sys_BDT;
  float xsec_weight;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<int> sig_N, bkg_N;
  std::vector<float> recordarr;

  // Access Event Tree
  sigTree->SetBranchAddress("photon_pt", &photon_pt);                       
  sigTree->SetBranchAddress("photon_eta", &photon_eta);                       
  sigTree->SetBranchAddress("photon_phi", &photon_phi);                         
  sigTree->SetBranchAddress("photon_e", &photon_e);                                             
  sigTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  sigTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  sigTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  sigTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  sigTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  sigTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  sigTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  sigTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  sigTree->SetBranchAddress("sys_invmass", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 1; cut+=0.01){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm < ptm_cut) continue;
      if(sys_costhetastar < cut){
	N++;
      }
    }
    sig_N.push_back(N);
    recordarr.push_back(record);
  }

  // Access Event Tree
  bkgTree->SetBranchAddress("photon_pt", &photon_pt);                       
  bkgTree->SetBranchAddress("photon_eta", &photon_eta);                       
  bkgTree->SetBranchAddress("photon_phi", &photon_phi);                         
  bkgTree->SetBranchAddress("photon_e", &photon_e);                                             
  bkgTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  bkgTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  bkgTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  bkgTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  bkgTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  bkgTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  bkgTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  bkgTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  bkgTree->SetBranchAddress("sys_invmass", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  bkgTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 1; cut+=0.01){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm < ptm_cut) continue;
      if(sys_costhetastar < cut){
#if mode == 1
	N = N+1*xsec_weight;
#else
	N++;
#endif
      }
    }
    bkg_N.push_back(N);
    recordarr.push_back(record);
  }

  TString Graphname ="cos(#theta*) cut on photon";

  //Non stacked plots
  TLegend *legend1 = new TLegend(0.15,0.75,0.3,0.85);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
    if(bkg_N[i] == 0)
      SBratio.push_back(0);
    else
      SBratio.push_back(float(sig_N[i])/sqrt(float(bkg_N[i])));
  }
  int dim = SBratio.size();
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(2);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,1200,900);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  xaxis->SetRangeUser(0,1);
  yaxis->SetRangeUser(0,250);
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  c01->Print(name+"_petajetaptmcos.png");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

double cut_petajetaptmcostau(TTree* sigTree, TTree* bkgTree, int sigmass, double peta_cut, double jeta_cut, double ptm_cut, double cos_cut){
  // Data structures to store info from produced flat ntuples
  float photon_pt;
  float photon_eta;
  float photon_phi;
  float photon_e;
  float photon_mvaval;
  float photon_mvacat;
  float ak8puppijet_pt;
  float ak8puppijet_eta;
  float ak8puppijet_phi;
  float ak8puppijet_e;
  float ak8puppijet_masssoftdropcorr;
  float ak8puppijet_tau21;
  float sys_seperation;
  float sys_costhetastar;
  float sys_ptoverm;
  float sys_invmass;
  float sys_BDT;
  float xsec_weight;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<int> sig_N, bkg_N;
  std::vector<float> recordarr;

  // Access Event Tree
  sigTree->SetBranchAddress("photon_pt", &photon_pt);                       
  sigTree->SetBranchAddress("photon_eta", &photon_eta);                       
  sigTree->SetBranchAddress("photon_phi", &photon_phi);                         
  sigTree->SetBranchAddress("photon_e", &photon_e);                                             
  sigTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  sigTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  sigTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  sigTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  sigTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  sigTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  sigTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  sigTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  sigTree->SetBranchAddress("sys_invmass", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 1; cut+=0.01){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_costhetastar > cos_cut) continue;
      if(ak8puppijet_tau21 < cut){
	N++;
      }
    }
    sig_N.push_back(N);
    recordarr.push_back(record);
  }

  // Access Event Tree
  bkgTree->SetBranchAddress("photon_pt", &photon_pt);                       
  bkgTree->SetBranchAddress("photon_eta", &photon_eta);                       
  bkgTree->SetBranchAddress("photon_phi", &photon_phi);                         
  bkgTree->SetBranchAddress("photon_e", &photon_e);                                             
  bkgTree->SetBranchAddress("ak8puppijet_pt", &ak8puppijet_pt);                           
  bkgTree->SetBranchAddress("ak8puppijet_eta", &ak8puppijet_eta);                         
  bkgTree->SetBranchAddress("ak8puppijet_phi", &ak8puppijet_phi);                          
  bkgTree->SetBranchAddress("ak8puppijet_e", &ak8puppijet_e);                               
  bkgTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &ak8puppijet_masssoftdropcorr);
  bkgTree->SetBranchAddress("ak8puppijet_tau21", &ak8puppijet_tau21);                     
  bkgTree->SetBranchAddress("sys_costhetastar", &sys_costhetastar);                       
  bkgTree->SetBranchAddress("sys_ptoverm", &sys_ptoverm);                                
  bkgTree->SetBranchAddress("sys_invmass", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  bkgTree->SetBranchAddress("xsec_weight", &xsec_weight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 1; cut+=0.01){
    int N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 65 || ak8puppijet_masssoftdropcorr > 105) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.05 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm < ptm_cut) continue;
      if(sys_costhetastar > cos_cut) continue;
      if(ak8puppijet_tau21 < cut){
#if mode == 1
	N = N+1*xsec_weight;
#else
	N++;
#endif
      }
    }
    bkg_N.push_back(N);
    recordarr.push_back(record);
  }

  TString Graphname ="tau21 cut on jet";

  //Non stacked plots
  TLegend *legend1 = new TLegend(0.15,0.75,0.3,0.85);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
    if(bkg_N[i] == 0)
      SBratio.push_back(0);
    else
      SBratio.push_back(float(sig_N[i])/sqrt(float(bkg_N[i])));
  }
  int dim = SBratio.size();
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(2);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,1200,900);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  xaxis->SetRangeUser(0,1);
  yaxis->SetRangeUser(0,250);
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  c01->Print(name+"_petajetaptmcostau21.png");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

void significance_count1D(){

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  int sigmass = 2000;
  TString sig_sample = "/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC"+std::to_string(sigmass)+"_WGamma_full_full.root";
  //TString bkg_sample = "/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/BackgroundCombinedMC_WGamma_full_full_weightedTo41p54_fitData.root";
  TString bkg_sample = "/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017_WGamma_full_full.root";

  TFile *infile_1 = TFile::Open(sig_sample);
  TTree* sigTree = (TTree*)infile_1->Get("Events");
  TFile *infile_2 = TFile::Open(bkg_sample);
  TTree* bkgTree = (TTree*)infile_2->Get("Events");

  double peta_cut = cut_peta(sigTree,bkgTree,sigmass);
  cout<<"peta_cut: "<<peta_cut<<endl;
  
  double jeta_cut = cut_petajeta(sigTree,bkgTree,sigmass,peta_cut);
  cout<<"jeta_cut: "<<jeta_cut<<endl;
  
  double ptm_cut = cut_petajetaptm(sigTree,bkgTree,sigmass,peta_cut,jeta_cut);
  cout<<"ptm_cut: "<<ptm_cut<<endl;
  
  double cos_cut = cut_petajetaptmcos(sigTree,bkgTree,sigmass,peta_cut,jeta_cut,ptm_cut);
  cout<<"cos_cut: "<<cos_cut<<endl;
  
  //double tau21_cut = cut_petajetaptmcostau(sigTree,bkgTree,sigmass,peta_cut,jeta_cut,ptm_cut,cos_cut);
  //cout<<"tau_cut: "<<tau21_cut<<endl;

  /*
  TGraph* gr1 = new TGraph(dim,&recordarr[0],&Sigeff[0]);
  gr1->SetTitle(Graphname+" efficiency");
  gr1->SetLineColor(6);
  gr1->SetLineWidth(4);
  //gr1->SetMarkerColor(4);
  //gr1->SetMarkerSize(1.5);

  TGraph* gr2 = new TGraph(dim,&recordarr[0],&Bkgeff[0]);
  //gr1->SetTitle(Graphname+" efficiency");
  gr2->SetLineColor(4);
  gr2->SetLineWidth(4);
  //gr2->SetMarkerColor(4);
  //gr2->SetMarkerSize(1.5);
  TCanvas *c02 = new TCanvas("c02",Graphname+" Efficiency",1600,900);
  xaxis = gr1->GetXaxis();
  yaxis = gr1->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("Eff");
  yaxis->SetTitleOffset(1.1);
  xaxis->SetTitleOffset(1.3);
  xaxis->SetRangeUser(0,3);
  yaxis->SetRangeUser(0,1);
  //yaxis->SetRangeUser(0.5,16000000);
  //c01->SetLogy();
  c02->cd();
  c02->SetGrid();
  cout<<"OK"<<endl;
  gr1->Draw("ACY+");
  gr2->Draw("SAME");
  legend1->Clear();
  legend1->AddEntry(gr1,"signal efficiency");
  legend1->AddEntry(gr2,"background efficiency");
  //legend1->Draw();
  c02->Print("1600_invptmcosphi_eff.png");
  */

}
