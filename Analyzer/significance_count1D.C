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
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#include <algorithm>
#include <map>
#endif
#define mode 2 //1 for MC W band, 2 for sideband

double cut_jmass(TTree* sigTree, TTree* bkgTree, int sigmass, int color, int iPeriod,int iPos, bool plot_CMS){
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
  float xsec_weight = 1;
  float xsec_kfactor = 1;
  float xsec_puweight = 1;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<std::vector<float>> sig_N, bkg_N;
  std::vector<std::vector<float>> highlimit;
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
  sigTree->SetBranchAddress("m", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cutlow = 60; cutlow < 75; cutlow+=1){
	cout<<"Running through: "<<cutlow<<endl;
	recordarr.push_back(cutlow);
	std::vector<float> sig_N_tmp;
	std::vector<float> highlimit_tmp;
	for(float cuthigh = 85; cuthigh < 105; cuthigh+=1){
		float N = 0;
        for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
          // Get Events
          sigTree->GetEntry(ientry);
          if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
	      if(ak8puppijet_masssoftdropcorr > cutlow && ak8puppijet_masssoftdropcorr < cuthigh)
	          N+=xsec_weight*xsec_kfactor*xsec_puweight;
        }
        sig_N_tmp.push_back(N);
        highlimit_tmp.push_back(cuthigh);
      }
	sig_N.push_back(sig_N_tmp);
	highlimit.push_back(highlimit_tmp);
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
  bkgTree->SetBranchAddress("m", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nbkg = bkgTree->GetEntries();
  
  for(float cutlow = 60; cutlow < 75; cutlow+=1){
	cout<<"Running through: "<<cutlow<<endl;
	recordarr.push_back(cutlow);
	std::vector<float> bkg_N_tmp;
	for(float cuthigh = 85; cuthigh < 105; cuthigh+=1){
		float N = 0;
        for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
          // Get Events
          bkgTree->GetEntry(ientry);
          if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
	      if(ak8puppijet_masssoftdropcorr > cutlow && ak8puppijet_masssoftdropcorr < cuthigh){
#if mode == 1
	        N+=xsec_weight*xsec_kfactor*xsec_puweight;
#else
	        N++;
#endif
          }
		}
		bkg_N_tmp.push_back(N);
	}
    bkg_N.push_back(bkg_N_tmp);
  }

  TString Graphname ="mass cut on jet";

  //Non stacked plots
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  TGraph2D* gr = new TGraph2D();
  for(int i=0; i<sig_N.size(); i++){
	  for(int j=0; j<sig_N.at(0).size(); j++){
	    float sbratio = 0;
	    if(bkg_N.at(i).at(j) != 0)
          sbratio = sig_N.at(i).at(j)/sqrt(bkg_N.at(i).at(j));
	    gr->SetPoint(i*sig_N.at(0).size()+j,recordarr.at(i),highlimit.at(i).at(j),sbratio);
		cout<<recordarr.at(i)<<" "<<highlimit.at(i).at(j)<<" "<<sbratio<<endl;
	  }
  }

  gr->SetTitle(Graphname);
  // gr->SetLineColor(color);
  // gr->SetLineWidth(4);
  // gr->SetMarkerColor(4);
  // gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,2400,1800);
  c01->SetLeftMargin(0.15);
  c01->SetRightMargin(0.15);
  c01->SetBottomMargin(0.15);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle("Lower boundary (GeV)");
  yaxis->SetTitle("Upper boundary (GeV)");
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("COLZ");
  TString name = std::to_string(sigmass);
  TLegend *legend1 = new TLegend(0.80,0.86,0.899,0.899);
  legend1->AddEntry(gr,"M-"+name,"l");
  legend1->Draw();
  if(plot_CMS) CMS_lumi(c01,iPeriod,iPos);
  c01->Print(name+"_jm.png");
  c01->Print(name+"_jm.pdf");
  c01->Print(name+"_jm.root");
  c01->Print(name+"_jm.svg");
  return 0;
}

double cut_peta(TTree* sigTree, TTree* bkgTree, int sigmass, int color, int iPeriod,int iPos, bool plot_CMS){
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
  float xsec_weight = 1;
  float xsec_kfactor = 1;
  float xsec_puweight = 1;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<float> sig_N, bkg_N;
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
  sigTree->SetBranchAddress("m", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 2; cut+=0.05){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      if(abs(photon_eta) < cut){
	    N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  bkgTree->SetBranchAddress("m", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 2; cut+=0.05){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      if(abs(photon_eta) < cut){
#if mode == 1
	N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
	float sbratio = 0;
	if(bkg_N[i] != 0)
      sbratio = float(sig_N[i])/sqrt(float(bkg_N[i]));
    SBratio.push_back(sbratio);
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<" "<<sbratio<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
  }
  int dim = SBratio.size();
  for(int i=0; i<dim; i++){
	  SBratio.at(i) = SBratio.at(i)/SBratio.back();
  }

  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(color);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,2400,1800);
  c01->SetBottomMargin(0.13);
  c01->SetLeftMargin(0.13);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.15);
  xaxis->SetTitleOffset(1.15);
  xaxis->SetRangeUser(0,2);
  yaxis->SetRangeUser(0,1.5);
  //c01->SetLogy();
  c01->cd();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  TLegend *legend1 = new TLegend(0.80,0.86,0.899,0.899);
  legend1->AddEntry(gr,"M-"+name,"l");
  legend1->Draw();
  if(plot_CMS) CMS_lumi(c01,iPeriod,iPos);
  c01->Print(name+"_peta.png");
  c01->Print(name+"_peta.pdf");
  c01->Print(name+"_peta.root");
  c01->Print(name+"_peta.svg");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

double cut_petajeta(TTree* sigTree, TTree* bkgTree, int sigmass, double peta_cut, int color, int iPeriod,int iPos, bool plot_CMS){
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
  float xsec_weight = 1;
  float xsec_kfactor = 1;
  float xsec_puweight = 1;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<float> sig_N, bkg_N;
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
  sigTree->SetBranchAddress("m", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 2; cut+=0.05){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) < cut){
	    N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  bkgTree->SetBranchAddress("m", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);;
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 2; cut+=0.05){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) < cut){
#if mode == 1
	N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
	float sbratio = 0;
	if(bkg_N[i] != 0)
      sbratio = float(sig_N[i])/sqrt(float(bkg_N[i]));
    SBratio.push_back(sbratio);
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<" "<<sbratio<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
  }
  int dim = SBratio.size();
    for(int i=0; i<dim; i++){
	  SBratio.at(i) = SBratio.at(i)/SBratio.back();
  }
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(color);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,2400,1800);
  c01->SetBottomMargin(0.13);
  c01->SetLeftMargin(0.13);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.15);
  xaxis->SetTitleOffset(1.15);
  xaxis->SetRangeUser(0,2);
  yaxis->SetRangeUser(0,1.5);
  //c01->SetLogy();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  TLegend *legend1 = new TLegend(0.80,0.86,0.899,0.899);
  legend1->AddEntry(gr,"M-"+name,"l");
  legend1->Draw();
  if(plot_CMS) CMS_lumi(c01,iPeriod,iPos);
  c01->Print(name+"_petajeta.png");
  c01->Print(name+"_petajeta.pdf");
  c01->Print(name+"_petajeta.root");
  c01->Print(name+"_petajeta.svg");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

double cut_petajetaptm(TTree* sigTree, TTree* bkgTree, int sigmass, double peta_cut, double jeta_cut, int color, int iPeriod,int iPos, bool plot_CMS){
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
  float xsec_weight = 1;
  float xsec_kfactor = 1;
  float xsec_puweight = 1;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<float> sig_N, bkg_N;
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
  sigTree->SetBranchAddress("m", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 1; cut+=0.01){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm > cut){
	    N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  bkgTree->SetBranchAddress("m", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 1; cut+=0.01){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm > cut){
#if mode == 1
	N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
	float sbratio = 0;
	if(bkg_N[i] != 0)
      sbratio = float(sig_N[i])/sqrt(float(bkg_N[i]));
    SBratio.push_back(sbratio);
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<" "<<sbratio<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
  }
  int dim = SBratio.size();
  float ori = SBratio.front();
  for(int i=0; i<dim; i++){
      cout<<SBratio.at(i)<<",";
	  SBratio.at(i) = SBratio.at(i)/ori;
	  cout<<SBratio.at(i)<<endl;
  }
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(color);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,2400,1800);
  c01->SetBottomMargin(0.13);
  c01->SetLeftMargin(0.13);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.15);
  xaxis->SetTitleOffset(1.15);
  xaxis->SetRangeUser(0,1);
  yaxis->SetRangeUser(0,1.5);
  //c01->SetLogy();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  TLegend *legend1 = new TLegend(0.80,0.86,0.899,0.899);
  legend1->AddEntry(gr,"M-"+name,"l");
  legend1->Draw();
  if(plot_CMS) CMS_lumi(c01,iPeriod,iPos);
  c01->Print(name+"_petajetaptm.png");
  c01->Print(name+"_petajetaptm.pdf");
  c01->Print(name+"_petajetaptm.root");
  c01->Print(name+"_petajetaptm.svg");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

double cut_petajetaptmcos(TTree* sigTree, TTree* bkgTree, int sigmass, double peta_cut, double jeta_cut, double ptm_cut, int color, int iPeriod,int iPos, bool plot_CMS){
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
  float xsec_weight = 1;
  float xsec_kfactor = 1;
  float xsec_puweight = 1;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<float> sig_N, bkg_N;
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
  sigTree->SetBranchAddress("m", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 1; cut+=0.01){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm < ptm_cut) continue;
      if(sys_costhetastar < cut){
	    N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  bkgTree->SetBranchAddress("m", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 1; cut+=0.01){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      if(abs(photon_eta) > peta_cut) continue;
      if(abs(ak8puppijet_eta) > jeta_cut) continue;
      if(sys_ptoverm < ptm_cut) continue;
      if(sys_costhetastar < cut){
#if mode == 1
	N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
	float sbratio = 0;
	if(bkg_N[i] != 0)
      sbratio = float(sig_N[i])/sqrt(float(bkg_N[i]));
    SBratio.push_back(sbratio);
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<" "<<sbratio<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
  }
  int dim = SBratio.size();
    for(int i=0; i<dim; i++){
	  SBratio.at(i) = SBratio.at(i)/SBratio.back();
  }
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(color);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,2400,1800);
  c01->SetBottomMargin(0.13);
  c01->SetLeftMargin(0.13);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.15);
  xaxis->SetTitleOffset(1.15);
  xaxis->SetRangeUser(0,1);
  yaxis->SetRangeUser(0,1.5);
  //c01->SetLogy();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  TLegend *legend1 = new TLegend(0.80,0.86,0.899,0.899);
  legend1->AddEntry(gr,"M-"+name,"l");
  legend1->Draw();
  if(plot_CMS) CMS_lumi(c01,iPeriod,iPos);
  c01->Print(name+"_petajetaptmcos.png");
  c01->Print(name+"_petajetaptmcos.pdf");
  c01->Print(name+"_petajetaptmcos.root");
  c01->Print(name+"_petajetaptmcos.svg");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

double cut_petajetaptmcostau(TTree* sigTree, TTree* bkgTree, int sigmass, double peta_cut, double jeta_cut, double ptm_cut, double cos_cut, int color, int iPeriod,int iPos, bool plot_CMS){
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
  float xsec_weight = 1;
  float xsec_kfactor = 1;
  float xsec_puweight = 1;
  Int_t Nbkg = 0, Nsig = 0;
  // S/sqrt(B) data
  std::vector<float> sig_N, bkg_N;
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
  sigTree->SetBranchAddress("m", &sys_invmass);
  sigTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nsig = sigTree->GetEntries();
  
  for(float cut = 0; cut < 1; cut+=0.01){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<sigTree->GetEntries(); ientry++) {
      // Get Events
      sigTree->GetEntry(ientry);
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      // if(abs(photon_eta) > peta_cut) continue;
      // if(abs(ak8puppijet_eta) > jeta_cut) continue;
      // if(sys_ptoverm < ptm_cut) continue;
      // if(sys_costhetastar > cos_cut) continue;
      if(ak8puppijet_tau21 < cut){
	    N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  bkgTree->SetBranchAddress("m", &sys_invmass);
  bkgTree->SetBranchAddress("sys_seperation", &sys_seperation);
#if mode == 1
  sigTree->SetBranchAddress("xsec_weight", &xsec_weight);
  sigTree->SetBranchAddress("xsec_kfactor", &xsec_kfactor);
  sigTree->SetBranchAddress("xsec_puweight", &xsec_puweight);
#endif
  Nbkg = bkgTree->GetEntries();

  for(float cut = 0; cut < 1; cut+=0.01){
    float N = 0;
    float record = cut;
    for(UInt_t ientry=0; ientry<bkgTree->GetEntries(); ientry++) {
      // Get Events
      bkgTree->GetEntry(ientry);
#if mode == 1
      if(ak8puppijet_masssoftdropcorr < 68 || ak8puppijet_masssoftdropcorr > 94) continue;
#else
      if(ak8puppijet_masssoftdropcorr < 40 || ak8puppijet_masssoftdropcorr > 65) continue;
#endif
      if(abs(sys_invmass - sigmass)  >  0.25 * sigmass) continue;
      // if(abs(photon_eta) > peta_cut) continue;
      // if(abs(ak8puppijet_eta) > jeta_cut) continue;
      // if(sys_ptoverm < ptm_cut) continue;
      // if(sys_costhetastar > cos_cut) continue;
      if(ak8puppijet_tau21 < cut){
#if mode == 1
	N+=xsec_weight*xsec_kfactor*xsec_puweight;
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
  
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  std::vector<float> SBratio;
  std::vector<float> Sigeff;
  std::vector<float> Bkgeff;
  for(int i=0; i<sig_N.size(); i++){
	float sbratio = 0;
	if(bkg_N[i] != 0)
      sbratio = float(sig_N[i])/sqrt(float(bkg_N[i]));
    SBratio.push_back(sbratio);
    cout<<sig_N[i]<<" "<<bkg_N[i]<<" "<<float(sig_N[i]) / float(Nsig)<<" "<<float(bkg_N[i]) / float(Nbkg)<<" "<<sbratio<<endl;
    Sigeff.push_back(float(sig_N[i]) / float(Nsig));
    Bkgeff.push_back(float(bkg_N[i]) / float(Nbkg));
  }
  int dim = SBratio.size();
    for(int i=0; i<dim; i++){
	  SBratio.at(i) = SBratio.at(i)/SBratio.back();
  }
  TGraph* gr = new TGraph(dim,&recordarr[0],&SBratio[0]);
  gr->SetTitle(Graphname);
  gr->SetLineColor(color);
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);

  TCanvas *c01 = new TCanvas("c01",Graphname,2400,1800);
  c01->SetBottomMargin(0.13);
  c01->SetLeftMargin(0.13);
  xaxis = gr->GetXaxis();
  yaxis = gr->GetYaxis();
  xaxis->SetTitle(Graphname);
  yaxis->SetTitle("S/#sqrt{B} a.u.");
  yaxis->SetTitleOffset(1.15);
  xaxis->SetTitleOffset(1.15);
  xaxis->SetRangeUser(0,1);
  yaxis->SetRangeUser(0,2.8);
  //c01->SetLogy();
  c01->SetGrid();
  cout<<"OK"<<endl;
  gr->Draw("AC");
  TString name = std::to_string(sigmass);
  if(plot_CMS) {
    extraText = "Simulation";
    CMS_lumi(c01,iPeriod,iPos);
  }
  TLegend *legend1 = new TLegend(0.80,0.86,0.899,0.899);
  legend1->AddEntry(gr,"M-"+name,"l");
  legend1->Draw();
  c01->Print(name+"_petajetaptmcostau21.png");
  c01->Print(name+"_petajetaptmcostau21.pdf");
  c01->Print(name+"_petajetaptmcostau21.root");
  c01->Print(name+"_petajetaptmcostau21.svg");

  cout << "max S/Sqrt(B): "<<*std::max_element(SBratio.begin(), SBratio.end())<<endl;
  return recordarr.at(std::max_element(SBratio.begin(),SBratio.end()) - SBratio.begin());
}

void significance_count1D(){

  gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.30;
  cmsTextSize = 0.40;
  int iPeriod = 5;
  int iPos = 0;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);
  
  TString width = "W";
  int sigmass = 700;
  int color = 2;
  // int sigmass = 1600;
  // int color = kGreen+3;
  // int sigmass = 2800;
  // int color = 4;
 
  TString sig_sample = "/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/presel_MC/SignalMC"+std::to_string(sigmass)+width+"_postproc_WGamma_full_full_presel_jmcorr_Mar17.root";
  //TString bkg_sample = "/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/presel_data/bkgMC_postproc_WGamma_full_full_presel_jmcorr_Mar17.root";
  TString bkg_sample = "/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/presel_data/SinglePhoton2017_postproc_WGamma_full_full_presel_jmcorr_Mar17.root";

  TFile *infile_1 = TFile::Open(sig_sample);
  TTree* sigTree = (TTree*)infile_1->Get("Events");
  TFile *infile_2 = TFile::Open(bkg_sample);
  TTree* bkgTree = (TTree*)infile_2->Get("Events");
  
  // double jmass_cut = cut_jmass(sigTree,bkgTree,sigmass,color,iPeriod,iPos,plot_CMS);
  // cout<<"jmass_cut: "<<jmass_cut<<endl;

  double peta_cut = cut_peta(sigTree,bkgTree,sigmass,color,iPeriod,iPos,plot_CMS);
  cout<<"peta_cut: "<<peta_cut<<endl;
  
  double jeta_cut = cut_petajeta(sigTree,bkgTree,sigmass,peta_cut,color,iPeriod,iPos,plot_CMS);
  cout<<"jeta_cut: "<<jeta_cut<<endl;
  
  double ptm_cut = cut_petajetaptm(sigTree,bkgTree,sigmass,peta_cut,jeta_cut,color,iPeriod,iPos,plot_CMS);
  cout<<"ptm_cut: "<<ptm_cut<<endl;
  
  double cos_cut = cut_petajetaptmcos(sigTree,bkgTree,sigmass,peta_cut,jeta_cut,ptm_cut,color,iPeriod,iPos,plot_CMS);
  cout<<"cos_cut: "<<cos_cut<<endl;
  
  //double tau21_cut = cut_petajetaptmcostau(sigTree,bkgTree,sigmass,peta_cut,jeta_cut,ptm_cut,cos_cut,color,iPeriod,iPos,plot_CMS);
  // double tau21_cut = cut_petajetaptmcostau(sigTree,bkgTree,sigmass,1.44,2,0.35,0.60,color,iPeriod,iPos,plot_CMS);
  // cout<<"tau_cut: "<<tau21_cut<<endl;

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
