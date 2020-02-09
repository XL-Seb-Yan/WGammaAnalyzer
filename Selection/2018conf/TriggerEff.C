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
#include "TH1F.h"
#include "TRandom.h"

void TriggerEff(TString era = "All"){
  gROOT->SetBatch(1);

  TFile file("Trigger_turnon.root");
  TH1F* hist00lpa = (TH1F*)file.Get("WGamma00lpa");
  TH1F* hist00lpb = (TH1F*)file.Get("WGamma00lpb");
  TH1F* hist00lpc = (TH1F*)file.Get("WGamma00lpc");
  TH1F* hist00hpa = (TH1F*)file.Get("WGamma00hpa");
  TH1F* hist00hpb = (TH1F*)file.Get("WGamma00hpb");
  TH1F* hist00hpc = (TH1F*)file.Get("WGamma00hpc");
  TH1F* hist01lpa = (TH1F*)file.Get("WGamma01lpa");
  TH1F* hist01lpb = (TH1F*)file.Get("WGamma01lpb");
  TH1F* hist01lpc = (TH1F*)file.Get("WGamma01lpc");
  TH1F* hist01hpa = (TH1F*)file.Get("WGamma01hpa");
  TH1F* hist01hpb = (TH1F*)file.Get("WGamma01hpb");
  TH1F* hist01hpc = (TH1F*)file.Get("WGamma01hpc");
  TH1F* hist02lpa = (TH1F*)file.Get("WGamma02lpa");
  TH1F* hist02lpb = (TH1F*)file.Get("WGamma02lpb");
  TH1F* hist02lpc = (TH1F*)file.Get("WGamma02lpc");
  TH1F* hist02hpa = (TH1F*)file.Get("WGamma02hpa");
  TH1F* hist02hpb = (TH1F*)file.Get("WGamma02hpb");
  TH1F* hist02hpc = (TH1F*)file.Get("WGamma02hpc");
  TH1F* hist03lpa = (TH1F*)file.Get("WGamma03lpa");
  TH1F* hist03lpb = (TH1F*)file.Get("WGamma03lpb");
  TH1F* hist03lpc = (TH1F*)file.Get("WGamma03lpc");
  TH1F* hist03hpa = (TH1F*)file.Get("WGamma03hpa");
  TH1F* hist03hpb = (TH1F*)file.Get("WGamma03hpb");
  TH1F* hist03hpc = (TH1F*)file.Get("WGamma03hpc");
  TH1F* hist04lpa = (TH1F*)file.Get("WGamma04lpa");
  TH1F* hist04lpb = (TH1F*)file.Get("WGamma04lpb");
  TH1F* hist04lpc = (TH1F*)file.Get("WGamma04lpc");
  TH1F* hist04hpa = (TH1F*)file.Get("WGamma04hpa");
  TH1F* hist04hpb = (TH1F*)file.Get("WGamma04hpb");
  TH1F* hist04hpc = (TH1F*)file.Get("WGamma04hpc");


  //Tirgger Efficiency
  TEfficiency *Eff11 = new TEfficiency(*hist01lpa, *hist00lpa); //pt
  TEfficiency *Eff12 = new TEfficiency(*hist02lpa, *hist00lpa);
  TEfficiency *Eff13 = new TEfficiency(*hist03lpa, *hist00lpa);
  TEfficiency *Eff14 = new TEfficiency(*hist04lpa, *hist00lpa);
  TEfficiency *Eff21 = new TEfficiency(*hist01hpb, *hist00hpb); //eta
  TEfficiency *Eff22 = new TEfficiency(*hist02hpb, *hist00hpb);
  TEfficiency *Eff23 = new TEfficiency(*hist03hpb, *hist00hpb);
  TEfficiency *Eff24 = new TEfficiency(*hist04hpb, *hist00hpb);
  TEfficiency *Eff31 = new TEfficiency(*hist01hpc, *hist00hpc); //m
  TEfficiency *Eff32 = new TEfficiency(*hist02hpc, *hist00hpc);
  TEfficiency *Eff33 = new TEfficiency(*hist03hpc, *hist00hpc);
  TEfficiency *Eff34 = new TEfficiency(*hist04hpc, *hist00hpc);

  //Non stacked plots
  TLegend *legend2 = new TLegend(0.62,0.16,0.87,0.28);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;

  TCanvas *ceff1 = new TCanvas("ceff1","",1200,900);
  ceff1->cd();
  Eff11->SetTitle("Trigger Efficiency vs photon pT ("+era+")");
  Eff11->SetLineWidth(2);
  Eff11->SetLineColor(46);
  Eff11->Draw();
  gPad->Update();
  Eff11->SetTitle("Trigger Efficiency on photon pT ("+era+"); pt_{#gamma} (GeV); Efficiency"); 
  Eff11->GetPaintedGraph()->GetYaxis()->SetRangeUser(-0.1,1.1);
  Eff11->Draw("SAME");
  Eff12->SetLineWidth(2);
  Eff12->SetLineColor(42);
  Eff12->Draw("SAME");
  Eff13->SetLineWidth(2);
  Eff13->SetLineColor(28);
  Eff13->Draw("SAME");
  Eff14->SetLineWidth(2);
  Eff14->SetLineColor(9);
  Eff14->Draw("SAME");
  legend2->Clear();
  legend2->AddEntry(Eff11,"HLT_Photon110EB_TightID_TightIso","f");
  legend2->AddEntry(Eff12,"HLT_Photon120EB_TightID_TightIso","f");
  legend2->AddEntry(Eff13,"HLT_Photon200","f");
  legend2->AddEntry(Eff14,"Combined","f");
  legend2->Draw();
  ceff1->Print("eff_pt.png");
  ceff1->Print("eff_pt.pdf");

  TCanvas *ceff2 = new TCanvas("ceff2","",1200,900);
  ceff2->cd();
  Eff21->SetTitle("Trigger Efficiency vs photon #eta ("+era+")");
  Eff21->SetLineWidth(2);
  Eff21->SetLineColor(46);
  Eff21->Draw();
  gPad->Update();
  Eff21->SetTitle("Trigger Efficiency on photon #eta ("+era+"); #eta_{#gamma}; Efficiency"); 
  Eff21->GetPaintedGraph()->GetYaxis()->SetRangeUser(-0.1,1.1);
  Eff21->Draw("SAME");
  Eff22->SetLineWidth(2);
  Eff22->SetLineColor(42);
  Eff22->Draw("SAME");
  Eff23->SetLineWidth(2);
  Eff23->SetLineColor(28);
  Eff23->Draw("SAME");
  Eff24->SetLineWidth(2);
  Eff24->SetLineColor(9);
  Eff24->Draw("SAME");
  legend2->Clear();
  legend2->AddEntry(Eff11,"HLT_Photon110EB_TightID_TightIso","f");
  legend2->AddEntry(Eff12,"HLT_Photon120EB_TightID_TightIso","f");
  legend2->AddEntry(Eff13,"HLT_Photon200","f");
  legend2->AddEntry(Eff14,"Combined","f");
  legend2->Draw();
  ceff2->Print("eff_eta.png");
  ceff2->Print("eff_eta.pdf");

  TCanvas *ceff3 = new TCanvas("ceff3","",1200,900);
  ceff3->cd();
  Eff31->SetTitle("Trigger Efficiency vs M_{j#gamma} ("+era+")");
  Eff31->SetLineWidth(2);
  Eff31->SetLineColor(46);
  Eff31->Draw();
  gPad->Update();
  Eff31->SetTitle("Trigger Efficiency on M_{j#gamma} ("+era+"); M_{j#gamma} (GeV); Efficiency"); 
  Eff31->GetPaintedGraph()->GetYaxis()->SetRangeUser(-0.1,1.1);
  Eff31->Draw("SAME");
  Eff32->SetLineWidth(2);
  Eff32->SetLineColor(42);
  Eff32->Draw("SAME");
  Eff33->SetLineWidth(2);
  Eff33->SetLineColor(28);
  Eff33->Draw("SAME");
  Eff34->SetLineWidth(2);
  Eff34->SetLineColor(9);
  Eff34->Draw("SAME");
  legend2->Clear();
  legend2->AddEntry(Eff11,"HLT_Photon110EB_TightID_TightIso","f");
  legend2->AddEntry(Eff12,"HLT_Photon120EB_TightID_TightIso","f");
  legend2->AddEntry(Eff13,"HLT_Photon200","f");
  legend2->AddEntry(Eff14,"Combined","f");
  legend2->Draw();
  ceff3->Print("eff_m.png");
  ceff3->Print("eff_m.pdf");

}
