#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooBukinPdf.h"
#include "RooIntegralMorph.h"
#include "RooNLLVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TH1.h"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"

using namespace RooFit;

void SignalInterpolationN(){

  gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 5;
  int iPos = 33;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.025);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);
  
  // Observable
  RooRealVar m("m","m",600,4000);
    
  const double step = 50;
  const int nMCpoints = 14;  
  RooAbsPdf* gMass[nMCpoints];   
  const double masses[nMCpoints] = {700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3500};
  //const double masses[nMCpoints] = {3000,3500};

  TFile *f[nMCpoints];
  RooWorkspace* xf[nMCpoints];

  for (int i = 0; i!=nMCpoints; ++i ){
    TString massname = std::to_string(int(masses[i]));
    TString name = "/afs/cern.ch/work/x/xuyan/work5/PROD17/Analyzer/CMSSW_9_4_13/src/WGammaAnalyzer/Analyzer/RooFit/biasstudy/"+massname+"N-shapes-Unbinned-CBGaus.root";
    if (!gSystem->AccessPathName(name)){
      f[i] = new TFile(name);
      xf[i] = (RooWorkspace*)f[i]->Get("w");
      xf[i]->var("m")->setRange(600,4000);
      gMass[i] = xf[i]->pdf("CBGaus");
    } else {
      std::cout<<"File is not found: "<<name<<std::endl;
      return;
    }
  }
  
  RooWorkspace w("w");
  w.import(*gMass[0],RooFit::RenameVariable("CBGaus","CBGaus_low"),RooFit::RenameAllVariablesExcept("low","m"));
  gMass[0] = w.pdf("CBGaus_low");
  for (int i = 1; i!=nMCpoints; ++i ) {
    TString name = Form("point_%d",i);
    w.import(*gMass[i],RooFit::RenameConflictNodes(name),RooFit::RenameAllVariablesExcept(name,"m"),RooFit::RenameVariable("CBGaus","CBGaus"+name));
        
    gMass[i] = w.pdf("CBGaus"+name);    
  }
    
  // C r e a t e   i n t e r p o l a t i n g   p d f
  // -----------------------------------------------
    
  // Create interpolation variable
  RooRealVar alpha("alpha","alpha",0,1.0) ;
    
  // Specify sampling density on observable and interpolation variable
  m.setBins(3400,"cache") ;
  alpha.setBins(1000,"cache") ;
    
  RooPlot* frame1[nMCpoints];
  RooPlot* frame2[nMCpoints];
  RooPlot* frame3[nMCpoints];
  TH1* hh[nMCpoints];
  
  for (int iPoint = 0; iPoint!=nMCpoints-1; ++iPoint) {
    cout<<"================================================================"<<endl;
    cout<<"===================Processing "<<iPoint<<" ====================="<<endl; 
    //for (int iPoint = 1; iPoint!=2; ++iPoint) {
            
    // Construct interpolating pdf in (m,a) represent g1(m) at a=a_min and g2(m) at a=a_max
    cout<<"Running on: "<<iPoint<<" "<<gMass[iPoint]<<" "<<gMass[iPoint+1]<<endl;
    RooIntegralMorph lmorph("lmorph","lmorph",*gMass[iPoint],*gMass[iPoint+1],m,alpha) ;
        
    // P l o t   i n t e r p o l a t i n g   p d f   a t   v a r i o u s   a l p h a
    // -----------------------------------------------------------------------------

    // Show end points as blue curves
    frame1[iPoint] = m.frame() ;
    gMass[iPoint]->plotOn(frame1[iPoint]) ;
    gMass[iPoint+1]->plotOn(frame1[iPoint]) ;
        
    int nPoints = int((masses[iPoint+1]-masses[iPoint])/step);
    for (int i=1; i!=nPoints; ++i) {
      cout<<"================================================================"<<endl;
      cout<<"===================Processing "<<i<<" ====================="<<endl; 
      //alpha=double(i)/double(nPoints);
      alpha.setVal(double(i)/double(nPoints)) ;
      lmorph.plotOn(frame1[iPoint],LineColor(kRed));
            
      RooDataSet* dataGen = lmorph.generate(m,30000);
      TH1D* distribs0 = (TH1D*)lmorph.createHistogram("distribs_5_10_0",m,Binning(170,600,4000));
      dataGen->fillHistogram(distribs0,m);
      /*
      double effcy = funcEff->Eval(masses[iPoint+1]-i*step);
      //double weight = effcy/30000.;
      double weight = effcy;
      //TH1D* distribs0 = new TH1D("distribs_5_10_0","",4000,0,4000);
      
            
      double effcy2 = distribs0->Integral(600,3500);
      cout<<"Integral is: "<<effcy2<<endl;
            
      
      //distribs0->Scale(weight/effcy2);
      */
      TFile* fileNew = new TFile(Form("PdfGenerateSignal/roodataset_signal-%d-narrow.root",int(masses[iPoint+1]-i*step)),"RECREATE");
      dataGen->Write();
      distribs0->Write();
      fileNew->Close();
            
      delete dataGen;
      delete distribs0;
    }
  }
  TCanvas *c = new TCanvas("","",2400,1800);
  c->cd();
  c->Clear();
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.105);
  frame1[0]->SetTitle("Signal Interpolation (Narrow)");
  TAxis *xaxis = frame1[0]->GetXaxis();
  TAxis *yaxis = frame1[0]->GetYaxis();
  xaxis->SetTitle("m_{W#gamma}");
  yaxis->SetTitle("Events (a.u.)");
  yaxis->SetTitleOffset(1.1);
  xaxis->SetRangeUser(0,4000);
  yaxis->SetRangeUser(0,0.5);
  frame1[0]->Draw();
  for (int iPoint = 1; iPoint!=nMCpoints-1; ++iPoint) {
    frame1[iPoint]->Draw("SAME");
  }
  CMS_lumi(c,iPeriod,iPos);
  c->Print("Signal_interpolation_N.png");
  c->Print("Signal_interpolation_N.pdf");
  c->Print("Signal_interpolation_N.root");
  c->Print("Signal_interpolation_N.svg");
    
  return; 
}
