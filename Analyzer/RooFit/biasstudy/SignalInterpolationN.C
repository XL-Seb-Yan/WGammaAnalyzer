#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooIntegralMorph.h"
#include "RooNLLVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TH1.h"

using namespace RooFit;

void SignalInterpolationN(){

  gROOT->SetBatch(1);
  // Observable
  RooRealVar m("m","m",600,4000);
    
  const double step = 50;
  const int nMCpoints = 14;  
  RooAbsPdf* gMass[nMCpoints];   
  const double masses[nMCpoints] = {700,800,900,1000,1200,1400,1800,2000,2200,2400,2600,2800,3000,3500};
  //const double masses[nMCpoints] = {3000,3500};

  TFile *f[nMCpoints];
  RooWorkspace* xf[nMCpoints];

  for (int i = 0; i!=nMCpoints; ++i ){
    TString massname = std::to_string(int(masses[i]));
    TString name = "/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/RooFitWorkspace/anchor/"+massname+"N-shapes-Unbinned-CB.root";
    if (!gSystem->AccessPathName(name)){
      f[i] = new TFile(name);
      xf[i] = (RooWorkspace*)f[i]->Get("w");
      xf[i]->var("m")->setRange(600,4000);
      gMass[i] = xf[i]->pdf("CBShape");
    } else {
      std::cout<<"File is not found: "<<name<<std::endl;
      return;
    }
  }
  
  RooWorkspace w("w");
  w.import(*gMass[0],RooFit::RenameVariable("CBShape","CBShape_low"),RooFit::RenameAllVariablesExcept("low","m"));
  gMass[0] = w.pdf("CBShape_low");
  for (int i = 1; i!=nMCpoints; ++i ) {
    TString name = Form("point_%d",i);
    w.import(*gMass[i],RooFit::RenameConflictNodes(name),RooFit::RenameAllVariablesExcept(name,"m"),RooFit::RenameVariable("CBShape","CBShape_low"+name));
        
    gMass[i] = w.pdf("CBShape_low"+name);    
  }

  TF1* funcEff = new TF1("effFunc","landau(0)",600,4000);
  funcEff->SetParameters(0.7992,-15.5926,2450.04);
  //funcEff->SetParameters(0.55203,803.672,337.859);
    
  // C r e a t e   i n t e r p o l a t i n g   p d f
  // -----------------------------------------------
    
  // Create interpolation variable
  RooRealVar alpha("alpha","alpha",0,1.0) ;
    
  // Specify sampling density on observable and interpolation variable
  m.setBins(3400,"cache") ;
  alpha.setBins(2000,"cache") ;
    
  RooPlot* frame1[nMCpoints];
  RooPlot* frame2[nMCpoints];
  RooPlot* frame3[nMCpoints];
  TH1* hh[nMCpoints];
    
  TCanvas* c[nMCpoints];
    
  for (int iPoint = 0; iPoint!=nMCpoints-1; ++iPoint) {
    cout<<"================================================================"<<endl;
    cout<<"===================Processing "<<iPoint<<" ====================="<<endl; 
    //for (int iPoint = 1; iPoint!=2; ++iPoint) {
            
    // Construct interpolating pdf in (m,a) represent g1(m) at a=a_min and g2(m) at a=a_max
    cout<<iPoint<<" "<<gMass[iPoint]<<" "<<gMass[iPoint+1]<<endl;
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
      TFile* fileNew = new TFile(Form("GenSignal/roodataset_signal-%d-narrow.root",int(masses[iPoint+1]-i*step)),"RECREATE");
      dataGen->Write();
      distribs0->Write();
      fileNew->Close();
            
      delete dataGen;
      delete distribs0;
    }
    TCanvas *c = new TCanvas("","",1200,900);
    c->cd();
    c->Clear();
    frame1[iPoint]->SetTitle("Signal Interpolation (Narrow)");
    TAxis *xaxis = frame1[iPoint]->GetXaxis();
    TAxis *yaxis = frame1[iPoint]->GetYaxis();
    xaxis->SetTitle("m_{W#gamma}");
    yaxis->SetTitle("Events (a.u.)");
    xaxis->SetRangeUser(0,4000);
    yaxis->SetRangeUser(0,0.5);
    frame1[iPoint]->Draw();
    TString pngname = std::to_string(int(masses[iPoint]));
    c->Print(pngname+"N.png");

    /*
    // S h o w   2 D   d i s t r i b u t i o n   o f   p d f ( m , a l p h a )
    // -----------------------------------------------------------------------
        
    // Create 2D histogram
    hh[iPoint] = lmorph.createHistogram("hh",m,Binning(40),YVar(alpha,Binning(40))) ;
    hh[iPoint]->SetLineColor(kBlue) ;
        
        
    // F i t   p d f   t o   d a t a s e t   w i t h   a l p h a = 0 . 8
    // -----------------------------------------------------------------
        
    // Generate a toy dataset at alpha = 0.8
    alpha=0.8 ;
    RooDataSet* data = lmorph.generate(m,1000) ;
        
    // Fit pdf to toy data
    lmorph.setCacheAlpha(kTRUE) ;
    lmorph.fitTo(*data,Verbose(kTRUE)) ;
        
    // Plot fitted pdf and data overlaid
    frame2[iPoint] = m.frame(Bins(1500)) ;
    data->plotOn(frame2[iPoint]) ;
    lmorph.plotOn(frame2[iPoint]) ;
        
        
    // S c a n   - l o g ( L )   v s   a l p h a
    // -----------------------------------------
        
    // Show scan -log(L) of dataset w.r.t alpha
    frame3[iPoint] = alpha.frame(Bins(100),Range(0.1,0.9)) ;
        
    // Make 2D pdf of histogram
    RooNLLVar nll("nll","nll",lmorph,*data) ;
    nll.plotOn(frame3[iPoint],ShiftToZero()) ;
        
    lmorph.setCacheAlpha(kFALSE) ;
        
        
        
    c[iPoint] = new TCanvas(Form("linearmorph_%d",iPoint),Form("linearmorph_%d",iPoint),700,700) ;
    c[iPoint]->Divide(2,2) ;
    c[iPoint]->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1[iPoint]->GetYaxis()->SetTitleOffset(1.6) ; frame1[iPoint]->Draw() ;
    c[iPoint]->cd(2) ; gPad->SetLeftMargin(0.20) ; hh[iPoint]->GetZaxis()->SetTitleOffset(2.5) ; hh[iPoint]->Draw("surf") ;
    c[iPoint]->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3[iPoint]->GetYaxis()->SetTitleOffset(1.4) ; frame3[iPoint]->Draw() ;
    c[iPoint]->cd(4) ; gPad->SetLeftMargin(0.15) ; frame2[iPoint]->GetYaxis()->SetTitleOffset(1.4) ; frame2[iPoint]->Draw() ;
    c[iPoint]->Modified(); c[iPoint]->Update();
    */
  }
    
  return; 
}
