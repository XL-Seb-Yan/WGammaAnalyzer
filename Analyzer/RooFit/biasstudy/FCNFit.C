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
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
 
#include <cassert>
#endif
 
using namespace ROOT::Math;

bool first = true; 
// function Object to be minimized
struct Chi2NegDriv {
   // the TGraph is a data member of the object
   TH1 *fHist;
   TString fun;
  
   Chi2NegDriv (TH1 *h, TString fun_){
	fHist = h;
	fun = fun_;
   }
 
   // Chi2NegDriv(TH1 *h) : fHist(h) {}
   
   std::vector<double> dijet_1(double m, const double *par){
	   TF1 f1("dijet_1","pow(x/13000,[0])",600,6000);
	   f1.SetParameter(0,par[0]);
	   double Norm = f1.Integral(600,6000);
	   TF1 f2("dijet_1_pdf","[0]*pow(x/13000,[1])",600,6000);
	   f2.SetParameters(1.0/Norm, par[0]);
	   std::vector<double> re;
	   re.push_back(f2.Eval(m));
	   re.push_back(f2.Derivative(m));
	   return re;
   }
   
   std::vector<double> VVdijet_2(double m, const double *par){
	   TF1 f1("VVdijet_2","pow(1-x/13000,[0])/pow(x/13000,[1]+[2]*log(x/13000))",600,6000);
	   f1.SetParameters(par[0],par[1],par[2]);
	   double Norm = f1.Integral(600,6000);
	   TF1 f2("VVdijet_2_pdf","[0]*pow(1-x/13000,[1])/pow(x/13000,[2]+[3]*log(x/13000))",600,6000);
	   f2.SetParameters(1.0/Norm,par[0],par[1],par[2]);
	   std::vector<double> re;
	   re.push_back(f2.Eval(m));
	   re.push_back(f2.Derivative(m));
	   return re;
   }
   
   std::vector<double> VVdijet_3(double m, const double *par){
	   TF1 f1("VVdijet_3","pow(1-x/13000,[0])/pow(x/13000,[1]+[2]*log(x/13000)+[3]*pow(log(x/13000),2))",600,6000);
	   f1.SetParameters(par[0],par[1],par[2],par[3]);
	   double Norm = f1.Integral(600,6000);
	   TF1 f2("VVdijet_3_pdf","[0]*pow(1-x/13000,[1])/pow(x/13000,[2]+[3]*log(x/13000)+[4]*pow(log(x/13000),2))",600,6000);
	   f2.SetParameters(1.0/Norm,par[0],par[1],par[2],par[3]);
	   std::vector<double> re;
	   re.push_back(f2.Eval(m));
	   re.push_back(f2.Derivative(m));
	   return re;
   }
   
   std::vector<double> ATLAS_2(double m, const double *par){
	   TF1 f1("ATLAS_2","pow(1-pow(x/13000,0.3333333),[0])/pow(x/13000,[1]+[2]*log(x/13000))",600,6000);
	   f1.SetParameters(par[0],par[1],par[2]);
	   double Norm = f1.Integral(600,6000);
	   //cout<<Norm<<endl;
	   TF1 f2("ATLAS_2_pdf","[0]*pow(1-pow(x/13000,0.3333333),[1])/pow(x/13000,[2]+[3]*log(x/13000))",600,6000);
	   f2.SetParameters(1.0/Norm,par[0],par[1],par[2]);
	   //cout<<f2.Integral(600,6000)<<endl;
	   std::vector<double> re;
	   re.push_back(f2.Eval(m));
	   re.push_back(f2.Derivative(m));
	   return re;
   }
 
   // implementation of the function to be minimized
   double operator() (const double *par) {
      assert(fHist != 0);
	  double sum = 0;
      for (int i  = 1; i < fHist->GetNbinsX()+1; ++i) {
		  double NData = fHist->GetBinContent(i);
		  double m = fHist->GetBinCenter(i);
		  double NPredict = -99;
		  double fun_driv = 0;
		  if(fun == "dijet_1"){
			NPredict = (dijet_1(m, par))[0] * 21435 * fHist->GetBinWidth(i);
			fun_driv = (dijet_1(m, par))[1]; 
		  }
		  else if (fun == "VVdijet_2"){
			NPredict = (VVdijet_2(m, par))[0] * 21435 * fHist->GetBinWidth(i);
			fun_driv = (VVdijet_2(m, par))[1]; 
		  }
		  else if (fun == "VVdijet_3"){
			NPredict = (VVdijet_3(m, par))[0] * 21435 * fHist->GetBinWidth(i);
			fun_driv = (VVdijet_3(m, par))[1]; 
		  }
		  else if (fun == "ATLAS_2"){
			NPredict = (ATLAS_2(m, par))[0] * 21435 * fHist->GetBinWidth(i);
			fun_driv = (ATLAS_2(m, par))[1]; 
		  }
		  if(NPredict == -99) {
			  cout<<"Invalid NPredict"<<endl;
			  break;
		  }
		  double error = fHist->GetBinError(i);
		  if (error == 0)
			  error = sqrt(NPredict);
		  
		  //cout<<m<<" "<<NData<<" "<<NPredict<<" "<<error<<" "<<fun_driv<<endl;
		  double multiplier = (fun_driv>0) ? fun_driv*10E14 : 0;
		  double chi2_final = pow((NData-NPredict)/error,2) + multiplier;
		  //double chi2_final = pow((NData-NPredict)/error,2);
		  sum += chi2_final;
      }
	  
      if(first)
         std::cout << "Total Initial distance square = " << sum << std::endl;
      first = false;
	  std::cout << "modified chi2 = " << sum << std::endl;
	  cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<endl;
      return sum;
   }
 
};
 
int FCNFit(){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();
   
	TH1F* DATAh = new TH1F("DATAh","DATAh",160,600,7000);
	float s_mass;
	TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/fullcut/Run2Data_postproc_WGammaRun2_SR_sigrange_fullcut_jmcorr_May22.root");
	TTree* tree = (TTree*)file.Get("Events");
	tree->SetBranchAddress("m", &s_mass);
	for (int ievt = 0; ievt<tree->GetEntries();ievt++) {
		tree->GetEntry(ievt);
		DATAh->Fill(s_mass, 1);
		//DATAh->Fill(s_mass, 2.4826);
	}
	DATAh->SetBinErrorOption(TH1::kPoisson);

   ROOT::Fit::Fitter  fitter;

   // make the functor objet
   Chi2NegDriv schi(DATAh, "VVdijet_3");
   ROOT::Math::Functor fcn(schi,4);
   // set the function and the initial parameter values
   double pStart[4] = {-23.3227,15.13,1.10,-0.14};
   fitter.SetFCN(fcn,pStart);
   fitter.Config().ParSettings(0).SetLimits(-35,-10);
   fitter.Config().ParSettings(1).SetLimits(2,27);
   fitter.Config().ParSettings(2).SetLimits(0,3);
   fitter.Config().ParSettings(3).SetLimits(-1,0.01);
   //set step sizes different than default ones (0.3 times parameter values)
   // for(int i = 0; i < 1; ++i) 
	   // fitter.Config().ParSettings(i).SetStepSize(0.01);
   bool ok = fitter.FitFCN();
   if (!ok) {
      Error("FCNFit","FCNFit failed");
      return 1;
   }
 
   const ROOT::Fit::FitResult & result = fitter.Result();
 
   std::cout << "Total chi2 " << result.MinFcnValue() << std::endl;
   result.Print(std::cout);
 /*
 
   gr->Draw("p0");
 
   // get fit parameters
   const double * parFit = result.GetParams();
 
   // draw the fitted line
   int n = 1000;
   double t0 = 0;
   double dt = 10;
   TPolyLine3D *l = new TPolyLine3D(n);
   for (int i = 0; i <n;++i) {
      double t = t0+ dt*i/n;
      double x,y,z;
      line(t,parFit,x,y,z);
      l->SetPoint(i,x,y,z);
   }
   l->SetLineColor(kRed);
   l->Draw("same");
 
   // draw original line
   TPolyLine3D *l0 = new TPolyLine3D(n);
   for (int i = 0; i <n;++i) {
      double t = t0+ dt*i/n;
      double x,y,z;
      line(t,p0,x,y,z);
      l0->SetPoint(i,x,y,z);
   }
   l0->SetLineColor(kBlue);
   l0->Draw("same");
   */
   return 0;
}