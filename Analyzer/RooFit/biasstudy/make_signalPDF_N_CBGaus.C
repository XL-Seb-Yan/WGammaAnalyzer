//make signal shapes
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
#include <RooMsgService.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooCBShape.h>
#include <RooGaussian.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooChi2Var.h>
#include <TLatex.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooPlotable.h>
#include <RooWorkspace.h>
#include <RooAddPdf.h>
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
#endif

std::string to_str_trim(const float a_value, const int n = 2)
{
    return std::to_string(a_value).substr(0,std::to_string(a_value).find(".") + n + 1);
}

void make_signalPDF_N_CBGaus()
{
  //gErrorIgnoreLevel = kInfo;
  using namespace std;
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37);
  
  gROOT->SetBatch(1);
  lumi_13TeV = "";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 5;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.045,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetBarWidth(1.03);
  gStyle->SetHistLineWidth(2);

  // --- Create obervable --- 
  RooRealVar *m = new RooRealVar("m","m",600,4000,""); //the name "m" will be used by RooDataSet to import data
  
  bool redoAnchor = false;
  
  //--- initial parameters
  float CB_mean_low, CB_mean_nominal, CB_mean_high;
  float CB_sigma_low, CB_sigma_nominal, CB_sigma_high;
  float CB_alpha_low, CB_alpha_nominal, CB_alpha_high;
  float CB_n_low, CB_n_nominal, CB_n_high;
  float Gaus_mean_low, Gaus_mean_nominal, Gaus_mean_high;
  float Gaus_sigma_low, Gaus_sigma_nominal, Gaus_sigma_high;
  float frac_low, frac_nominal, frac_high;
  
  // for(int i=700; i<3501; i+=100){
	  // if(!redoAnchor) break;
	  // std::string mass_str = std::to_string(i);
	  // TFile *f = NULL;
	  // TString filename = "/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/RooFitWorkspace/N/"+mass_str+"N-shapes-Unbinned-CBGaus.root";
	  // f = TFile::Open(filename);
	  // if(f == NULL) {
		  // cout<<"root file not found!"<<endl;
		  // continue;
		  // }
	  // RooWorkspace *w = (RooWorkspace*)f->Get("w");
	  // CB_mean_nominal = w->var("CB_mean")->getValV();
	  // CB_sigma_nominal = w->var("CB_sigma")->getValV();
	  // CB_alpha_nominal = w->var("CB_alpha")->getValV();
	  // CB_n_nominal = w->var("CB_n")->getValV();
	  // Gaus_mean_nominal = w->var("Gaus_mean")->getValV();
	  // Gaus_sigma_nominal = w->var("Gaus_sigma")->getValV();
	  // frac_nominal = w->var("frac")->getValV();
	  
	  // TString fun_name = "CBGaus";
	  // RooRealVar* CB_mean = new RooRealVar("CB_mean","CB_mean",CB_mean_nominal);
	  // RooRealVar* CB_sigma = new RooRealVar("CB_sigma","CB_sigma",CB_sigma_nominal);
	  // RooRealVar* CB_alpha = new RooRealVar("CB_alpha","CB_alpha",CB_alpha_nominal);
	  // RooRealVar* CB_n = new RooRealVar("CB_n","CB_n",CB_n_nominal);
	  // RooCBShape* CB_model = new RooCBShape("CBShape","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
	  
	  // RooRealVar* Gaus_mean = new RooRealVar("Gaus_mean","Gaus_mean",Gaus_mean_nominal);
	  // RooRealVar* Gaus_sigma = new RooRealVar("Gaus_sigma","Gaus_sigma",Gaus_sigma_nominal);
	  // RooGaussian* Gaus_model = new RooGaussian("Gaussian","Gaussian Function",*m,*Gaus_mean,*Gaus_sigma);
	  // RooRealVar* frac = new RooRealVar("frac","frac",frac_nominal);
	  // RooAddPdf* com_model = new RooAddPdf("CBGaus","CBGaus",RooArgList(*CB_model,*Gaus_model),RooArgList(*frac));
	  
	  // RooDataSet* dataGen = com_model->generate(*m,30000,RooFit::Name("ImorphData"));
	  // TH1D* distribs0 = (TH1D*)com_model->createHistogram("GeneratedData",*m,Binning(170,600,4000));
	  // dataGen->fillHistogram(distribs0,*m);
	  // TFile* fileNew = new TFile(Form("roodataset_signal-%d-narrow.root",int(i)),"RECREATE");
	  // distribs0->Write();
	  // fileNew->Close();
	  // delete fileNew;
  // }
  
  for(int i=3150; i<3251; i+=50){

	  if(redoAnchor) break;
	  std::string mass_str = std::to_string(i);
	  TFile *f = NULL;
	  TString filename = "/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/2017/RooFitWorkspace/N/PdfGenerateSignal/roodataset_signal-"+mass_str+"-narrow.root";
	  f = TFile::Open(filename);
	  if(f == NULL) {
		  cout<<"root file not found: "<<mass_str<<endl;
		  continue;
	  }
	  TString fun_name = "CBGaus";
	  // nominal settings
	  // RooRealVar* CB_mean = new RooRealVar("CB_mean","CB_mean",i,i-100,i+100,"");
	  // RooRealVar* CB_sigma = new RooRealVar("CB_sigma","CB_sigma",50,20,110,"");
	  // RooRealVar* CB_alpha = new RooRealVar("CB_alpha","CB_alpha",1.5,0.1,3,"");
	  // RooRealVar* CB_n = new RooRealVar("CB_n","CB_n",1,0.1,5,"");
	  // RooCBShape* CB_model = new RooCBShape("CBShape","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
	  
	  // RooRealVar* Gaus_mean = new RooRealVar("Gaus_mean","Gaus_mean",i,i-100,i+100,"");
	  // RooRealVar* Gaus_sigma = new RooRealVar("Gaus_sigma","Gaus_sigma",70,20,150,"");
	  // RooGaussian* Gaus_model = new RooGaussian("Gaussian","Gaussian Function",*m,*Gaus_mean,*Gaus_sigma);
	  // RooRealVar* frac = new RooRealVar("frac","frac",0.7,0.5,1);
	  // RooAddPdf* com_model = new RooAddPdf("CBGaus","CBGaus",RooArgList(*CB_model,*Gaus_model),RooArgList(*frac));
	  
	  RooRealVar* CB_mean = new RooRealVar("CB_mean","CB_mean",i,i-100,i+100,"");
	  RooRealVar* CB_sigma = new RooRealVar("CB_sigma","CB_sigma",70,40,90,"");
	  RooRealVar* CB_alpha = new RooRealVar("CB_alpha","CB_alpha",2,1,3,"");
	  RooRealVar* CB_n = new RooRealVar("CB_n","CB_n",3,0.5,5,"");
	  RooCBShape* CB_model = new RooCBShape("CBShape","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
	  
	  RooRealVar* Gaus_mean = new RooRealVar("Gaus_mean","Gaus_mean",i,i-100,i+100,"");
	  RooRealVar* Gaus_sigma = new RooRealVar("Gaus_sigma","Gaus_sigma",100,50,120,"");
	  RooGaussian* Gaus_model = new RooGaussian("Gaussian","Gaussian Function",*m,*Gaus_mean,*Gaus_sigma);
	  RooRealVar* frac = new RooRealVar("frac","frac",0.75,0.5,1);
	  RooAddPdf* com_model = new RooAddPdf("CBGaus","CBGaus",RooArgList(*CB_model,*Gaus_model),RooArgList(*frac));
	  
	  TH1F *MChist = (TH1F*)f->Get("GeneratedData__m");
	  
	  cout<<MChist<<endl;

  
	RooDataHist datah("Signal MC (W band)","Signal MC (W band)",RooArgSet(*m),MChist);

	RooFitResult *r = com_model->fitTo(datah,Range(i*0.75,i*1.25),RooFit::Minimizer("Minuit2"),SumW2Error(false),Save()); //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)


	// --- Plot ---
	gStyle->SetOptStat(111111);
	RooPlot *frame = m->frame();
	datah.plotOn(frame,RooFit::Name("data"));
	frame->getAttMarker()->SetMarkerSize(2);
	com_model->plotOn(frame,LineColor(2),RooFit::Name("fit"),Range(i*0.75,i*1.25));
	//com_model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
	com_model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"),Range(i*0.75,i*1.25));
	com_model->plotOn(frame,Components("CBShape"),LineStyle(kDashed),LineColor(1),RooFit::Name("CB"),Range(i*0.75,i*1.25));
	com_model->plotOn(frame,Components("Gaussian"),LineStyle(kDashed),RooFit::Name("Gaus"),Range(i*0.75,i*1.25));
	com_model->plotOn(frame,LineColor(2),RooFit::Name("fit"),Range(i*0.75,i*1.25));
	datah.plotOn(frame);
	frame->getAttMarker()->SetMarkerSize(2);

	TString chi2txt = "#chi^{2}/ndf: "+to_str_trim(frame->chiSquare("fit","data", 7));
	TString NLLtxt = "NLL: "+to_str_trim(com_model->createNLL(datah)->getVal());
	TString CBmean = "#mu_{CB}: "+to_str_trim(CB_mean->getVal())+" #pm "+to_str_trim(CB_mean->getError())+" (GeV)";
	TString CBsigma = "#sigma_{CB}: "+to_str_trim(CB_sigma->getVal())+" #pm "+to_str_trim(CB_sigma->getError())+" (GeV)";
	TString CBalpha = "#alpha_{CB}: "+to_str_trim(CB_alpha->getVal())+" #pm "+to_str_trim(CB_alpha->getError());
	TString CBn = "n_{CB}: "+to_str_trim(CB_n->getVal())+" #pm "+to_str_trim(CB_n->getError());
	TString Gaussmean = "#mu_{Gaus}: "+to_str_trim(Gaus_mean->getVal())+" #pm "+to_str_trim(Gaus_mean->getError())+" (GeV)";
	TString Gaussigma = "#sigma_{Gaus}: "+to_str_trim(Gaus_sigma->getVal())+" #pm "+to_str_trim(Gaus_sigma->getError())+" (GeV)";
	TString Fraction = "CB frac: "+to_str_trim(frac->getVal())+" #pm "+to_str_trim(frac->getError());
	TLatex *chi2lax = new TLatex(0.15,0.7-0*0.04,chi2txt);
	//TLatex *NLLlax = new TLatex(0.15,0.8-0.06,NLLtxt);
	TLatex *CBmeanlax = new TLatex(0.15,0.7-1*0.04,CBmean);
	TLatex *CBsigmalax = new TLatex(0.15,0.7-2*0.04,CBsigma);
	TLatex *CBalphalax = new TLatex(0.15,0.7-3*0.04,CBalpha);
	TLatex *CBnlax = new TLatex(0.15,0.7-4*0.04,CBn);
	TLatex *Gaussigmalax = new TLatex(0.15,0.7-5*0.04,Gaussigma);
	TLatex *Gaussmeanlax = new TLatex(0.15,0.7-6*0.04,Gaussmean);
	TLatex *Fraclax = new TLatex(0.15,0.7-7*0.04,Fraction);

	CBmeanlax->SetNDC();
	CBsigmalax->SetNDC();
	CBalphalax->SetNDC();
	CBnlax->SetNDC();
	Gaussigmalax->SetNDC();
	Gaussmeanlax->SetNDC();
	Fraclax->SetNDC();

	CBmeanlax->SetTextSize(0.026);
	CBsigmalax->SetTextSize(0.026);
	CBalphalax->SetTextSize(0.026);
	CBnlax->SetTextSize(0.026);
	Gaussigmalax->SetTextSize(0.026);
	Gaussmeanlax->SetTextSize(0.026);
	Fraclax->SetTextSize(0.026);

	frame->addObject(CBmeanlax);
	frame->addObject(CBsigmalax);
	frame->addObject(CBalphalax);
	frame->addObject(CBnlax);
	frame->addObject(Gaussigmalax);
	frame->addObject(Gaussmeanlax);
	frame->addObject(Fraclax);

	TCanvas *c01 = new TCanvas("c01","c01",2100,2000);
	TPad *p01a = new TPad("p01a","p01a",0.05,0.27,0.95,1.0);
	TPad *p01b = new TPad("p01b","p01b",0.05,0.10,0.95,0.315);
	p01a->Draw();
	p01b->Draw();
	p01a->cd();
	p01a->SetLeftMargin(0.11);
	p01a->SetRightMargin(0.1);
	p01a->SetBottomMargin(0.11);
	//axis,log scale and range setting functions must be called after all plotOn functions being called
	TAxis* xaxis = frame->GetXaxis();
	TAxis* yaxis = frame->GetYaxis();
	xaxis->SetTitle("");
	xaxis->SetTitleOffset(1.2);
	yaxis->SetTitle("Events / 20 GeV");
	yaxis->SetTitleOffset(1.35);
	//yaxis->SetRangeUser(0,15000);
	xaxis->SetRangeUser(i*0.75,i*1.25);
	frame->Draw();
	CMS_lumi(p01a,iPeriod,iPos);
	TLegend *l =  new TLegend(0.60,0.75,0.88,0.88);
	l->AddEntry(frame->findObject("data"),"2017 Signal MC","lep");
	l->AddEntry(frame->findObject("fit"),"Signal fit "+fun_name,"l");
	l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
	//l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
	l->Draw("same");

	p01b->cd();
	p01b->SetLeftMargin(0.11);
	p01b->SetRightMargin(0.1);
	p01b->SetTopMargin(0.2);
	p01b->SetBottomMargin(0.25);
	RooHist* hpull = frame->pullHist();
	RooPlot* pull_frame = m->frame();
	pull_frame->addPlotable(hpull,"B") ;
	hpull->SetMarkerStyle(8);
	hpull->SetMarkerSize(0);
	hpull->SetLineWidth(0);
	hpull->SetFillColor(kAzure+1);
	hpull->SetFillColorAlpha(kAzure+1,0.7);
	xaxis = pull_frame->GetXaxis();
	yaxis = pull_frame->GetYaxis();
	xaxis->SetTitle("M_{j#gamma} (GeV)");
	yaxis->SetTitle("#frac{MC-fit}{#sigma_{stat}}");
	yaxis->SetTitleOffset(0.4);
	yaxis->SetRangeUser(-5,5);
	xaxis->SetLabelSize(0.12);
	xaxis->SetTitleSize(0.12);
	yaxis->SetLabelSize(0.12);
	yaxis->SetTitleSize(0.12);
	yaxis->SetNdivisions(5);
	xaxis->SetRangeUser(i*0.75,i*1.25);
	p01b->SetGrid();
	pull_frame->Draw();
	p01b->Update();
	TString plotname = "CBGaus"+mass_str+".png";
	c01->Print(plotname);
	
	// --- Output root file ---
	RooWorkspace *w = new RooWorkspace("w","w");
	w->import(*m);
	w->import(datah,Rename("signal_MC"));
	w->import(*com_model);
	TString outname = mass_str+"N-shapes-Unbinned-CBGaus"+".root";
	w->writeToFile(outname);
	}

}
