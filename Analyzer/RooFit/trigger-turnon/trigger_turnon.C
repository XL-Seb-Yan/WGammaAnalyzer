#include <TMath.h>
#include <TLegend.h>
void trigger_turnon(int seed=37)
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37); 

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("m","m",400,2000,""); //the name "m" will be used by RooDataSet to import data, normalization range is 600-3500 but plot range can be defined to like 600-3000

  //--- background PDF ---
  TString fun_name = "ErfExp";
  RooRealVar* p0 = new RooRealVar("p0","p0",-0.0005,-0.1,-0.0001);
  RooRealVar* p1 = new RooRealVar("p1","p1",0,-0.1,1);
  RooRealVar* offset = new RooRealVar("offset","offset",600,0,1000);
  RooRealVar* width  = new RooRealVar("width","width",60,0,100);
  RooGenericPdf * model = new RooGenericPdf("ErfExp","ErfExp","TMath::Exp(p0*m+p1*pow(m,2))*(1.+TMath::Erf((m-offset)/width))/2.",RooArgList(*x, *p0, *p1, *offset, *width));

  // --- Import Binned dataset ---
  float s_mass;
  float x_weight = 1;
  float x_kfactor = 1;
  float x_puweight = 1;
  TH1F MChist("MC","MC",40,400,2000);
  TFile file("/afs/cern.ch/user/x/xuyan/WGProj/PROD17/DATA/2017/fullcut/BkgMC_postproc_WGamma_full_68-94_fullcut_jmcorr_Mar17.root");
  TTree *tree = (TTree*)file.Get("Events");
  tree->SetBranchAddress("m", &s_mass);
  tree->SetBranchAddress("xsec_weight", &x_weight);
  tree->SetBranchAddress("xsec_kfactor", &x_kfactor);
  tree->SetBranchAddress("xsec_puweight", &x_puweight);
  for (int ievt = 0; ievt<tree->GetEntries();ievt++) {
    tree->GetEntry(ievt);
    MChist.Fill(s_mass, x_weight*x_kfactor*x_puweight);
  }
  RooDataHist datah("Background MC (W band)","Background MC (W band)",RooArgSet(*x),&MChist);
  cout<<"number of weighted entries: "<<datah.sum(false)<<endl;
  
  // --- Perform extended ML fit of composite PDF to toy data ---
    RooFitResult *r = model->fitTo(datah,Range(400,2000),RooFit::Minimizer("Minuit2"),SumW2Error(false),Save());  //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)
    
  // --- plot for chi2 calculation and visualization ---
  //x->setBins(100); //fit is unbinned but chi2 is calculated by binning data with this value
  //x->setBins(95); //fit is unbinned but chi2 is calculated by binning data with this value
  RooPlot *frame = x->frame();
  frame->SetTitle("Background MC (W band)");
  datah.plotOn(frame,RooFit::Name("datah"),Binning(40,400,2000),DataError(RooAbsData::SumW2)); //for weighted data
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  datah.plotOn(frame,RooFit::Name("datah"),Binning(40,400,2000),DataError(RooAbsData::SumW2)); //for weighted data
  
  // --- Visualization ---
  gStyle->SetOptStat(111111);
  TCanvas *c01 = new TCanvas("c01","c01",1200,900);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 20 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0,2500);
  xaxis->SetLimits(400,1500);
  //c01->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.6,0.7,0.8,0.78);
  l->AddEntry(frame->findObject(fun_name),"MC signal region fit "+fun_name,"l");
  l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
  c01->Print(fun_name+".png");
  c01->Print(fun_name+".pdf");
  c01->Print(fun_name+".svg");
  
  // TCanvas *c02 = new TCanvas("c02","c02",1200,300);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  // c02->cd();
  // RooHist* hpull = frame->pullHist();
  // RooPlot* pull_frame = x->frame();
  // pull_frame->addPlotable(hpull,"P") ;
  // hpull->SetMarkerStyle(8);
  // xaxis = pull_frame->GetXaxis();
  // yaxis = pull_frame->GetYaxis();
  // xaxis->SetTitle("M_{j#gamma}");
  // yaxis->SetTitle("#frac{MC-fit}{#sigma_{stat}}");
  // yaxis->SetTitleOffset(0.5);
  // yaxis->SetRangeUser(-5,5);
  // xaxis->SetLabelSize(0.08);
  // xaxis->SetTitleSize(0.08);
  // xaxis->SetRangeUser(0,2000);
  // yaxis->SetLabelSize(0.08);
  // yaxis->SetTitleSize(0.08);
  // yaxis->SetNdivisions(5);
  // c02->SetBottomMargin(0.18);
  // c02->SetTopMargin(0.18);
  // c02->SetGrid();
  // pull_frame->Draw();
  // c02->Update();
  // c02->Print("pull"+fun_name+".png");
  // c02->Print("pull"+fun_name+".pdf"); 
  // c02->Print("pull"+fun_name+".svg"); 
  
  TCanvas *c03 = new TCanvas("c03","c03",1200,900);
  c03->cd();
  TF1 *f1 = new TF1("erf","2250*(1.+TMath::Erf((x-[0])/[1]))/2.",400,1500);
  f1->SetParameters(offset->getValV(), width->getValV());
  TF1 *f2 = new TF1("exp","300000*TMath::Exp([0]*x+[1]*pow(x,2))",400,1500);
  f2->SetParameters(p0->getValV(), p1->getValV());
  cout<<"Turn-on end point is: "<<f1->GetX(0.99)<<endl;
  f1->GetYaxis()->SetRangeUser(0,2500);
  f1->Draw();
  f2->Draw("SAME");
  c03->Print("pdf.png");
  c03->Print("pdf.svg");
}
