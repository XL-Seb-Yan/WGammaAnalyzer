#define fun_type 7
#define isNorm 2 //1: Norm 2: Orig
#include <TMath.h>
#include <TLegend.h>
void make_sideband_shapes(int seed=37)
{
  gErrorIgnoreLevel = kInfo;
  gROOT->SetBatch(1);
  using namespace RooFit;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRandom::randomGenerator()->SetSeed(37); 

  // --- Create obervable --- 
  RooRealVar *x = new RooRealVar("m","m",600,3500,""); //the name "m" will be used by RooDataSet to import data, normalization range is 600-3500 but plot range can be defined to like 600-3000

  //--- background PDF ---
#if fun_type == 1
  //-----------------------------dijet 2-----------------------------------
  TString fun_name = "dijet2";
  RooRealVar *dijet2_p0 = new RooRealVar("dijet2_p0","dijet2_p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *dijet2_p1 = new RooRealVar("dijet2_p1","dijet2_p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(m/13000,dijet2_p0+dijet2_p1*log(m/13000)))",RooArgList(*x,*dijet2_p0,*dijet2_p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 2
  //-----------------------------ATLAS1-----------------------------------
  TString fun_name = "ATLAS1";
  RooRealVar *ATLAS1_p0 = new RooRealVar("ATLAS1_p0","ATLAS1_p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *ATLAS1_p1 = new RooRealVar("ATLAS1_p1","ATLAS1_p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(m/13000,1/3),ATLAS1_p0)/pow(m/13000,ATLAS1_p1))",RooArgList(*x,*ATLAS1_p0,*ATLAS1_p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 3
  //-----------------------------VVdijet2-----------------------------------
  TString fun_name = "VVdijet2";
  /*
  //Normalized fit
  RooRealVar *VVdijet2_p0 = new RooRealVar("VVdijet2_p0","VVdijet2_p_0",-15,-20,20,""); //-10.5471 +- 3.06516
  RooRealVar *VVdijet2_p1 = new RooRealVar("VVdijet2_p1","VVdijet2_p_1",15,-20,20,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet2_p2 = new RooRealVar("VVdijet2_p2","VVdijet2_p_2",0,-20,20,""); //-0.799413 +- 0.560435
  */
  RooRealVar *VVdijet2_p0 = new RooRealVar("VVdijet2_p0","VVdijet2_p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *VVdijet2_p1 = new RooRealVar("VVdijet2_p1","VVdijet2_p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet2_p2 = new RooRealVar("VVdijet2_p2","VVdijet2_p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-m/13000,VVdijet2_p0)/pow(m/13000,VVdijet2_p1+VVdijet2_p2*log(m/13000)))",RooArgList(*x,*VVdijet2_p0,*VVdijet2_p1,*VVdijet2_p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 4
  //-----------------------------exp-----------------------------------
  TString fun_name = "exp";
  RooRealVar *exp_p0 = new RooRealVar("exp_p0","exp_p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *exp_p1 = new RooRealVar("exp_p1","exp_p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"exp(exp_p0*(m/13000)+exp_p1*pow(m/13000,2))",RooArgList(*x,*exp_p0,*exp_p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 5
  //-----------------------------dijet 3-----------------------------------
  TString fun_name = "dijet3";
  RooRealVar *dijet3_p0 = new RooRealVar("dijet3_p0","dijet3_p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *dijet3_p1 = new RooRealVar("dijet3_p1","dijet3_p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *dijet3_p2 = new RooRealVar("dijet3_p2","dijet3_p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(m/13000,dijet3_p0+dijet3_p1*log(m/13000)+dijet3_p2*pow(log(m/13000),2)))",RooArgList(*x,*dijet3_p0,*dijet3_p1,*dijet3_p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 6
  //-----------------------------ATLAS2-----------------------------------
  TString fun_name = "ATLAS2";
  RooRealVar *ATLAS2_p0 = new RooRealVar("ATLAS2_p0","ATLAS2_p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *ATLAS2_p1 = new RooRealVar("ATLAS2_p1","ATLAS2_p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *ATLAS2_p2 = new RooRealVar("ATLAS2_p2","ATLAS2_p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(m/13000,1/3),ATLAS2_p0)/pow(m/13000,ATLAS2_p1+ATLAS2_p2*log(m/13000)))",RooArgList(*x,*ATLAS2_p0,*ATLAS2_p1,*ATLAS2_p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 7
  //-----------------------------VVdijet3-----------------------------------
  TString fun_name = "VVdijet3";
  /*
  //Normalized fit
  RooRealVar *VVdijet3_p0 = new RooRealVar("VVdijet3_p0","VVdijet3_p_0",-15,-20,20,""); //-10.5471 +- 3.06516
  RooRealVar *VVdijet3_p1 = new RooRealVar("VVdijet3_p1","VVdijet3_p_1",15,-20,20,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet3_p2 = new RooRealVar("VVdijet3_p2","VVdijet3_p_2",0,-20,20,""); //-0.799413 +- 0.560435
  */
  RooRealVar *VVdijet3_p0 = new RooRealVar("VVdijet3_p0","VVdijet3_p_0",0,-100,100,""); //-10.5471 +- 3.06516
  RooRealVar *VVdijet3_p1 = new RooRealVar("VVdijet3_p1","VVdijet3_p_1",0,-10,10,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet3_p2 = new RooRealVar("VVdijet3_p2","VVdijet3_p_2",0,-1,1,""); //-0.799413 +- 0.560435
  RooRealVar *VVdijet3_p3 = new RooRealVar("VVdijet3_p3","VVdijet3_p_3",0,-1,1,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-m/13000,VVdijet3_p0)/pow(m/13000,VVdijet3_p1+VVdijet3_p2*log(m/13000)+VVdijet3_p3*pow(log(m/13000),2)))",RooArgList(*x,*VVdijet3_p0,*VVdijet3_p1,*VVdijet3_p2,*VVdijet3_p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#endif

  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/fullcutdataset/SinglePhoton2017_WGamma_Wsideband_full_finalcut.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooDataSet data("Data sideband","Data sideband",RooArgSet(*x),Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset
  RooRealVar weight("weight","weight",double(5653)/double(1076));
  RooDataSet data_norm("Data sideband norm.","Data sideband (nomalized to yield in W band)",&data,RooArgSet(*x,weight),"",weight.GetName());
  cout<<"number of data: "<<data.numEntries()<<endl;
  
  // --- Perform extended ML fit of composite PDF to toy data ---
#if isNorm == 1
    RooFitResult *r = model->fitTo(data_norm,Range(600,3500),RooFit::Minimizer("Minuit2"),SumW2Error(false),Save()); //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)
#else
    RooFitResult *r = model->fitTo(data,Range(600,3500),RooFit::Minimizer("Minuit2"),Save()); //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)
#endif
    
  // --- plot for chi2 calculation and visualization ---
  x->setBins(145); //fit is unbinned but chi2 is calculated by binning data with this value
  RooPlot *frame = x->frame();
#if isNorm == 1
    frame->SetTitle("Data Sideband nomalized to W band");
    RooDataHist datah("dh","binned data",RooArgSet(*x),data_norm);
    datah.plotOn(frame,RooFit::Name("datah"),Binning(145,600,3500),DataError(RooAbsData::SumW2)); //for weighted data
#else
    frame->SetTitle("Data Sideband");
    RooDataHist datah("dh","binned data",RooArgSet(*x),data);
    datah.plotOn(frame,RooFit::Name("datah"),Binning(145,600,3500),DataError(RooAbsData::Poisson)); //for unweighted data
#endif
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
#if isNorm == 1
    datah.plotOn(frame,RooFit::Name("datah"),Binning(145,600,3500),DataError(RooAbsData::SumW2)); //for weighted data
#else
    datah.plotOn(frame,RooFit::Name("datah"),Binning(145,600,3500),DataError(RooAbsData::Poisson)); //for unweighted data
#endif

    frame->Print("V");
  RooAbsReal* nll = NULL;
#if isNorm == 1
    nll = model->createNLL(data_norm);
#else
    nll = model->createNLL(data);
#endif
  
  cout<<"NLL is: "<<nll->getVal()<<endl;

  cout<<datah<<endl;

  cout<<"data bins: "<<datah.numEntries()<<endl;
  int n_0 = 0;
  for(int i=0; i<datah.numEntries(); i++){
    datah.get(i) ;
    if(datah.weight() == 0)
      n_0++;
  }
  cout<<"Number of empty bins: "<<n_0<<endl;
  
  RooArgSet *flparams = model->getParameters(*x);
  int nfloparam = (flparams->selectByAttrib("Constant",kFALSE))->getSize();
  cout<<"# of floating params: "<<nfloparam<<endl;

  // Chi2 problem solved:
  // Unweighted data:
  // RooChi2Var uses RooAbsData::Possion to obtain Possion interval (which can treat 0 entry bins), only skip bins with both PDF and data = 0
  // chi2square() obtains chi2 by skipping 0 entry bins, divide by numbr of non-0 entry bins
  // Weighted data:
  // RooChi2Var uses RooAbsData::SumW2, fail when arrive 0-entry bins
  // chi2square() obtains chi2 by skipping 0 entry bins, divide by numbr of non-0 entry bins
  RooChi2Var chi2 ("chi2", "chi2", *model,datah,DataError(RooAbsData::SumW2));//Default: SumW2
  TString chi2txt = "Chi2/Ndof: "+std::to_string(frame->chiSquare(fun_name,"datah",nfloparam));
  
  cout<<"chi2 from RooRealVar with possion is: "<<chi2.getVal()<<endl;
  cout<<"reduced chi2/ndof is: "<<frame->chiSquare(fun_name,"datah", nfloparam)<<endl;
  cout<<"chi2/ndof is: "<<frame->chiSquare(fun_name,"datah")<<endl;
  double n_non0 = frame->chiSquare(fun_name,"datah", nfloparam)*2/(frame->chiSquare(fun_name,"datah", nfloparam)-frame->chiSquare(fun_name,"datah"));
  cout<<"Estimated number of bins from chi2/ndof: "<<n_non0<<endl;
  cout<<"compare: reduced chi2: "<<frame->chiSquare(fun_name,"datah", nfloparam)*(n_non0-nfloparam)<<" non-reduced chi2: "<<frame->chiSquare(fun_name,"datah")*n_non0<<endl;
  frame->Print("V");
  
  TString NLLtxt = "NLL: "+std::to_string(nll->getVal());
  TLatex *NLLlax = new TLatex(1830,80,NLLtxt);
  NLLlax->SetTextSize(0.025);
  NLLlax->SetTextColor(kBlack);
  frame->addObject(NLLlax);

  // access residuals
  double SSR = 0;
  RooCurve* curve = (RooCurve*)frame->findObject(fun_name,RooCurve::Class());
  RooHist* hist = (RooHist*)frame->findObject("datah",RooHist::Class());
  Double_t xstart,xstop,y;
  curve->GetPoint(0,xstart,y);
  curve->GetPoint(curve->GetN()-1,xstop,y);
  for(Int_t i=0 ; i<hist->GetN() ; i++) {
    Double_t x,point;
    hist->GetPoint(i,x,point);
    // Only calculate pull for bins inside curve range
    if (x<xstart || x>xstop) continue;
    Double_t yy;
    yy = point - curve->interpolate(x);
    SSR += pow(yy,2);
  }
  cout<<"SSR is: "<<SSR<<" NdoF is: "<<hist->GetN()<<endl;
  TString SSRtxt = "SSR: "+std::to_string(SSR);
  TLatex *SSRlax = new TLatex(1830,60,SSRtxt);
  SSRlax->SetTextSize(0.025);
  SSRlax->SetTextColor(kBlack);
  frame->addObject(SSRlax);

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
  yaxis->SetRangeUser(0.2,10000);
  xaxis->SetRangeUser(600,3000);
  //frame->SetMaximum(600);
  //frame->SetMinimum(0);
  c01->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.6,0.7,0.8,0.78);
#if isNorm == 1
    l->AddEntry(frame->findObject(fun_name),"Data SB norm. fit "+fun_name,"l");
#else
    l->AddEntry(frame->findObject(fun_name),"Data SB fit "+fun_name,"l");
#endif
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err1"),"Fit Error 1 #sigma","f");
  l->AddEntry(frame->findObject("err2"),"Fit Error 2 #sigma","f");
  l->Draw("same");
#if isNorm == 1
    c01->Print(fun_name+"_norm.png");
#else
    c01->Print(fun_name+".png");
#endif
  TCanvas *c02 = new TCanvas("c02","c02",1200,300);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  c02->cd();
  RooHist* hpull = frame->pullHist();
  RooPlot* pull_frame = x->frame();
  pull_frame->addPlotable(hpull,"P") ;
  hpull->SetMarkerStyle(8);
  xaxis = pull_frame->GetXaxis();
  yaxis = pull_frame->GetYaxis();
  xaxis->SetTitle("M_{j#gamma}");
  yaxis->SetTitle("#frac{data_norm-fit}{#sigma_{stat}}");
  yaxis->SetTitleOffset(0.5);
  yaxis->SetRangeUser(-5,5);
  xaxis->SetLabelSize(0.08);
  xaxis->SetTitleSize(0.08);
  xaxis->SetRangeUser(600,3000);
  yaxis->SetLabelSize(0.08);
  yaxis->SetTitleSize(0.08);
  yaxis->SetNdivisions(5);
  c02->SetBottomMargin(0.18);
  c02->SetTopMargin(0.18);
  c02->SetGrid();
  pull_frame->Draw();
  c02->Update();
#if isNorm == 1
  c02->Print("pull"+fun_name+"_norm.png"); //pull hist would not be accurated for weighted data
#else
  c02->Print("pull"+fun_name+".png"); 
#endif
  
  // --- Output root file ---
#if isNorm == 1
    RooWorkspace *w = new RooWorkspace("w","w");
    w->import(*x);
    w->import(data_norm,Rename("data_norm_sideband"));
    w->import(*model);
    w->writeToFile("sideband_norm-shapes-Unbinned-"+fun_name+".root");
#else
    RooWorkspace *w = new RooWorkspace("w","w");
    w->import(*x);
    w->import(data,Rename("data_sideband"));
    w->import(*model);
    w->writeToFile("sideband-shapes-Unbinned-"+fun_name+".root");
#endif
    
  

  /*
  // --- Perform extended ML fit of composite PDF to toy data ---
  RooFitResult *ex_r = NULL;
  if(isNorm)
    ex_r = ex_model->fitTo(data_norm,Range(600,3000),RooFit::Minimizer("Minuit2"),Extended(true),SumW2Error(false),Save());
  else
    ex_r = ex_model->fitTo(data,Range(600,3000),RooFit::Minimizer("Minuit2"),Extended(true),Save());
  cout<<"Normalization is: "<<bkg_norm->getVal()<<endl;
  */
  
}
