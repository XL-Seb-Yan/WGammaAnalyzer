#define fun_type 9
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
  RooRealVar *x = new RooRealVar("sys_invmass","invmass",600,2500,""); //the name "sys_invmass" will be used by RooDataSet to import data

  //--- background PDF ---
#if fun_type == 1
  //-----------------------------dijet 1-----------------------------------
  TString fun_name = "dijet1";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0))",RooArgList(*x,*p0));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 2
  //-----------------------------dijet 2-----------------------------------
  TString fun_name = "dijet2";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 3
  //-----------------------------dijet3-----------------------------------
  TString fun_name = "dijet3";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)+p2*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
  #elif fun_type == 4
  //-----------------------------dijet4----------------------------------- Problematic
  TString fun_name = "dijet4";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-10,10,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-10,10,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-10,10,""); //-0.799413 +- 0.560435
  RooRealVar *p3 = new RooRealVar("p3","p_3",0,-10,10,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(sys_invmass/13000,p0+p1*log(sys_invmass/13000)+p2*pow(log(sys_invmass/13000),2)+p3*pow(log(sys_invmass/13000),3)))",RooArgList(*x,*p0,*p1,*p2,*p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 5
  //-----------------------------ATLAS1-----------------------------------
  TString fun_name = "ATLAS1";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1))",RooArgList(*x,*p0,*p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 6 
  //-----------------------------ATLAS2-----------------------------------
  TString fun_name = "ATLAS2";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1,*p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 7
  //-----------------------------ATLAS3-----------------------------------
  TString fun_name = "ATLAS3";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-pow(sys_invmass/13000,1/3),p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)+p3*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2,*p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 8
  //-----------------------------VVdijet1-----------------------------------
  TString fun_name = "VVdijet1";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1))",RooArgList(*x,*p0,*p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 9
  //-----------------------------VVdijet2-----------------------------------
  TString fun_name = "VVdijet2";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)))",RooArgList(*x,*p0,*p1,*p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 10
  //-----------------------------VVdijet3-----------------------------------
  TString fun_name = "VVdijet3";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"(pow(1-sys_invmass/13000,p0)/pow(sys_invmass/13000,p1+p2*log(sys_invmass/13000)+p3*pow(log(sys_invmass/13000),2)))",RooArgList(*x,*p0,*p1,*p2,*p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 11
  //-----------------------------exp1-----------------------------------
  TString fun_name = "expow1";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1*pow(sys_invmass/13000,2))",RooArgList(*x,*p0,*p1));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 12
  //-----------------------------exp2-----------------------------------
  TString fun_name = "expow2";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",-60,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1*pow(sys_invmass/13000,2)+p2*pow(sys_invmass/13000,3))",RooArgList(*x,*p0,*p1,*p2));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#elif fun_type == 13
  //-----------------------------exp3-----------------------------------
  TString fun_name = "expow3";
  RooRealVar *p0 = new RooRealVar("p0","p_0",0,-1000,1000,""); //-10.5471 +- 3.06516
  RooRealVar *p1 = new RooRealVar("p1","p_1",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p2 = new RooRealVar("p2","p_2",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooRealVar *p3 = new RooRealVar("p3","p_3",0,-1000,1000,""); //-0.799413 +- 0.560435
  RooGenericPdf *model = new RooGenericPdf(fun_name,"pow(sys_invmass/13000,p0)*exp(p1*pow(sys_invmass/13000,2)+p2*pow(sys_invmass/13000,3)+p3*pow(sys_invmass/13000,4))",RooArgList(*x,*p0,*p1,*p2,*p3));
  RooRealVar *bkg_norm = new RooRealVar("bkg_norm","bkg_norm",5000,0,100000,""); 
  RooAddPdf *ex_model = new RooAddPdf(fun_name+"extended",fun_name+"extended",RooArgList(*model),RooArgList(*bkg_norm));
#endif


  // --- Import unBinned dataset ---
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/Singlephoton2017_Wsideband_full_finalcut.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(*x);
  RooArgSet imargset(imarglist);
  RooDataSet data("Data sideband","Data sideband",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset
  RooRealVar weight("weight","weight",double(5653)/double(1076));
  RooDataSet data_norm("Data sideband norm.","Data sideband (nomalized to yield in W band)",&data,RooArgSet(*x,weight),"",weight.GetName());
  cout<<"number of data: "<<data.numEntries()<<endl;
 
#define isNorm 2 //1: Norm 2: Orig
  
  // --- Perform extended ML fit of composite PDF to toy data ---
#if isNorm == 1
    RooFitResult *r = model->fitTo(data_norm,Range(600,2500),RooFit::Minimizer("Minuit2"),SumW2Error(false),Save()); //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)
#else
    RooFitResult *r = model->fitTo(data,Range(600,2500),RooFit::Minimizer("Minuit2"),Save()); //SumW2Error(false) for weighted data, see how to choose this with same calling without SumW2Error(false)
#endif
  // --- plot for chi2 calculation and visualization ---
  x->setBins(95); //fit is unbinned but chi2 is calculated by binning data with this value
  RooPlot *frame = x->frame();
#if isNorm == 1
    frame->SetTitle("Data Sideband nomalized to W band");
    RooDataHist datah("dh","binned data",RooArgSet(*x),data_norm);
    datah.plotOn(frame,RooFit::Name("datah"),Binning(95,600,2500),DataError(RooAbsData::SumW2)); //for weighted data
#else
    frame->SetTitle("Data Sideband");
    RooDataHist datah("dh","binned data",RooArgSet(*x),data);
    datah.plotOn(frame,RooFit::Name("datah"),Binning(95,600,2500),DataError(RooAbsData::Poisson)); //for unweighted data
#endif
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  model->plotOn(frame,VisualizeError(*r,2,kFALSE),FillColor(kYellow),LineColor(0),RooFit::Name("err2"));
  model->plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kGreen),LineColor(0),RooFit::Name("err1"));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
#if isNorm == 1
    datah.plotOn(frame,RooFit::Name("datah"),Binning(95,600,2500),DataError(RooAbsData::SumW2)); //for weighted data
#else
    datah.plotOn(frame,RooFit::Name("datah"),Binning(95,600,2500),DataError(RooAbsData::Poisson)); //for unweighted data
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

   //chi2 calculation is problematic, waiting for anwser on Oct.18 2019
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
  /*
  // --- Output root file ---
  if(isNorm){
    RooWorkspace *w = new RooWorkspace("w","w");
    w->import(*x);
    w->import(data_norm,Rename("data_norm_sideband"));
    w->import(*model);
    w->writeToFile("sideband_norm-shapes-Unbinned-"+fun_name+".root");
  }
  else{
    RooWorkspace *w = new RooWorkspace("w","w");
    w->import(*x);
    w->import(data,Rename("data_sideband"));
    w->import(*model);
    w->writeToFile("sideband-shapes-Unbinned-"+fun_name+".root");
  }
  */

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
