void chi2ana(){
  using namespace RooFit;
  gROOT->SetBatch(1);
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/chi2check/sideband-shapes-Unbinned-dijet2chi2test.root");
  RooWorkspace *w = (RooWorkspace*)file.Get("w");
  w->Print();
  TString fun_name = "dijet2";
  RooAbsPdf *model = w->pdf("dijet2");
  RooAbsData *data = w->data("data_sideband");
  RooRealVar *x = w->var("sys_invmass");
  gStyle->SetOptStat(111111);
  RooPlot *frame = x->frame(Title("Data Sideband"));
  x->setBins(30);
  RooDataHist datah("dh","binned data",RooArgSet(*x),*data);
  datah.plotOn(frame,RooFit::Name("datah"),Binning(30,600,3000),DataError(RooAbsData::Poisson));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  datah.plotOn(frame,RooFit::Name("datah"),Binning(30,600,3000),DataError(RooAbsData::Poisson));
  
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
  RooChi2Var chi2 ("chi2", "chi2", *model,datah,DataError(RooAbsData::Poisson));//Default: SumW2
  cout<<"chi2 from RooRealVar with possion is: "<<chi2.getVal()<<endl;
  cout<<"reduced chi2/ndof is: "<<frame->chiSquare(fun_name,"datah", nfloparam)<<endl;
  cout<<"chi2/ndof is: "<<frame->chiSquare(fun_name,"datah")<<endl;
  double n_non0 = frame->chiSquare(fun_name,"datah", nfloparam)*2/(frame->chiSquare(fun_name,"datah", nfloparam)-frame->chiSquare(fun_name,"datah"));
  cout<<"Estimated number of bins from chi2/ndof: "<<n_non0<<endl;
  cout<<"compare: reduced chi2: "<<frame->chiSquare(fun_name,"datah", nfloparam)*(n_non0-nfloparam)<<" non-reduced chi2: "<<frame->chiSquare(fun_name,"datah")*n_non0<<endl;
  frame->Print("V");

  TCanvas *c01 = new TCanvas("c01","c01",1200,900);
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{j#gamma} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 20 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.1,10000);
  //frame->SetMaximum(600);
  //frame->SetMinimum(0);
  c01->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.6,0.7,0.8,0.78);
  l->AddEntry(frame->findObject(fun_name),"Data SB fit "+fun_name,"l");
  l->Draw("same");
  c01->Print(fun_name+".png");
}
