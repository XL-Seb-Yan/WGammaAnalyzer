void chi2ana(){
  using namespace RooFit;
  gROOT->SetBatch(1);
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/sideband-shapes-Unbinned-dijet2chi2test.root");
  RooWorkspace *w = (RooWorkspace*)file.Get("w");
  w->Print();
  TString fun_name = "dijet2";
  RooAbsPdf *model = w->pdf("dijet2");
  RooAbsData *data = w->data("data_sideband");
  RooRealVar *x = w->var("sys_invmass");
  gStyle->SetOptStat(111111);
  RooPlot *frame = x->frame(Title("Data Sideband"));
  data->plotOn(frame,RooFit::Name("data"),Binning(120,600,3000),DataError(RooAbsData::Poisson));
  model->plotOn(frame,LineStyle(kDashed),RooFit::Name(fun_name));
  data->plotOn(frame);
  RooDataHist datah("dh","binned data",RooArgSet(*x),*data);
  cout<<"data bins: "<<datah.numEntries()<<endl;
  RooArgSet *flparams = model->getParameters(*x);
  int nfloparam = (flparams->selectByAttrib("Constant",kFALSE))->getSize();
  cout<<"# of floating params: "<<nfloparam<<endl;
  RooChi2Var chi2 ("chi2", "chi2", *model,datah,DataError(RooAbsData::Poisson));//Default: SumW2
  TString chi2txt = "Chi2/Ndof: "+std::to_string(frame->chiSquare(fun_name,"data",nfloparam));
  RooAbsReal *chi2Exp = model->createChi2(datah,DataError(RooAbsData::Expected));//Using RooChi2Var with default: Expected error
  cout<<"chi2 with possion is: "<<chi2.getVal()<<endl;
  cout<<"chi2 with expected is: "<<chi2Exp->getVal()<<endl;
  cout<<"reduced chi2 is: "<<frame->chiSquare(fun_name,"data", nfloparam)<<endl;
  cout<<"chi2 is: "<<frame->chiSquare(fun_name,"data")<<endl;
  cout<<"compare: reduced: "<<frame->chiSquare(fun_name,"data", nfloparam)*(datah.numEntries()-nfloparam)<<" non-reduced: "<<frame->chiSquare(fun_name,"data")*datah.numEntries()<<endl;

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
