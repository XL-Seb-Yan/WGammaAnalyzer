void ReadWorkspace(){
  TFile file("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/RooFit/sideband-shapes-Unbinned-dijet2.root");
  RooWorkspace *w = (RooWorkspace*)file.Get("w");
  w->Print();
  RooAbsPdf *pdf = w->pdf("dijet2");
  RooRealVar *p0 = w->var("p0");
  RooRealVar *p1 = w->var("p1");
  cout<<p0->getVal()<<endl;
  cout<<p1->getVal()<<endl;
}
