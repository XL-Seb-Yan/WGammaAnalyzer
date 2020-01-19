void runtreeplot()
{
  gROOT->SetBatch(1);
  TFile infile("SinglePhoton2017_WGamma_full_full_Jan12.root");
  TTree* tree = (TTree*)infile.Get("Events");
  tree->Process("treeplot.C");
}
