void runtreeplot()
{
  TFile infile("SinglePhoton2017_WGamma_50105_full.root");
  TTree* tree = (TTree*)infile.Get("Events");
  tree->Process("treeplot.C");
}
