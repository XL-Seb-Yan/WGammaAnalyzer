void runtreeplot()
{
  TFile infile("SinglePhoton2017_sideband.root");
  TTree* tree = (TTree*)infile.Get("Events");
  tree->Process("treeplot.C");
}
