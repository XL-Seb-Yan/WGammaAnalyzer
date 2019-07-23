void runtreeplot()
{
  TFile infile("SinglePhoton2017B_WGamma_select.root");
  TTree* tree = (TTree*)infile.Get("Events");
  tree->Process("treeplot.C");
}
