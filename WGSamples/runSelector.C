void runSelector()
{
  TFile infile("flatTuple_949_final_v6.root");
  TDirectory *folder;
  folder = (TDirectory*)infile.Get("ntuplizer");
  TTree* tree = (TTree*)folder->Get("tree");
  tree->Process("Selector.C");
}
