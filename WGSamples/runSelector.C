void runSelector()
{
  TFile infile("flatTuple_muon_johnB.root");
  TDirectory *folder;
  folder = (TDirectory*)infile.Get("ntuplizer");
  TTree* tree = (TTree*)folder->Get("tree");
  tree->Process("Selector.C");
}
