void runSelector()
{
  TFile infile("flatTuple_muon_full_B20.root");
  TDirectory *folder;
  folder = (TDirectory*)infile.Get("ntuplizer");
  TTree* tree = (TTree*)folder->Get("tree");
  tree->Process("Selector.C");
}
