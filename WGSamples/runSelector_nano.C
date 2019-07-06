void runSelector()
{
  TFile infile("MC1600_nanoAOD.root");
  //TDirectory *folder;
  //folder = (TDirectory*)infile.Get("ntuplizer");
  //TTree* tree = (TTree*)folder->Get("tree");
  TTree* tree = (TTree*)infile.Get("Events");
  tree->Process("Selector.C");
}
