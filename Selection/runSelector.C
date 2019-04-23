void runSelector(TString fname)
{
  TFile* infile(fname);
  TDirectory *folder;
  folder = (TDirectory*)infile->Get("ntuplizer");
  TTree* tree = (TTree*)folder->Get("tree");
  tree->Process("Selector.C");
}
