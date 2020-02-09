void runSelector_nano()
{
  gROOT->SetBatch(1);
  TFile infile("Signal_Nano_16_M1000N.root");
  TTree* tree = (TTree*)infile.Get("Events");
  cout<<tree<<endl;
  tree->Process("Selector_nano.C");
}
