void runtreeplot()
{
  gROOT->SetBatch(1);
  TFile infile("BackgroundMCCombined_Wwindow_presel_weightedTo41p54_fitData.root");
  TTree* tree = (TTree*)infile.Get("Events");
  tree->Process("treeplot.C");
}
