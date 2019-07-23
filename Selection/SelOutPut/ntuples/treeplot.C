#define treeplot_cxx
// The class definition in treeplot.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("treeplot.C")
// root> T->Process("treeplot.C","some options")
// root> T->Process("treeplot.C+")
//


#include "treeplot.h"
#include <TROOT.h> 
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TAxis.h>


Int_t count1 = 0;
Int_t count2 = 0;
TH1 *hist1 = new TH1F("1","pt_{#gamma}",50,0,2400);
TH1 *hist2 = new TH1F("2","eta_{#gamma}",50,-5,5);
TH1 *hist3 = new TH1F("3","pt_{j}",50,0,2400);
TH1 *hist4 = new TH1F("4","eta_{j}",50,-5,5);
TH1 *hist5 = new TH1F("5","E_{j}",50,0,2400);
TH1 *hist6 = new TH1F("6","masssoftdrop_{j}",60,50,110);
TH1 *hist7 = new TH1F("7","tau21_{j}",50,0,1);
TH1 *hist8 = new TH1F("8","cos(#theta*)_{p}",50,0,1);
TH1 *hist9 = new TH1F("9","pt/M",50,0,1);
TH1 *hist10 = new TH1F("10","invariant mass",50,600,1000);
TH1 *hist11 = new TH1F("11","seperation",50,0,8);

void treeplot::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
  gROOT->SetBatch(1);

   TString option = GetOption();
}

void treeplot::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t treeplot::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);

   if(entry%1000==0)
     cout<<"Processing "<<entry<<endl;
   
   if((*sys_invmass > 600 && *sys_invmass < 1000) && (*ak8puppijet_masssoftdropcorr > 90 && *ak8puppijet_masssoftdropcorr < 105)){
     hist1->Fill(*photon_pt);
     hist2->Fill(*photon_eta);
     hist3->Fill(*ak8puppijet_pt);
     hist4->Fill(*ak8puppijet_eta);
     hist5->Fill(*ak8puppijet_e);
     hist6->Fill(*ak8puppijet_masssoftdropcorr);
     hist7->Fill(*ak8puppijet_tau21);
     hist8->Fill(*sys_costhetastar);
     hist9->Fill(*sys_ptoverm);
     hist10->Fill(*sys_invmass);
     hist11->Fill(*sys_seperation);
   }
   

   return kTRUE;
}

void treeplot::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void treeplot::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  hist1->SetLineColor(8);
  hist2->SetLineColor(8);
  hist3->SetLineColor(8);
  hist4->SetLineColor(8);
  hist5->SetLineColor(8);
  hist6->SetLineColor(8);
  hist7->SetLineColor(8);
  hist8->SetLineColor(8);
  hist9->SetLineColor(8);
  hist10->SetLineColor(8);
  hist11->SetLineColor(8);
  

  gStyle->SetOptStat(0);
  
  //Non stacked plots

  TLegend *legend = new TLegend(0.65,0.8,0.9,0.9);

  TCanvas *c01 = new TCanvas("c01","pt_{#gamma}",1200,900);
  TAxis *xaxis = hist1->GetXaxis();
  TAxis *yaxis = hist1->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c01->SetLogy();
  c01->cd();
  hist1->SetLineWidth(2);
  hist1->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist1,"2017 SinglePhoton C","f");
  //legend->Draw();
  c01->Print("p_pt.png");

  TCanvas *c02 = new TCanvas("c02","eta_{#gamma}",1200,900);
  xaxis = hist2->GetXaxis();
  yaxis = hist2->GetYaxis();
  xaxis->SetTitle("eta_{#gamma}");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c02->SetLogy();
  c02->cd();
  hist2->SetLineWidth(2);
  hist2->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist2,"2017 signal MC , pass EleVeto","f");
  //legend->Draw();
  c02->Print("p_eta.png");

  TCanvas *c03 = new TCanvas("c03","pt AK8Jet",1200,900);
  xaxis = hist3->GetXaxis();
  yaxis = hist3->GetYaxis();
  xaxis->SetTitle("pt AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  //yaxis->SetRangeUser(0,2000);
  c03->SetLogy();
  c03->cd();
  hist3->SetLineWidth(2);
  hist3->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist5,"2017 signal MC ","f");
  //legend->Draw();
  c03->Print("j_pt.png");

  TCanvas *c04 = new TCanvas("c04","eta AK8Jet",1200,900);
  xaxis = hist4->GetXaxis();
  yaxis = hist4->GetYaxis();
  xaxis->SetTitle("eta AK8Jet");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c04->SetLogy();
  c04->cd();
  hist4->SetLineWidth(2);
  hist4->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist4,"2017 signal MC ","f");
  //legend->Draw();
  c04->Print("j_eta.png");

  TCanvas *c05 = new TCanvas("c05","E AK8Jet",1200,900);
  xaxis = hist5->GetXaxis();
  yaxis = hist5->GetYaxis();
  xaxis->SetTitle("E AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c05->SetLogy();
  c05->cd();
  hist5->SetLineWidth(2);
  hist5->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist5,"2017 signal MC ","f");
  //legend->Draw();
  c05->Print("j_e.png");

  TCanvas *c06 = new TCanvas("c06","mass softdrop AK8Jet",1200,900);
  xaxis = hist6->GetXaxis();
  yaxis = hist6->GetYaxis();
  xaxis->SetTitle("mass softdrop AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  //yaxis->SetRangeUser(0,400);
  c06->SetLogy();
  c06->cd();
  hist6->SetLineWidth(2);
  hist6->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist6,"2017 signal MC ","f");
  //legend->Draw();
  c06->Print("j_masssd.png");

  TCanvas *c07 = new TCanvas("c07","tau21 AK8Jet",1200,900);
  xaxis = hist7->GetXaxis();
  yaxis = hist7->GetYaxis();
  xaxis->SetTitle("tau21 AK8Jet");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c07->SetLogy();
  c07->cd();
  hist7->SetLineWidth(2);
  hist7->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist7,"2017 signal MC ","f");
  //legend->Draw();
  c07->Print("j_tau21.png");

  TCanvas *c08 = new TCanvas("c8","cos AK8Jet",1200,900);
  xaxis = hist8->GetXaxis();
  yaxis = hist8->GetYaxis();
  xaxis->SetTitle("cos(#theta*)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c08->SetLogy();
  c08->cd();
  hist8->SetLineWidth(2);
  hist8->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist8,"2017 signal MC ","f");
  //legend->Draw();
  c08->Print("s_cos.png");

  TCanvas *c09 = new TCanvas("c9","ptm AK8Jet",1200,900);
  xaxis = hist9->GetXaxis();
  yaxis = hist9->GetYaxis();
  xaxis->SetTitle("pt/M");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c09->SetLogy();
  c09->cd();
  hist9->SetLineWidth(2);
  hist9->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist9,"2017 signal MC ","f");
  //legend->Draw();
  c09->Print("s_ptm.png");

  TCanvas *c10 = new TCanvas("c10","invmass",1200,900);
  xaxis = hist10->GetXaxis();
  yaxis = hist10->GetYaxis();
  xaxis->SetTitle("pt/M");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c10->SetLogy();
  c10->cd();
  hist10->SetLineWidth(2);
  hist10->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist10,"2017 signal MC ","f");
  //legend->Draw();
  c10->Print("s_invmass.png");

  TCanvas *c11 = new TCanvas("c11","seperation",1200,900);
  xaxis = hist11->GetXaxis();
  yaxis = hist11->GetYaxis();
  xaxis->SetTitle("#Delta R");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c11->SetLogy();
  c11->cd();
  hist11->SetLineWidth(2);
  hist11->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist11,"2017 signal MC ","f");
  //legend->Draw();
  c11->Print("s_seperation.png");

}
