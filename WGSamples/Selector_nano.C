#define Selector_cxx
// The class definition in Selector.h has been generated automatically
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
// root> T->Process("Selector.C")
// root> T->Process("Selector.C","some options")
// root> T->Process("Selector.C+")
//


#include "Selector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TMath.h>

Int_t count1 = 0;
Int_t count2 = 0;
TH1 *hist1 = new TH1F("1","pt_{#gamma}",50,0,2400);
TH1 *hist2 = new TH1F("2","eta_{#gamma}",50,-1,1);
TH1 *hist3 = new TH1F("3","mvaID_{#gamma}",50,-1,1);
TH1 *hist4 = new TH1F("4","mvaIDv1_{#gamma}",50,-1,1);
TH1 *hist5 = new TH1F("5","pt_{j}",50,0,2400);
TH1 *hist6 = new TH1F("6","eta_{j}",50,-5,5);
TH1 *hist7 = new TH1F("7","mass_{j}",50,0,500);
TH1 *hist8 = new TH1F("8","masssoftdrop_{j}",60,50,110);
TH1 *hist9 = new TH1F("9","tau1_{j}",50,0,1);
TH1 *hist10 = new TH1F("10","tau2_{j}",50,0,1);


void Selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void Selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t Selector::Process(Long64_t entry)
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

   if(entry%100000==0)
     cout<<"Processing "<<entry<<endl;
  

   for(int i=0; i<Photon_eta.GetSize(); i++){
     if(Photon_electronVeto[i] == true && Photon_pt[i] >= 200){
       hist1->Fill(Photon_pt[i]);
       hist2->Fill(Photon_eta[i]);
       hist3->Fill(Photon_mvaID[i]);
       hist4->Fill(Photon_mvaIDV1[i]);
     }
   }
   for(int i=0; i<FatJet_eta.GetSize(); i++){
     if(FatJet_pt[i] >= 200 && (FatJet_eta[i] < 2.4 && FatJet_eta[i] > -2.4)){
       hist5->Fill(FatJet_pt[i]);
       hist6->Fill(FatJet_eta[i]);
       hist7->Fill(FatJet_mass[i]);
       hist8->Fill(FatJet_msoftdrop[i]);
       hist9->Fill(FatJet_tau1[i]);
       hist10->Fill(FatJet_tau2[i]);
     }
   }
   

   return kTRUE;
}

void Selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Selector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  gStyle->SetOptStat(0);
    //Tirgger Efficiency
  TEfficiency *Eff = new TEfficiency(*hist2, *hist3);
  
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
  legend->AddEntry(hist1,"2017 signal MC , pass EleVeto","f");
  legend->Draw();
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
  legend->Draw();
  c02->Print("p_eta.png");

  TCanvas *c03 = new TCanvas("c03","mvaID_{#gamma}",1200,900);
  xaxis = hist3->GetXaxis();
  yaxis = hist3->GetYaxis();
  xaxis->SetTitle("mvaID_{#gamma} (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c03->SetLogy();
  c03->cd();
  hist3->SetLineWidth(2);
  hist3->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist3,"2017 signal MC , pass EleVeto","f");
  legend->Draw();
  c03->Print("p_mvaID.png");

  TCanvas *c04 = new TCanvas("c04","mvaID1_{#gamma}",1200,900);
  xaxis = hist4->GetXaxis();
  yaxis = hist4->GetYaxis();
  xaxis->SetTitle("mvaID1_{#gamma} (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c04->SetLogy();
  c04->cd();
  hist4->SetLineWidth(2);
  hist4->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist4,"2017 signal MC , pass EleVeto","f");
  legend->Draw();
  c04->Print("p_mvaID1.png");

  TCanvas *c05 = new TCanvas("c05","pt AK8Jet",1200,900);
  xaxis = hist5->GetXaxis();
  yaxis = hist5->GetYaxis();
  xaxis->SetTitle("pt AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0.1,100000);
  yaxis->SetRangeUser(0,2000);
  //c05->SetLogy();
  c05->cd();
  hist5->SetLineWidth(2);
  hist5->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist5,"2017 signal MC ","f");
  legend->Draw();
  c05->Print("j_pt.png");

  TCanvas *c06 = new TCanvas("c06","eta AK8Jet",1200,900);
  xaxis = hist6->GetXaxis();
  yaxis = hist6->GetYaxis();
  xaxis->SetTitle("eta AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c06->SetLogy();
  c06->cd();
  hist6->SetLineWidth(2);
  hist6->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist6,"2017 signal MC ","f");
  legend->Draw();
  c06->Print("j_eta.png");

  TCanvas *c07 = new TCanvas("c07","mass AK8Jet",1200,900);
  xaxis = hist7->GetXaxis();
  yaxis = hist7->GetYaxis();
  xaxis->SetTitle("mass AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c07->SetLogy();
  c07->cd();
  hist7->SetLineWidth(2);
  hist7->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist7,"2017 signal MC ","f");
  legend->Draw();
  c07->Print("j_mass.png");

  TCanvas *c08 = new TCanvas("c08","mass softdrop AK8Jet",1200,900);
  xaxis = hist8->GetXaxis();
  yaxis = hist8->GetYaxis();
  xaxis->SetTitle("mass softdrop AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0.1,100000);
  yaxis->SetRangeUser(0,400);
  //c08->SetLogy();
  c08->cd();
  hist8->SetLineWidth(2);
  hist8->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist8,"2018 signal MC ","f");
  legend->Draw();
  c08->Print("j_masssd.png");

  TCanvas *c09 = new TCanvas("c09","tau1 AK8Jet",1200,900);
  xaxis = hist9->GetXaxis();
  yaxis = hist9->GetYaxis();
  xaxis->SetTitle("tau1 AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c09->SetLogy();
  c09->cd();
  hist9->SetLineWidth(2);
  hist9->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist9,"2018 signal MC ","f");
  legend->Draw();
  c09->Print("j_tau1.png");

  TCanvas *c10 = new TCanvas("c10","tau2 AK8Jet",1200,900);
  xaxis = hist10->GetXaxis();
  yaxis = hist10->GetYaxis();
  xaxis->SetTitle("tau2 AK8Jet (GeV)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.1,100000);
  c10->SetLogy();
  c10->cd();
  hist10->SetLineWidth(2);
  hist10->Draw("HIST");
  legend->Clear();
  legend->AddEntry(hist9,"2018 signal MC ","f");
  legend->Draw();
  c10->Print("j_tau2.png");
}
