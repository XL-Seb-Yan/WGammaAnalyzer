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

Int_t count1 = 0;
Int_t count2 = 0;
TH1 *hist1 = new TH1F("1","pt_{#gamma}",100,0,1000);
TH1 *hist2 = new TH1F("2","pt_{#gamma}",100,0,1000);
TH1 *hist3 = new TH1F("3","pt_{#gamma}",100,0,1000);

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
   
   if (*ph_N > 0){

     bool passTrig_175 = false;
     bool passTrig_165 = false;
     for(int i=0; i<HLT_isFired.GetSize(); i++){
       if (HLT_isFired[i].first.find("HLT_Photon175_") != std::string::npos)
	 passTrig_175 = HLT_isFired[i].second;
       //if (HLT_isFired[i].first.find("HLT_Photon165_HE10_") != std::string::npos)
       //passTrig_165 = HLT_isFired[i].second;
     }
     bool passTrig = (passTrig_175 || passTrig_165);
     hist1->Fill(ph_pt[0]);
     if(ph_passLooseId[0] == 1){
       hist3->Fill(ph_pt[0]);
       if(passTrig){
	 hist2->Fill(ph_pt[0]);
       }
     }
       /*
	 else{
	 if (entry < 10000){
	 for(int i=0; i<HLT_isFired.GetSize(); i++){
	 cout<<HLT_isFired[i].second<<" "<<HLT_isFired[i].first<<endl;
	 }
	 cout<<endl;
	 }
	 }
       */
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
  yaxis->SetTitle("Entries / 50 GeV");
  yaxis->SetTitleOffset(1.3);
  //yaxis->SetRangeUser(0,0.14);
  c01->SetLogy();
  c01->cd();
  hist1->SetLineWidth(2);
  hist1->Draw("HIST");
  hist3->SetLineWidth(2);
  hist3->SetLineColor(2);
  hist3->Draw("SAMEHIST");
  hist2->SetFillColor(30);
  hist2->SetLineColor(30);
  hist2->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hist1,"2016 SingleMuon C, with photons","f");
  legend->AddEntry(hist3,"2016 SingleMuon C, pass loose ID","f");
  legend->AddEntry(hist2,"2016 SingleMuon C, fired trigger","f");
  legend->Draw();
  c01->Print("p_pt.png");

  TCanvas *c02 = new TCanvas("c02","Trigger Efficiency",1200,900);
  c02->cd();
  Eff->Draw();
  c02->Print("eff.png");
}
