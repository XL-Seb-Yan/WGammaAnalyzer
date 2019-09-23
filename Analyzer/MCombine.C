void MCombine()
{

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  // Plots
  TH1 *hist1 = new TH1F("1","pt_{#gamma}",50,0,3000);
  TH1 *hist2 = new TH1F("2","eta_{#gamma}",50,-2,2);
  TH1 *hist3 = new TH1F("3","pt_{j}",50,0,3000);
  TH1 *hist4 = new TH1F("4","eta_{j}",50,-2,2);
  TH1 *hist5 = new TH1F("5","E_{j}",50,0,3000);
  TH1 *hist6 = new TH1F("6","masssoftdrop_{j}",40,40,120);
  TH1 *hist7 = new TH1F("7","tau21_{j}",50,0,1);
  TH1 *hist8 = new TH1F("8","cos(#theta*)_{p}",50,0,1);
  TH1 *hist9 = new TH1F("9","pt/M",50,0,2);
  TH1 *hist10 = new TH1F("10","invariant mass",50,0,4000);
  TH1 *hist11 = new TH1F("11","seperation",50,0,8);
  
  // Local variables to store to outfile
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21;
  // System
  float sys_costhetastar, sys_ptoverm, sys_invmass;

  // Open input file
  Float_t p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass;
  std::vector<TFile*> file_v;
  // Data
  //file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/SinglePhoton2017_full_presel.root"));
  // GJets
  /*
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/GJets100To200_full_presel.root"));
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/GJets200To400_full_presel.root"));
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/GJets400To600_full_presel.root"));
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/GJets600ToInf_full_presel.root"));
  */
  //QCD
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/QCD300To500_full_presel.root"));
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/QCD500To700_full_presel.root"));
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/QCD700To1000_full_presel.root"));
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/QCD1000To1500_full_presel.root"));
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/QCD1500To2000_full_presel.root"));
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/QCD2000ToInf_full_presel.root"));
  std::vector<double> scale_v;
  //scale_v.push_back(1);
  /*
  scale_v.push_back((double)8300/(double)9959190*41.54*1000);
  scale_v.push_back((double)2200/(double)18536907*41.54*1000);
  scale_v.push_back((double)260/(double)5088564*41.54*1000);
  scale_v.push_back((double)84/(double)3289629*41.54*1000);
  */
  scale_v.push_back((double)311900/(double)60316577*41.54*1000);
  scale_v.push_back((double)29070/(double)56207744*41.54*1000);
  scale_v.push_back((double)5962/(double)47724800*41.54*1000);
  scale_v.push_back((double)1005/(double)16595628*41.54*1000);
  scale_v.push_back((double)101.8/(double)11634434*41.54*1000);
  scale_v.push_back((double)20.54/(double)5941306*41.54*1000);
  
  for(int i = 0; i<file_v.size(); i++){
    TTree* theTree = (TTree*)file_v.at(i)->Get("Events");
    // Improt variables for cutting
    theTree->SetBranchAddress("photon_pt", &p_pt);
    theTree->SetBranchAddress("photon_eta", &p_eta);
    theTree->SetBranchAddress("photon_phi", &p_phi);
    theTree->SetBranchAddress("photon_e", &p_e);
    theTree->SetBranchAddress("ak8puppijet_pt", &j_pt);
    theTree->SetBranchAddress("ak8puppijet_eta", &j_eta);
    theTree->SetBranchAddress("ak8puppijet_phi", &j_phi);
    theTree->SetBranchAddress("ak8puppijet_e", &j_e);
    theTree->SetBranchAddress("ak8puppijet_masssoftdropcorr", &j_mass);
    theTree->SetBranchAddress("ak8puppijet_tau21", &j_tau21);
    theTree->SetBranchAddress("sys_costhetastar", &s_cos);
    theTree->SetBranchAddress("sys_ptoverm", &s_ptm);
    theTree->SetBranchAddress("sys_invmass", &s_mass);
  
    for (int ievt = 0; ievt<theTree->GetEntries();ievt++){
      theTree->GetEntry(ievt);

      double kfac = 1;
      // !!!! k factor only applied to GJets !!!!!
      // Comment following if processing Data and QCD
      /*
      if(p_pt >= 100 && p_pt < 400)
	kfac = 1.8;
      else if(p_pt >= 400 && p_pt < 600)
	kfac = 1.4;
      */
      
      hist1->Fill(p_pt, scale_v.at(i)*kfac);
      hist2->Fill(p_eta, scale_v.at(i)*kfac);
      hist3->Fill(j_pt, scale_v.at(i)*kfac);
      hist4->Fill(j_eta, scale_v.at(i)*kfac);
      hist5->Fill(j_e, scale_v.at(i)*kfac);
      hist6->Fill(j_mass, scale_v.at(i)*kfac);
      hist7->Fill(j_tau21, scale_v.at(i)*kfac);
      hist8->Fill(s_cos, scale_v.at(i)*kfac);
      hist9->Fill(s_ptm, scale_v.at(i)*kfac);
      hist10->Fill(s_mass, scale_v.at(i)*kfac);
    }
  }

  // Write histos to root file
  TFile *outFile = TFile::Open("Histogram_QCD.root","RECREATE");
  hist1->Write();
  hist2->Write();
  hist3->Write();
  hist4->Write();
  hist5->Write();
  hist6->Write();
  hist7->Write();
  hist8->Write();
  hist9->Write();
  hist10->Write();
  outFile->Close();
  
  int color = 2;
  hist1->SetLineColor(color);
  hist2->SetLineColor(color);
  hist3->SetLineColor(color);
  hist4->SetLineColor(color);
  hist5->SetLineColor(color);
  hist6->SetLineColor(color);
  hist7->SetLineColor(color);
  hist8->SetLineColor(color);
  hist9->SetLineColor(color);
  hist10->SetLineColor(color);

  //Non stacked plots

  TLegend *legend = new TLegend(0.65,0.8,0.9,0.9);

  TCanvas *c01 = new TCanvas("c01","pt_{#gamma}",1200,900);
  TAxis *xaxis = hist1->GetXaxis();
  TAxis *yaxis = hist1->GetYaxis();
  xaxis->SetTitle("pt_{#gamma} / M(j#gamma)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  c01->SetLogy();
  c01->cd();
  hist1->SetLineWidth(2);
  hist1->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist1,"2017 SinglePhoton C","f");
  //legend->Draw();
  c01->Print("p_pt.png");

  TCanvas *c02 = new TCanvas("c02","eta_{#gamma}",1200,900);
  xaxis = hist2->GetXaxis();
  yaxis = hist2->GetYaxis();
  xaxis->SetTitle("eta_{#gamma}");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  c02->SetLogy();
  c02->cd();
  hist2->SetLineWidth(2);
  hist2->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist2,"2017 signal MC , pass EleVeto","f");
  //legend->Draw();
  c02->Print("p_eta.png");

  TCanvas *c03 = new TCanvas("c03","pt AK8Jet",1200,900);
  xaxis = hist3->GetXaxis();
  yaxis = hist3->GetYaxis();
  xaxis->SetTitle("pt AK8Jet / M(j#gamma)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  //yaxis->SetRangeUser(0,2000);
  c03->SetLogy();
  c03->cd();
  hist3->SetLineWidth(2);
  hist3->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist5,"2017 signal MC ","f");
  //legend->Draw();
  c03->Print("j_pt.png");

  TCanvas *c04 = new TCanvas("c04","eta AK8Jet",1200,900);
  xaxis = hist4->GetXaxis();
  yaxis = hist4->GetYaxis();
  xaxis->SetTitle("eta AK8Jet");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  c04->SetLogy();
  c04->cd();
  hist4->SetLineWidth(2);
  hist4->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist4,"2017 signal MC ","f");
  //legend->Draw();
  c04->Print("j_eta.png");

  TCanvas *c05 = new TCanvas("c05","E AK8Jet",1200,900);
  xaxis = hist5->GetXaxis();
  yaxis = hist5->GetYaxis();
  xaxis->SetTitle("E AK8Jet / M(j#gamma)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  c05->SetLogy();
  c05->cd();
  hist5->SetLineWidth(2);
  hist5->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist5,"2017 signal MC ","f");
  //legend->Draw();
  c05->Print("j_e.png");

  TCanvas *c06 = new TCanvas("c06","mass softdrop AK8Jet",1200,900);
  xaxis = hist6->GetXaxis();
  yaxis = hist6->GetYaxis();
  xaxis->SetTitle("mass softdrop AK8Jet (GeV)");
  yaxis->SetTitle("a.u.");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  //yaxis->SetRangeUser(0,400);
  c06->SetLogy();
  c06->cd();
  hist6->SetLineWidth(2);
  hist6->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist6,"2017 signal MC ","f");
  //legend->Draw();
  c06->Print("j_masssd.png");

  TCanvas *c07 = new TCanvas("c07","tau21 AK8Jet",1200,900);
  xaxis = hist7->GetXaxis();
  yaxis = hist7->GetYaxis();
  xaxis->SetTitle("tau21 AK8Jet");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  c07->SetLogy();
  c07->cd();
  hist7->SetLineWidth(2);
  hist7->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist7,"2017 signal MC ","f");
  //legend->Draw();
  c07->Print("j_tau21.png");

  TCanvas *c08 = new TCanvas("c8","cos AK8Jet",1200,900);
  xaxis = hist8->GetXaxis();
  yaxis = hist8->GetYaxis();
  xaxis->SetTitle("cos(#theta*)");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  c08->SetLogy();
  c08->cd();
  hist8->SetLineWidth(2);
  hist8->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist8,"2017 signal MC ","f");
  //legend->Draw();
  c08->Print("s_cos.png");

  TCanvas *c09 = new TCanvas("c9","ptm AK8Jet",1200,900);
  xaxis = hist9->GetXaxis();
  yaxis = hist9->GetYaxis();
  xaxis->SetTitle("pt/M");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  c09->SetLogy();
  c09->cd();
  hist9->SetLineWidth(2);
  hist9->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist9,"2017 signal MC ","f");
  //legend->Draw();
  c09->Print("s_ptm.png");

  TCanvas *c10 = new TCanvas("c10","invmass",1200,900);
  xaxis = hist10->GetXaxis();
  yaxis = hist10->GetYaxis();
  xaxis->SetTitle("pt/M");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0.001,10000);
  c10->SetLogy();
  c10->cd();
  hist10->SetLineWidth(2);
  hist10->Draw("HIST");
  legend->Clear();
  //legend->AddEntry(hist10,"2017 signal MC ","f");
  //legend->Draw();
  c10->Print("s_invmass.png");
}
