void TMVApplication( )
{

  //gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  // Plots

  TH1F *hist1 = new TH1F("1","BDT response",50,-0.5,0.5);
  TH1F *hist2 = new TH1F("2","BDT response",50,-1,1);
  TH1F *hist3 = new TH1F("3","BDT response",50,-1,1);
  
  // Local variables to store to outfile
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21;
  // System
  float sys_costhetastar, sys_ptoverm, sys_invmass, sys_seperation, sys_BDT;

  // Create output file
  TFile *outFile = TFile::Open("SinglePhoton2017_WGamma_BDT1000_SwindowWsideband.root", "RECREATE");
  //TFile *outFile = TFile::Open("Signal1000_WGamma_BDT1000_SwindowWwindow.root", "RECREATE");
  TTree *outTree = new TTree("Events","Events"); 
  outTree->Branch("photon_pt",       &photon_pt,      "photon_pt/F");
  outTree->Branch("photon_eta",      &photon_eta,      "photon_eta/F");
  outTree->Branch("photon_phi",      &photon_phi,      "photon_phi/F");
  outTree->Branch("photon_e",        &photon_e,      "photon_e/F");
  outTree->Branch("ak8puppijet_pt",       &ak8puppijet_pt,      "ak8puppijet_pt/F");
  outTree->Branch("ak8puppijet_eta",      &ak8puppijet_eta,      "ak8puppijet_eta/F");
  outTree->Branch("ak8puppijet_phi",      &ak8puppijet_phi,      "ak8puppijet_phi/F");
  outTree->Branch("ak8puppijet_e",        &ak8puppijet_e,      "ak8puppijet_e/F");
  outTree->Branch("ak8puppijet_masssoftdropcorr",   &ak8puppijet_masssoftdropcorr,  "ak8puppijet_masssoftdropcorr/F");
  outTree->Branch("ak8puppijet_tau21",              &ak8puppijet_tau21,             "ak8puppijet_tau21/F");
  outTree->Branch("sys_costhetastar",        &sys_costhetastar,      "sys_costhetastar/F");
  outTree->Branch("sys_ptoverm",             &sys_ptoverm,           "sys_ptoverm/F");
  outTree->Branch("sys_invmass",             &sys_invmass,           "sys_invmass/F");
  outTree->Branch("sys_seperation",          &sys_seperation,        "sys_seperation/F");
  outTree->Branch("sys_BDT",                 &sys_BDT,        "sys_BDT/F");
  
  TMVA::Reader *reader = new TMVA::Reader("Color");
  // Link local values and weight file values
  Float_t p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, s_seperation; 
  reader->AddVariable( "sys_costhetastar", &s_cos );
  reader->AddVariable( "sys_ptoverm", &s_ptm );
  reader->AddVariable( "sys_seperation", &s_seperation );
  reader->AddVariable( "photon_pt", &p_pt );
  reader->AddVariable( "photon_eta", &p_eta );
  reader->AddVariable( "ak8puppijet_pt", &j_pt );
  reader->AddVariable( "ak8puppijet_eta", &j_eta );
  reader->AddVariable( "ak8puppijet_e", &j_e );

  
  reader->BookMVA( "BDT classifier", "for_81_study/MC800/weights/MVAnalysis_BDT.weights.xml" );
  TFile *input = TFile::Open("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017_WGamma_50105_full.root");
  //TFile *input = TFile::Open("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC1000_WGamma_50105_full.root");
  TTree* theTree = (TTree*)input->Get("Events");
  // Improt variables for BDT evaluation
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
  theTree->SetBranchAddress("sys_seperation", &s_seperation);
  theTree->SetBranchAddress("sys_invmass", &s_mass);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);

    if(s_mass > 750 && s_mass < 1250 && j_mass > 50 && j_mass < 65) {
    //if(s_mass > 600 && s_mass < 1000) {
      // BDT evaluation
      Float_t response = reader->EvaluateMVA( "BDT classifier" );

      photon_pt = p_pt;
      photon_eta = p_eta;
      photon_phi = p_phi;
      photon_e = p_e;
      ak8puppijet_pt = j_pt;
      ak8puppijet_eta = j_eta;
      ak8puppijet_phi = j_phi;
      ak8puppijet_e = j_e;
      ak8puppijet_masssoftdropcorr = j_mass;
      ak8puppijet_tau21 = j_tau21;
      sys_costhetastar = s_cos;
      sys_ptoverm = s_ptm;
      sys_seperation = s_seperation;
      sys_invmass = s_mass;
      sys_BDT = response;
      hist1->Fill(response);
  
      outTree->Fill();
    }
  }
  outFile->Write();
  outFile->Close();
  delete reader;

  hist1->SetLineColor(2);
  hist1->SetLineWidth(3);
  hist1->Scale(1/double(hist1->GetEntries()));

  TLegend *legend = new TLegend(0.7,0.75,0.85,0.85);
  TCanvas *c01 = new TCanvas("c01","",1200,900);
  TAxis *xaxis = hist1->GetXaxis();
  TAxis *yaxis = hist1->GetYaxis();
  xaxis->SetTitle("BDT");
  yaxis->SetTitle("Entries / 0.1");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.0001,1);
  c01->SetLogy();
  c01->cd();
  hist1->Draw("HIST");
  c01->Print("BDT.png");
  
  /*
  hist1->SetLineColor(2);
  hist2->SetLineColor(4);
  hist1->SetLineWidth(3);
  hist2->SetLineWidth(3);
  hist1->Scale(1/double(hist1->GetEntries()));
  hist2->Scale(1/double(hist2->GetEntries()));

  TLegend *legend = new TLegend(0.7,0.75,0.85,0.85);
  TCanvas *c01 = new TCanvas("c01","",1200,900);
  TAxis *xaxis = hist2->GetXaxis();
  TAxis *yaxis = hist2->GetYaxis();
  xaxis->SetTitle("pt_{#gamma}/M");
  yaxis->SetTitle("Entries / 0.1");
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.2);
  //yaxis->SetRangeUser(0.0001,1);
  //c01->SetLogy();
  c01->cd();
  hist2->Draw("HIST");
  hist1->Draw("SAMEHIST");
  legend->Clear();
  legend->AddEntry(hist1,"65-70 GeV" ,"f");
  legend->AddEntry(hist2,"70-100 GeV","f");
  legend->Draw();
  c01->Print("BDT.png");
  */
}
