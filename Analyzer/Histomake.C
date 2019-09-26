//Combined subrange MCs to a single MC file
void Histomake()
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

  // Open input file
  Float_t p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight;
  std::vector<TFile*> file_v;
  // Data
  file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/SinglePhoton2017_full_presel.root"));
  // GJets
  //file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/GJetsCombined_full_presel_weightedTo41p54_fitData.root"));
  // QCD
  //file_v.push_back(TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Analyzer/QCDCombined_full_presel_weightedTo41p54_fitData.root"));
  
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
    //theTree->SetBranchAddress("xsec_weight", &x_weight);
  
    for (int ievt = 0; ievt<theTree->GetEntries();ievt++){
      theTree->GetEntry(ievt);

      double weight = 1;
      hist1->Fill(p_pt, weight);
      hist2->Fill(p_eta, weight);
      hist3->Fill(j_pt, weight);
      hist4->Fill(j_eta, weight);
      hist5->Fill(j_e, weight);
      hist6->Fill(j_mass, weight);
      hist7->Fill(j_tau21, weight);
      hist8->Fill(s_cos, weight);
      hist9->Fill(s_ptm, weight);
      hist10->Fill(s_mass, weight);
    }
  }

  // Write histos to root file
  TFile *outFile = TFile::Open("Histogram_Data_fullyweighted.root","RECREATE");
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
}
