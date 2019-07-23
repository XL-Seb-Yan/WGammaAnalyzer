void TMVApplication( )
{

  // Plots
  // Local variables to store to outfile
  // Photon
  float photon_pt, photon_eta, photon_phi, photon_e, photon_mvaval, photon_mvacat;
  // Jet
  float ak8puppijet_pt, ak8puppijet_eta, ak8puppijet_phi, ak8puppijet_e, ak8puppijet_masssoftdropcorr, ak8puppijet_tau21;
  // System
  float sys_costhetastar, sys_ptoverm, sys_invmass, sys_seperation, sys_BDT;

  // Create output file
  //TFile *outFile = TFile::Open("SinglePhoton2017BCD_WGamma_BDT800_SwindowWsideband.root", "RECREATE");
  TFile *outFile = TFile::Open("Signal800_WGamma_BDT800_SwindowWwindow.root", "RECREATE");
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

  
  reader->BookMVA( "BDT classifier", "MC800/weights/MVAnalysis_BDT.weights.xml" );
  //TFile *input = TFile::Open("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017BCD_WGamma_50105_full.root");
  TFile *input = TFile::Open("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC800_WGamma_50105_full.root");
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

    if(s_mass > 600 && s_mass < 1000 && j_mass > 65 && j_mass < 95) {
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
  
      outTree->Fill();
    }
  }
  outFile->Write();
  outFile->Close();
  delete reader;
}
