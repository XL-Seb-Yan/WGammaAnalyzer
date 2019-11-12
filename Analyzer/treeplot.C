void treeplot()
{

  gROOT->SetBatch(1);
  //gStyle->SetOptStat(0);
  // Plots

  TH1F *h = new TH1F("","",100,0,4000);

  // Open input file
  Float_t p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight; 
  //TFile *input = TFile::Open("BackgroundCombinedMC_WGamma_full_full_weightedTo41p54_fitData.root");
  TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017_WGamma_full_full.root");
  TTree* theTree = (TTree*)input->Get("Events");
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
  theTree->SetBranchAddress("xsec_weight", &x_weight);

  double sumW = 0;
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++){
    theTree->GetEntry(ievt);
    
    if(j_mass < 40 || j_mass > 65) continue;
    //if(abs(p_eta) > 1.44) continue;
    //if(abs(j_eta) > 2) continue;
    //if(j_tau21 > 0.3) continue;
    //if(s_ptm < 0.37) continue;
    //if(s_cos > 0.6) continue;
    
    //h->Fill(s_mass,x_weight);
    //sumW += x_weight;
    h->Fill(s_mass);
  }
  cout<<sumW<<endl;

  TCanvas *c01 = new TCanvas("c01","",1200,900);
  TAxis *xaxis = h->GetXaxis();
  TAxis *yaxis = h->GetYaxis();
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.1,10000);
  c01->SetLogy();
  c01->cd();
  h->Draw("HIST");
  c01->Print("test.png");
}
