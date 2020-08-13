#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
void kfactorfit()
{

   gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Preliminary";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 4;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);
  
  // Plots
  TH1 *hist1_1 = new TH1F("1_1","pt_{#gamma}",100,0,3000);
  TH1 *hist1_2 = new TH1F("1_2","eta_{#gamma}",50,-2,2);
  TH1 *hist1_3 = new TH1F("1_3","pt_{j}",100,0,3000);
  TH1 *hist1_4 = new TH1F("1_4","eta_{j}",50,-2,2);
  TH1 *hist1_5 = new TH1F("1_5","E_{j}",100,0,3000);
  TH1 *hist1_6 = new TH1F("1_6","masssoftdrop_{j}",50,40,140);
  TH1 *hist1_7 = new TH1F("1_7","tau21_{j}",50,0,1);
  TH1 *hist1_8 = new TH1F("1_8","cos(#theta*)_{p}",50,0,1);
  TH1 *hist1_9 = new TH1F("1_9","pt/M",50,0,2);
  TH1 *hist1_10 = new TH1F("1_10","invariant mass",100,0,4000);
  TH1 *hist1_11 = new TH1F("1_11","PV_N",50,0,100);
  TH1 *hist2_1 = new TH1F("2_1","pt_{#gamma}",100,0,3000);
  TH1 *hist2_2 = new TH1F("2_2","eta_{#gamma}",50,-2,2);
  TH1 *hist2_3 = new TH1F("2_3","pt_{j}",100,0,3000);
  TH1 *hist2_4 = new TH1F("2_4","eta_{j}",50,-2,2);
  TH1 *hist2_5 = new TH1F("2_5","E_{j}",100,0,3000);
  TH1 *hist2_6 = new TH1F("2_6","masssoftdrop_{j}",50,40,140);
  TH1 *hist2_7 = new TH1F("2_7","tau21_{j}",50,0,1);
  TH1 *hist2_8 = new TH1F("2_8","cos(#theta*)_{p}",50,0,1);
  TH1 *hist2_9 = new TH1F("2_9","pt/M",50,0,2);
  TH1 *hist2_10 = new TH1F("2_10","invariant mass",100,0,4000);
  TH1 *hist2_11 = new TH1F("2_11","PV_N",50,0,100);
  TH1 *hist3_1 = new TH1F("3_1","pt_{#gamma}",100,0,3000);
  TH1 *hist3_2 = new TH1F("3_2","eta_{#gamma}",50,-2,2);
  TH1 *hist3_3 = new TH1F("3_3","pt_{j}",100,0,3000);
  TH1 *hist3_4 = new TH1F("3_4","eta_{j}",50,-2,2);
  TH1 *hist3_5 = new TH1F("3_5","E_{j}",100,0,3000);
  TH1 *hist3_6 = new TH1F("3_6","masssoftdrop_{j}",50,40,140);
  TH1 *hist3_7 = new TH1F("3_7","tau21_{j}",50,0,1);
  TH1 *hist3_8 = new TH1F("3_8","cos(#theta*)_{p}",50,0,1);
  TH1 *hist3_9 = new TH1F("3_9","pt/M",50,0,2);
  TH1 *hist3_10 = new TH1F("3_10","invariant mass",100,0,4000);
  TH1 *hist3_11 = new TH1F("3_11","PV_N",50,0,100);

  // Open input file
  float p_pt, p_eta, p_phi, p_e, j_pt, j_eta, j_phi, j_e, j_mass, j_tau21, s_cos, s_ptm, s_mass, x_weight, x_puweight, x_sf;
  int s_PV;
  
  // Data
  TFile* f_data = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/Run2Data_postproc_WGammaRun2_full_full_jmcorr_May22.root");
  TFile* f_gjets = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/GJets_postproc_WGammaRun2_full_full_jmcorr_May22.root");
  TFile* f_qcd = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/DATA/Run2/presel/QCD_postproc_WGammaRun2_full_full_jmcorr_May22.root");
  
  cout<<"Processing data"<<endl;
  
    TTree* theTree = (TTree*)f_data->Get("Events");
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
    theTree->SetBranchAddress("m", &s_mass);
    theTree->SetBranchAddress("xsec_weight", &x_weight);
	theTree->SetBranchAddress("xsec_puweight", &x_puweight);
	theTree->SetBranchAddress("xsec_sf", &x_sf);
	theTree->SetBranchAddress("sys_pvn", &s_PV);
  
  
    for (int ievt = 0; ievt<theTree->GetEntries();ievt++){
      theTree->GetEntry(ievt);
	  // if(p_pt < 225) continue;
	  // if(j_pt < 225) continue;
      hist1_1->Fill(p_pt);
      hist1_2->Fill(p_eta);
      hist1_3->Fill(j_pt);
      hist1_4->Fill(j_eta);
      hist1_5->Fill(j_e);
      hist1_6->Fill(j_mass);
      hist1_7->Fill(j_tau21);
      hist1_8->Fill(s_cos);
      hist1_9->Fill(s_ptm);
      hist1_10->Fill(s_mass);
	  hist1_11->Fill(s_PV);
    }
	
	cout<<"Processing GJets"<<endl;
	theTree = (TTree*)f_gjets->Get("Events");
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
    theTree->SetBranchAddress("m", &s_mass);
    theTree->SetBranchAddress("xsec_weight", &x_weight);
	theTree->SetBranchAddress("xsec_puweight", &x_puweight);
	theTree->SetBranchAddress("xsec_sf", &x_sf);
	theTree->SetBranchAddress("sys_pvn", &s_PV);
  
  
    for (int ievt = 0; ievt<theTree->GetEntries();ievt++){
      theTree->GetEntry(ievt);
	  // if(p_pt < 225) continue;
	  // if(j_pt < 225) continue;
	  x_weight = x_weight * 3.303395;
      hist2_1->Fill(p_pt, x_weight*x_puweight*x_sf);
      hist2_2->Fill(p_eta, x_weight*x_puweight*x_sf);
      hist2_3->Fill(j_pt, x_weight*x_puweight*x_sf);
      hist2_4->Fill(j_eta, x_weight*x_puweight*x_sf);
      hist2_5->Fill(j_e, x_weight*x_puweight*x_sf);
      hist2_6->Fill(j_mass, x_weight*x_puweight*x_sf);
      hist2_7->Fill(j_tau21, x_weight*x_puweight*x_sf);
      hist2_8->Fill(s_cos, x_weight*x_puweight*x_sf);
      hist2_9->Fill(s_ptm, x_weight*x_puweight*x_sf);
      hist2_10->Fill(s_mass, x_weight*x_puweight*x_sf);
	  hist2_11->Fill(s_PV, x_weight*x_puweight*x_sf);
    }
	
	cout<<"Processing QCD"<<endl;
    theTree = (TTree*)f_qcd->Get("Events");
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
    theTree->SetBranchAddress("m", &s_mass);
    theTree->SetBranchAddress("xsec_weight", &x_weight);
	theTree->SetBranchAddress("xsec_puweight", &x_puweight);
	theTree->SetBranchAddress("xsec_sf", &x_sf);
	theTree->SetBranchAddress("sys_pvn", &s_PV);
  
  
    for (int ievt = 0; ievt<theTree->GetEntries();ievt++){
      theTree->GetEntry(ievt);
	  // if(p_pt < 225) continue;
	  // if(j_pt < 225) continue;
	  x_weight = x_weight * 3.303395;
      hist3_1->Fill(p_pt, x_weight*x_puweight*x_sf);
      hist3_2->Fill(p_eta, x_weight*x_puweight*x_sf);
      hist3_3->Fill(j_pt, x_weight*x_puweight*x_sf);
      hist3_4->Fill(j_eta, x_weight*x_puweight*x_sf);
      hist3_5->Fill(j_e, x_weight*x_puweight*x_sf);
      hist3_6->Fill(j_mass, x_weight*x_puweight*x_sf);
      hist3_7->Fill(j_tau21, x_weight*x_puweight*x_sf);
      hist3_8->Fill(s_cos, x_weight*x_puweight*x_sf);
      hist3_9->Fill(s_ptm, x_weight*x_puweight*x_sf);
      hist3_10->Fill(s_mass, x_weight*x_puweight*x_sf);
	  hist3_11->Fill(s_PV, x_weight*x_puweight*x_sf);
    }
	
  // Least Square fit result: GJets: 1.41, QCD: 0.57
  // Least Square fitting to invariant mass histograms
  TGraph2D *g = new TGraph2D();
  int NBins = hist1_1->GetNbinsX();
  cout<<"Number of bins: "<<NBins<<endl;
  int point = 0;
  double minls = 99999999999;
  double sGJets = -99;
  double sQCD = -99;
  for(int i=20; i<150; i++){
    for(int j=20; j<150; j++){
      double ls = 0;
      double scaleGJets = i*0.01;
      double scaleQCD = j*0.01;
      for(int k=7; k<NBins+1; k++){ //starts from 210. GeV
	    double NData = hist1_1->GetBinContent(k);
	    double NGjets = hist2_1->GetBinContent(k);
	    double NQCD = hist3_1->GetBinContent(k);
	    //cout<<NData<<" "<<NGjets<<" "<<NQCD<<" "<<pow((NData-(scaleGJets*NGjets + scaleQCD*NQCD))/NData,2)<<endl;
	    if(NData == 0) continue;//exclude first several empty bins lower than trigger turn on
	    ls += pow((NData-(scaleGJets*NGjets + scaleQCD*NQCD))/1000,2);
      }
      if(ls < minls){
	    minls = ls;
	    sGJets = scaleGJets;
	    sQCD = scaleQCD;
      }
      cout<<"Runnig through "<<scaleGJets<<" "<<scaleQCD<<" Least Square is: "<<ls<<endl;
      g->SetPoint(point,scaleGJets,scaleQCD,ls);
      point++;
    }
  }
   cout<<sGJets<<" "<<sQCD<<" "<<minls<<endl;
   
   sGJets = 1.35;
   sQCD = 1.2;
   
  TCanvas *c = new TCanvas("c","scale",1200,900);
  TAxis *xaxis = g->GetXaxis();
  TAxis *yaxis = g->GetYaxis();
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  xaxis->SetTitle("GJets scale");
  yaxis->SetTitle("QCD scale");
  //yaxis->SetRangeUser(0.5,10000000);
  c->SetLogz();
  c->cd();
  c->SetRightMargin(0.15);
  c->SetGrid();
  g->Draw("COLZ");
  cout<<"OK"<<endl;
  c->Print("MCFit.png");
   
  hist2_1->Scale(sGJets);
  hist2_2->Scale(sGJets);
  hist2_3->Scale(sGJets);
  hist2_4->Scale(sGJets);
  hist2_5->Scale(sGJets);
  hist2_6->Scale(sGJets);
  hist2_7->Scale(sGJets);
  hist2_8->Scale(sGJets);
  hist2_9->Scale(sGJets);
  hist2_10->Scale(sGJets);
  hist2_11->Scale(sGJets);
  hist3_1->Scale(sQCD);
  hist3_2->Scale(sQCD);
  hist3_3->Scale(sQCD);
  hist3_4->Scale(sQCD);
  hist3_5->Scale(sQCD);
  hist3_6->Scale(sQCD);
  hist3_7->Scale(sQCD);
  hist3_8->Scale(sQCD);
  hist3_9->Scale(sQCD);
  hist3_10->Scale(sQCD);
  hist3_11->Scale(sQCD);
  
  int color1 = 2;
  int color2 = 8;
  int color3 = 7;

  hist1_1->SetLineColor(1);
  hist1_2->SetLineColor(1);
  hist1_3->SetLineColor(1);
  hist1_4->SetLineColor(1);
  hist1_5->SetLineColor(1);
  hist1_6->SetLineColor(1);
  hist1_7->SetLineColor(1);
  hist1_8->SetLineColor(1);
  hist1_9->SetLineColor(1);
  hist1_10->SetLineColor(1);
  hist1_11->SetLineColor(1);
  hist2_1->SetLineColor(1);
  hist2_2->SetLineColor(1);
  hist2_3->SetLineColor(1);
  hist2_4->SetLineColor(1);
  hist2_5->SetLineColor(1);
  hist2_6->SetLineColor(1);
  hist2_7->SetLineColor(1);
  hist2_8->SetLineColor(1);
  hist2_9->SetLineColor(1);
  hist2_10->SetLineColor(1);
  hist2_11->SetLineColor(1);
  hist3_1->SetLineColor(1);
  hist3_2->SetLineColor(1);
  hist3_3->SetLineColor(1);
  hist3_4->SetLineColor(1);
  hist3_5->SetLineColor(1);
  hist3_6->SetLineColor(1);
  hist3_7->SetLineColor(1);
  hist3_8->SetLineColor(1);
  hist3_9->SetLineColor(1);
  hist3_10->SetLineColor(1);
  hist3_11->SetLineColor(1);

  hist2_1->SetFillColor(color2);
  hist2_2->SetFillColor(color2);
  hist2_3->SetFillColor(color2);
  hist2_4->SetFillColor(color2);
  hist2_5->SetFillColor(color2);
  hist2_6->SetFillColor(color2);
  hist2_7->SetFillColor(color2);
  hist2_8->SetFillColor(color2);
  hist2_9->SetFillColor(color2);
  hist2_10->SetFillColor(color2);
  hist2_11->SetFillColor(color2);
  hist3_1->SetFillColor(color3);
  hist3_2->SetFillColor(color3);
  hist3_3->SetFillColor(color3);
  hist3_4->SetFillColor(color3);
  hist3_5->SetFillColor(color3);
  hist3_6->SetFillColor(color3);
  hist3_7->SetFillColor(color3);
  hist3_8->SetFillColor(color3);
  hist3_9->SetFillColor(color3);
  hist3_10->SetFillColor(color3);
  hist3_11->SetFillColor(color3);
  
  hist1_1->SetMarkerStyle(8);
  hist1_2->SetMarkerStyle(8);
  hist1_3->SetMarkerStyle(8);
  hist1_4->SetMarkerStyle(8);
  hist1_5->SetMarkerStyle(8);
  hist1_6->SetMarkerStyle(8);
  hist1_7->SetMarkerStyle(8);
  hist1_8->SetMarkerStyle(8);
  hist1_9->SetMarkerStyle(8);
  hist1_10->SetMarkerStyle(8);
  hist1_11->SetMarkerStyle(8);
  hist1_1->SetMarkerSize(1.5);
  hist1_2->SetMarkerSize(1.5);
  hist1_3->SetMarkerSize(1.5);
  hist1_4->SetMarkerSize(1.5);
  hist1_5->SetMarkerSize(1.5);
  hist1_6->SetMarkerSize(1.5);
  hist1_7->SetMarkerSize(1.5);
  hist1_8->SetMarkerSize(1.5);
  hist1_9->SetMarkerSize(1.5);
  hist1_10->SetMarkerSize(1.5);
  hist1_11->SetMarkerSize(1.5);
  
  //Stacked plot
  THStack *stack1 = new THStack("stack1","pt_{p}");
  THStack *stack2 = new THStack("stack2","eta_{p}");
  THStack *stack3 = new THStack("stack3","pt_{j}");
  THStack *stack4 = new THStack("stack4","eta_{j}");
  THStack *stack5 = new THStack("stack5","e_{j}");
  THStack *stack6 = new THStack("stack6","Jet softdrop mass");
  THStack *stack7 = new THStack("stack7","Jet tau21");
  THStack *stack8 = new THStack("stack8","cos(#theta*)");
  THStack *stack9 = new THStack("stack9","pt_{j} / M");
  THStack *stack10 = new THStack("stack10","Invariant mass");
  THStack *stack11 = new THStack("stack11","PV");
  stack1->Add(hist3_1); stack1->Add(hist2_1);
  stack2->Add(hist3_2); stack2->Add(hist2_2);
  stack3->Add(hist3_3); stack3->Add(hist2_3);
  stack4->Add(hist3_4); stack4->Add(hist2_4);
  stack5->Add(hist3_5); stack5->Add(hist2_5);
  stack6->Add(hist3_6); stack6->Add(hist2_6);
  stack7->Add(hist3_7); stack7->Add(hist2_7);
  stack8->Add(hist3_8); stack8->Add(hist2_8);
  stack9->Add(hist3_9); stack9->Add(hist2_9);
  stack10->Add(hist3_10); stack10->Add(hist2_10);
  stack11->Add(hist3_11); stack11->Add(hist2_11);

  TLegend *legend = new TLegend(0.54,0.78,0.9,0.9);
  // Residual plot
  TH1 *data = NULL;
  TH1 *bkg = NULL;
  TH1 *pull = NULL;
  TAxis *xaxis1 = NULL;
  TAxis *yaxis1 = NULL;
  TAxis *xaxis2 = NULL;
  TAxis *yaxis2 = NULL;

  //===========================================================
  //gStyle->SetHistMinimumZero();
  TCanvas *c01 = new TCanvas("c01","",2400,2200);
  c01->cd();
  TPad *p01a = new TPad("p01a","p01a",0.1,0.30,0.9,1.0);
  TPad *p01b = new TPad("p01b","p01b",0.1,0.1,0.9,0.315);
  p01a->Draw();
  p01b->Draw();
  p01a->cd();
  p01a->SetBottomMargin(0.11);
  p01a->SetLogy();
  xaxis1 = hist1_1->GetXaxis();
  yaxis1 = hist1_1->GetYaxis();
  xaxis1->SetTitle("pt_{#gamma} (GeV)");
  yaxis1->SetTitle("Entries / 60 GeV");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_1->Draw("E1");
  stack1->Draw("SAMEHIST");
  hist1_1->Draw("E1SAME");
  hist1_1->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_1,"Run2 Data","lep");
  legend->AddEntry(hist2_1,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_1,"2017 MC, QCD(weighted)","f");
  legend->Draw();
  
  p01b->cd();
  p01b->SetTopMargin(0.037);
  p01b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("pt_{#gamma} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p01b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c01->Print("p_pt.png");
  c01->Print("p_pt.pdf");
  c01->Print("p_pt.svg");
  c01->Print("p_pt.root");
  //==========================================================

  //===========================================================
  TCanvas *c02 = new TCanvas("c02","",2400,2200);
  c02->cd();
  TPad *p02a = new TPad("p02a","p02a",0.1,0.30,0.9,1.0);
  TPad *p02b = new TPad("p02b","p02b",0.1,0.1,0.9,0.315);
  p02a->Draw();
  p02b->Draw();
  p02a->cd();
  p02a->SetBottomMargin(0.11);
  p02a->SetLogy();
  xaxis1 = hist1_2->GetXaxis();
  yaxis1 = hist1_2->GetYaxis();
  xaxis1->SetTitle("#eta_{#gamma}");
  yaxis1->SetTitle("Entries / 0.08");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_2->Draw("E1");
  stack2->Draw("SAMEHIST");
  hist1_2->Draw("E1SAME");
  hist1_2->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_2,"Run2 Data","lep");
  legend->AddEntry(hist2_2,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_2,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p02b->cd();
  p02b->SetTopMargin(0.037);
  p02b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_2->Clone();
  bkg = (TH1*)hist2_2->Clone();
  bkg->Add((TH1*)hist3_2->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("#eta_{#gamma}");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p02b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c02->Print("p_eta.png");
  c02->Print("p_eta.pdf");
  c02->Print("p_eta.svg");
  c02->Print("p_eta.root");
  //==========================================================

  //===========================================================
  TCanvas *c03 = new TCanvas("c03","",2400,2200);
  c03->cd();
  TPad *p03a = new TPad("p03a","p03a",0.1,0.30,0.9,1.0);
  TPad *p03b = new TPad("p03b","p03b",0.1,0.1,0.9,0.315);
  p03a->Draw();
  p03b->Draw();
  p03a->cd();
  p03a->SetBottomMargin(0.11);
  p03a->SetLogy();
  xaxis1 = hist1_3->GetXaxis();
  yaxis1 = hist1_3->GetYaxis();
  xaxis1->SetTitle("pt_{AK8jet} (GeV)");
  yaxis1->SetTitle("Entries / 60 (GeV)");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_3->Draw("E1");
  stack3->Draw("SAMEHIST");
  hist1_3->Draw("E1SAME");
  hist1_3->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_3,"Run2 Data","lep");
  legend->AddEntry(hist2_3,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_3,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p03b->cd();
  p03b->SetTopMargin(0.037);
  p03b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_3->Clone();
  bkg = (TH1*)hist2_3->Clone();
  bkg->Add((TH1*)hist3_3->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("pt_{AK8jet} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p03b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c03->Print("j_pt.png");
  c03->Print("j_pt.pdf");
  c03->Print("j_pt.svg");
  c03->Print("j_pt.root");
  //==========================================================

  //===========================================================
  TCanvas *c04 = new TCanvas("c04","",2400,2200);
  c04->cd();
  TPad *p04a = new TPad("p04a","p04a",0.1,0.30,0.9,1.0);
  TPad *p04b = new TPad("p04b","p04b",0.1,0.1,0.9,0.315);
  p04a->Draw();
  p04b->Draw();
  p04a->cd();
  p04a->SetBottomMargin(0.11);
  p04a->SetLogy();
  xaxis1 = hist1_4->GetXaxis();
  yaxis1 = hist1_4->GetYaxis();
  xaxis1->SetTitle("#eta_{AK8jet}");
  yaxis1->SetTitle("Entries / 0.08");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_4->Draw("E1");
  stack4->Draw("SAMEHIST");
  hist1_4->Draw("E1SAME");
  hist1_4->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_4,"Run2 Data","lep");
  legend->AddEntry(hist2_4,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_4,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p04b->cd();
  p04b->SetTopMargin(0.037);
  p04b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_4->Clone();
  bkg = (TH1*)hist2_4->Clone();
  bkg->Add((TH1*)hist3_4->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("#eta_{AK8jet}");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p04b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c04->Print("j_eta.png");
  c04->Print("j_eta.pdf");
  c04->Print("j_eta.svg");
  c04->Print("j_eta.root");
  //==========================================================

  //===========================================================
  TCanvas *c05 = new TCanvas("c05","",2400,2200);
  c05->cd();
  TPad *p05a = new TPad("p05a","p05a",0.1,0.30,0.9,1.0);
  TPad *p05b = new TPad("p05b","p05b",0.1,0.1,0.9,0.315);
  p05a->Draw();
  p05b->Draw();
  p05a->cd();
  p05a->SetBottomMargin(0.11);
  p05a->SetLogy();
  xaxis1 = hist1_5->GetXaxis();
  yaxis1 = hist1_5->GetYaxis();
  xaxis1->SetTitle("E_{AK8jet} (GeV)");
  yaxis1->SetTitle("Entries / 60 (GeV)");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_5->Draw("E1");
  stack5->Draw("SAMEHIST");
  hist1_5->Draw("E1SAME");
  hist1_5->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_5,"Run2 Data","lep");
  legend->AddEntry(hist2_5,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_5,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p05b->cd();
  p05b->SetTopMargin(0.037);
  p05b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_5->Clone();
  bkg = (TH1*)hist2_5->Clone();
  bkg->Add((TH1*)hist3_5->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("E_{AK8jet} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p05b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c05->Print("j_e.png");
  c05->Print("j_e.pdf");
  c05->Print("j_e.svg");
  c05->Print("j_e.root");
  //==========================================================

  //===========================================================
  TCanvas *c06 = new TCanvas("c06","",2400,2200);
  c06->cd();
  TPad *p06a = new TPad("p06a","p06a",0.1,0.30,0.9,1.0);
  TPad *p06b = new TPad("p06b","p06b",0.1,0.1,0.9,0.315);
  p06a->Draw();
  p06b->Draw();
  p06a->cd();
  p06a->SetBottomMargin(0.11);
  //p06a->SetLogy();
  xaxis1 = hist1_6->GetXaxis();
  yaxis1 = hist1_6->GetYaxis();
  xaxis1->SetTitle("Softdrop m_{AK8jet} (GeV)");
  yaxis1->SetTitle("Entries / 2 (GeV)");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0,25000);
  xaxis1->SetLimits(30,150);
  hist1_6->Draw("E1");
  stack6->Draw("SAMEHIST");
  hist1_6->Draw("E1SAME");
  hist1_6->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_6,"Run2 Data","lep");
  legend->AddEntry(hist2_6,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_6,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p06b->cd();
  p06b->SetTopMargin(0.037);
  p06b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_6->Clone();
  bkg = (TH1*)hist2_6->Clone();
  bkg->Add((TH1*)hist3_6->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("Softdrop m_{AK8jet} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  xaxis2->SetLimits(30,120);
  p06b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c06->Print("j_msd.png");
  c06->Print("j_msd.pdf");
  c06->Print("j_msd.svg");
  c06->Print("j_msd.root");
  //==========================================================

  //===========================================================
  TCanvas *c07 = new TCanvas("c07","",2400,2200);
  c07->cd();
  TPad *p07a = new TPad("p07a","p07a",0.1,0.30,0.9,1.0);
  TPad *p07b = new TPad("p07b","p07b",0.1,0.1,0.9,0.315);
  p07a->Draw();
  p07b->Draw();
  p07a->cd();
  p07a->SetBottomMargin(0.11);
  p07a->SetLogy();
  xaxis1 = hist1_7->GetXaxis();
  yaxis1 = hist1_7->GetYaxis();
  xaxis1->SetTitle("AK8 jet #tau_{21}");
  yaxis1->SetTitle("Entries / 0.02");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_7->Draw("E1");
  stack7->Draw("SAMEHIST");
  hist1_7->Draw("E1SAME");
  hist1_7->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_7,"Run2 Data","lep");
  legend->AddEntry(hist2_7,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_7,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p07b->cd();
  p07b->SetTopMargin(0.037);
  p07b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_7->Clone();
  bkg = (TH1*)hist2_7->Clone();
  bkg->Add((TH1*)hist3_7->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("AK8 jet #tau_{21}");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p07b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c07->Print("j_tau21.png");
  c07->Print("j_tau21.pdf");
  c07->Print("j_tau21.svg");
  c07->Print("j_tau21.root");
  //==========================================================

   //===========================================================
  TCanvas *c08 = new TCanvas("c08","",2400,2200);
  c08->cd();
  TPad *p08a = new TPad("p08a","p08a",0.1,0.30,0.9,1.0);
  TPad *p08b = new TPad("p08b","p08b",0.1,0.1,0.9,0.315);
  p08a->Draw();
  p08b->Draw();
  p08a->cd();
  p08a->SetBottomMargin(0.11);
  p08a->SetLogy();
  xaxis1 = hist1_8->GetXaxis();
  yaxis1 = hist1_8->GetYaxis();
  xaxis1->SetTitle("cos(#theta*)");
  yaxis1->SetTitle("Entries / 0.02");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_8->Draw("E1");
  stack8->Draw("SAMEHIST");
  hist1_8->Draw("E1SAME");
  hist1_8->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_8,"Run2 Data","lep");
  legend->AddEntry(hist2_8,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_8,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p08b->cd();
  p08b->SetTopMargin(0.037);
  p08b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_8->Clone();
  bkg = (TH1*)hist2_8->Clone();
  bkg->Add((TH1*)hist3_8->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("cos(#theta*)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p08b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c08->Print("s_cos.png");
  c08->Print("s_cos.pdf");
  c08->Print("s_cos.svg");
  c08->Print("s_cos.root");
  //==========================================================

     //===========================================================
  TCanvas *c09 = new TCanvas("c09","",2400,2200);
  c09->cd();
  TPad *p09a = new TPad("p09a","p09a",0.1,0.30,0.9,1.0);
  TPad *p09b = new TPad("p09b","p09b",0.1,0.1,0.9,0.315);
  p09a->Draw();
  p09b->Draw();
  p09a->cd();
  p09a->SetBottomMargin(0.11);
  p09a->SetLogy();
  xaxis1 = hist1_9->GetXaxis();
  yaxis1 = hist1_9->GetYaxis();
  xaxis1->SetTitle("pt_{#gamma}/M_{j#gamma}");
  yaxis1->SetTitle("Entries / 0.04");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_9->Draw("E1");
  stack9->Draw("SAMEHIST");
  hist1_9->Draw("E1SAME");
  hist1_9->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_9,"Run2 Data","lep");
  legend->AddEntry(hist2_9,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_9,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p09b->cd();
  p09b->SetTopMargin(0.037);
  p09b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_9->Clone();
  bkg = (TH1*)hist2_9->Clone();
  bkg->Add((TH1*)hist3_9->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("pt_{#gamma}/M_{j#gamma}");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p09b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c09->Print("s_ptm.png");
  c09->Print("s_ptms.pdf");
  c09->Print("s_ptm.svg");
  c09->Print("s_ptm.root");
  //==========================================================

  //===========================================================
  TCanvas *c10 = new TCanvas("c10","",2400,2200);
  c10->cd();
  TPad *p10a = new TPad("p10a","p10a",0.1,0.30,0.9,1.0);
  TPad *p10b = new TPad("p10b","p10b",0.1,0.1,0.9,0.315);
  p10a->Draw();
  p10b->Draw();
  p10a->cd();
  p10a->SetBottomMargin(0.11);
  p10a->SetLogy();
  xaxis1 = hist1_10->GetXaxis();
  yaxis1 = hist1_10->GetYaxis();
  xaxis1->SetTitle("M_{j#gamma} (GeV)");
  yaxis1->SetTitle("Entries / 80 GeV");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_10->Draw("E1");
  stack10->Draw("SAMEHIST");
  hist1_10->Draw("E1SAME");
  hist1_10->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_10,"Run2 Data","lep");
  legend->AddEntry(hist2_10,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_10,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p10b->cd();
  p10b->SetTopMargin(0.037);
  p10b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_10->Clone();
  bkg = (TH1*)hist2_10->Clone();
  bkg->Add((TH1*)hist3_10->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("M_{j#gamma} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p10b->SetGrid();
  p10b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c10->Print("s_m.png");
  c10->Print("s_m.pdf");
  c10->Print("s_m.svg");
  c10->Print("s_m.root");
  //==========================================================
  
  //===========================================================
  //get MC sum of weight
  double w_MC = 0;
  w_MC = hist2_11->Integral() + hist3_11->Integral();
  TFile* pileup_central = TFile::Open("/afs/cern.ch/user/x/xuyan/WGProj/PROD17/DATA/pileup/Pileup_17.root", "READ");
  TH1F* pileup_Data_central = (TH1F*)pileup_central->Get("pileup");
  pileup_Data_central->Scale(w_MC/(double)pileup_Data_central->Integral());
  pileup_Data_central->SetLineColor(1);
  pileup_Data_central->SetMarkerStyle(8);
  TCanvas *c11 = new TCanvas("c11","",2400,2200);
  c11->cd();
  TPad *p11a = new TPad("p11a","p11a",0.1,0.30,0.9,1.0);
  TPad *p11b = new TPad("p11b","p11b",0.1,0.1,0.9,0.315);
  p11a->Draw();
  p11b->Draw();
  p11a->cd();
  p11a->SetBottomMargin(0.11);
  p11a->SetLogy();
  xaxis1 = pileup_Data_central->GetXaxis();
  yaxis1 = pileup_Data_central->GetYaxis();
  xaxis1->SetTitle("PV");
  yaxis1->SetTitle("Entries / 1");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,100000000);
  pileup_Data_central->SetMarkerSize(2);
  pileup_Data_central->Draw("E1");
  stack11->Draw("SAMEHIST");
  pileup_Data_central->Draw("E1SAME");
  pileup_Data_central->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(pileup_Data_central,"2018 Data","lep");
  legend->AddEntry(hist2_11,"2018 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_11,"2018 MC, QCD(weighted)","f");
  legend->Draw();

  p11b->cd();
  p11b->SetTopMargin(0.037);
  p11b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)pileup_Data_central->Clone();
  bkg = (TH1*)hist2_11->Clone();
  bkg->Add((TH1*)hist3_11->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("PV");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p11b->SetGrid();
  p11b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c11->Print("s_PV.png");
  c11->Print("s_PV.pdf");
  c11->Print("s_PV.svg");
  c11->Print("s_PV.root");
  //==========================================================
}