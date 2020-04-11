void Histofit()
{

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);

  // Write histos to root file
  TFile *file1 = TFile::Open("Histogram_Data.root");
  TFile *file2 = TFile::Open("Histogram_GJets.root");
  TFile *file3 = TFile::Open("Histogram_QCD.root");
  TH1* hist1_1 = (TH1*)file1->Get("1"); //p_pt
  TH1* hist1_2 = (TH1*)file1->Get("2"); //p_eta
  TH1* hist1_3 = (TH1*)file1->Get("3"); //j_pt
  TH1* hist1_4 = (TH1*)file1->Get("4"); //j_eta
  TH1* hist1_5 = (TH1*)file1->Get("5"); //j_e
  TH1* hist1_6 = (TH1*)file1->Get("6"); //j_mass
  TH1* hist1_7 = (TH1*)file1->Get("7"); //j_tau21
  TH1* hist1_8 = (TH1*)file1->Get("8"); //s_cos
  TH1* hist1_9 = (TH1*)file1->Get("9"); //s_ptm
  TH1* hist1_10 = (TH1*)file1->Get("10"); //s_invmass
  TH1* hist2_1 = (TH1*)file2->Get("1");
  TH1* hist2_2 = (TH1*)file2->Get("2");
  TH1* hist2_3 = (TH1*)file2->Get("3");
  TH1* hist2_4 = (TH1*)file2->Get("4");
  TH1* hist2_5 = (TH1*)file2->Get("5");
  TH1* hist2_6 = (TH1*)file2->Get("6");
  TH1* hist2_7 = (TH1*)file2->Get("7");
  TH1* hist2_8 = (TH1*)file2->Get("8");
  TH1* hist2_9 = (TH1*)file2->Get("9");
  TH1* hist2_10 = (TH1*)file2->Get("10");
  TH1* hist3_1 = (TH1*)file3->Get("1");
  TH1* hist3_2 = (TH1*)file3->Get("2");
  TH1* hist3_3 = (TH1*)file3->Get("3");
  TH1* hist3_4 = (TH1*)file3->Get("4");
  TH1* hist3_5 = (TH1*)file3->Get("5");
  TH1* hist3_6 = (TH1*)file3->Get("6");
  TH1* hist3_7 = (TH1*)file3->Get("7");
  TH1* hist3_8 = (TH1*)file3->Get("8");
  TH1* hist3_9 = (TH1*)file3->Get("9");
  TH1* hist3_10 = (TH1*)file3->Get("10");

  /*
  // Maximum likelihood fit result: GJets: 1.45, QCD: 0.01
  // Maximum likelihood fitting to invariant mass histograms
  TGraph2D *g = new TGraph2D();
  int NBins = hist1_10->GetNbinsX();
  cout<<"Number of bins: "<<NBins<<endl;
  int point = 0;
  double maxlh = -99999999999999;
  double sGJets = -99;
  double sQCD = -99;
  for(int i=1; i<150; i++){
    for(int j=1; j<150; j++){
      double lh = 0;
      double scaleGJets = i*0.01;
      double scaleQCD = j*0.01;
      for(int k=1; k<NBins+1; k++){
	double NData = hist1_10->GetBinContent(k);
	double NGjets = hist2_10->GetBinContent(k);
	double NQCD = hist3_10->GetBinContent(k);
	if(scaleGJets*NGjets + scaleQCD*NQCD == 0) continue; //exclude first several empty bins lower than trigger turn on
	lh += (NData*log(scaleGJets*NGjets + scaleQCD*NQCD) - (scaleGJets*NGjets + scaleQCD*NQCD));
      }
      if(lh > maxlh){
	maxlh = lh;
	sGJets = i*0.01;
	sQCD = j*0.01;
      }
      cout<<"Runnig through "<<i*0.01<<" "<<j*0.01<<" Likelihood is: "<<lh<<endl;
      g->SetPoint(point,scaleGJets,scaleQCD,lh);
      point++;
    }
  }
  */

  
  // Least Square fit result: GJets: 1.41, QCD: 0.57
  // Least Square fitting to invariant mass histograms
  TGraph2D *g = new TGraph2D();
  int NBins = hist1_10->GetNbinsX();
  cout<<"Number of bins: "<<NBins<<endl;
  int point = 0;
  double minls = 99999999999;
  double sGJets = -99;
  double sQCD = -99;
  for(int i=400; i<1000; i++){
    for(int j=0; j<150; j++){
      double ls = 0;
      double scaleGJets = i*0.002;
      double scaleQCD = j*0.01;
      for(int k=1; k<NBins+1; k++){
	double NData = hist1_10->GetBinContent(k);
	double NGjets = hist2_10->GetBinContent(k);
	double NQCD = hist3_10->GetBinContent(k);
	//cout<<NData<<" "<<NGjets<<" "<<NQCD<<" "<<pow((NData-(scaleGJets*NGjets + scaleQCD*NQCD))/NData,2)<<endl;
	if(NData == 0) continue;//exclude first several empty bins lower than trigger turn on
	ls += pow((NData-(scaleGJets*NGjets + scaleQCD*NQCD))/1000,2);
      }
      if(ls < minls){
	minls = ls;
	sGJets = scaleGJets;
	sQCD = scaleQCD;
      }
      cout<<"Runnig through "<<i*0.002<<" "<<j*0.01<<" Least Square is: "<<ls<<endl;
      g->SetPoint(point,scaleGJets,scaleQCD,ls);
      point++;
    }
  }
  

  cout<<sGJets<<" "<<sQCD<<" "<<minls<<endl;

  TCanvas *c01 = new TCanvas("c01","scale",1200,900);
  TAxis *xaxis = g->GetXaxis();
  TAxis *yaxis = g->GetYaxis();
  yaxis->SetTitleOffset(1.3);
  xaxis->SetTitleOffset(1.3);
  xaxis->SetTitle("GJets scale");
  yaxis->SetTitle("QCD scale");
  //yaxis->SetRangeUser(0.5,10000000);
  c01->SetLogz();
  c01->cd();
  c01->SetGrid();
  g->Draw("COLZ");
  cout<<"OK"<<endl;
  c01->Print("MCFit.png");
}
