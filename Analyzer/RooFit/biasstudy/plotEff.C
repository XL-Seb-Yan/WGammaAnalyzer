#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
void plotEff(){
  gROOT->SetBatch(1);
  lumi_13TeV = "41.53fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 5;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.025);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(3);

  //16 madrgaph
  // double narrow[15] = {0.15595,0.15545,0.1574,0.15345,0.15465,0.15245,0.1474,0.1408,0.13385,0.1378,0.1296,0.1284,0.1224,0.121,0.1225};
  // double wide[14] = {0.14575,0.14786,0.1489,0.1508,0.1451,0.14405,0.1374,0.1383,0.1265,0.12637,0.12565,0.11244,0.10935,0.1043};  
  // double massn[15] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500};
  // double massw[14] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2800,3000,3500};
  
  //17 pythia
  //double narrow[15] = {0.1365,0.1346,0.1385,0.1352,0.1312,0.1324,0.1296,0.1233,0.1236,0.1190,0.1145,0.1065,0.1085,0.1092,0.1058};
  //double wide[15] = {0.0963,0.0951,0.0889,0.0829,0.0658,0.0524,0.0408,0.0320,0.0244,0.0177,0.0140,0.0102,0.0077,0.0060,0.0029};

  //17 madrgaph
  // double narrow[14] = {0.1604,0.1635,0.1647,0.16915,0.1636,0.16085,0.1589,0.1461,0.14065,0.13755,0.13215,0.13525,0.12675,0.12615};
  // double wide[14] = {0.15835,0.15715,0.16165,0.1652,0.15315,0.15405,0.14605,0.1465,0.1366,0.13415,0.12685,0.1229,0.11805,0.1044};
  // double massn[14] = {700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3500};
  // double massw[14] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3500};

  //17 spin1
  double narrow[5] = {0.1847,0.2071,0.1972,0.1784,0.1554};
  double wide[5] = {0.1533,0.1740,0.1806,0.1687,0.1651};
  double massn[5] = {700,1200,2000,2800,3500};
  double massw[5] = {700,1200,2000,2800,3500};

  //17 spin0spin1 combine
  //double narrow[5] = {0.1746,0.1868,0.1724,0.1583,0.1419};
  //double wide[5] = {0.1580,0.1650,0.1596,0.1444,0.1357};
  //double mass[5] = {700,1200,2000,2800,3500};
  
  
  TGraph *gr1 = new TGraph(5,&massn[0],&narrow[0]);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(43);
  //gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(3);
  gr1->SetLineWidth(2);
  TF1 *f1 = new TF1("fun1","pol4",500,3500);
  f1->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  f1->SetLineColor(2);
  f1->SetLineStyle(2);
  gr1->Fit(f1,"F");
  gr1->Fit(f1,"F");
  f1->SetLineColor(2);
  f1->SetLineStyle(2);
  f1->SetLineWidth(2);

  TGraph *gr2 = new TGraph(5,&massw[0],&wide[0]);
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(43);
  //gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(3);
  gr2->SetLineWidth(2);
  TF1 *f2 = new TF1("fun2","pol4",500,3500);
  f2->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11),-2.86354*pow(10,-15));
  f2->SetLineColor(4);
  f2->SetLineStyle(2);
  gr2->Fit(f2,"F");
  gr2->Fit(f2,"F");
  f2->SetLineColor(4);
  f2->SetLineStyle(2);
  f2->SetLineWidth(2);
  
  // double STDN = 0;
  // for(int i=0; i<14; i++){
	  // cout<<(narrow[i] - f1->Eval(massn[i])) / narrow[i]<<",";
	  // STDN += pow((narrow[i] - f1->Eval(massn[i])) / narrow[i],2);
  // }
  // cout<<endl;
  // STDN = sqrt(STDN);
  // cout<<STDN<<endl;
  
  // double STDW = 0;
  // for(int i=0; i<15; i++){
	  // cout<<(wide[i] - f2->Eval(massw[i])) / wide[i]<<",";
	  // STDW += pow((wide[i] - f2->Eval(massw[i])) / wide[i],2);
  // }
  // cout<<endl;
  // STDW = sqrt(STDW);
  // cout<<STDW<<endl;

  TLegend *leg = new TLegend(0.15,0.15,0.35,0.25);
  leg->AddEntry(gr1,"2017 FullSIM narrow","lp");
  leg->AddEntry(gr1,"2017 FullSIM wide","lp");
  TCanvas *c0 = new TCanvas("c0","",2400,1800);
  c0->cd();
  c0->SetBottomMargin(0.11);
  c0->SetLeftMargin(0.13);
  gr1->SetTitle("Total Efficiency");
  TAxis *xaxis = gr1->GetXaxis();
  TAxis *yaxis = gr1->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  c0->SetLogy();
  xaxis->SetLimits(500,3500); 
  yaxis->SetRangeUser(0.01,1);
  yaxis->SetTitle("Efficiency #times Acceptance");
  yaxis->SetTitleOffset(1.1);
  f1->Draw("APL");
  gr1->Draw("AP");
  gr2->Draw("SAMEP");
  f1->Draw("APLSAME");
  CMS_lumi(c0,iPeriod,iPos);
  leg->Draw();
  c0->Print("Eff17spin-1.png");
  c0->Print("Eff17spin-1.pdf");
  c0->Print("Eff17spin-1.root");
  c0->Print("Eff17spin-1.svg");
}
