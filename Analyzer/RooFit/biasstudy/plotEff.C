#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
void plotEff(){
  gROOT->SetBatch(1);
  lumi_13TeV = "35.92 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 4;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.025);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(3);

  //16 madgraph
  //double narrow[13] = {0.1579,0.1562,0.1571,0.1537,0.1559,0.1530,0.1477,0.1343,0.1379,0.1299,0.1291,0.1250,0.1204};
  //double wide[15] = {0.1473,0.1488,0.1500,0.1514,0.1466,0.1448,0.1387,0.1389,0.1265,0.1265,0.1257,0.1183,0.1106,0.1101,0.1041};
  
  //17 pythia
  //double narrow[15] = {0.1365,0.1346,0.1385,0.1352,0.1312,0.1324,0.1296,0.1233,0.1236,0.1190,0.1145,0.1065,0.1085,0.1092,0.1058};
  //double wide[15] = {0.0963,0.0951,0.0889,0.0829,0.0658,0.0524,0.0408,0.0320,0.0244,0.0177,0.0140,0.0102,0.0077,0.0060,0.0029};

  //17 madrgaph
  //double narrow[14] = {0.1645,0.1680,0.1690,0.1717,0.1665,0.1635,0.1619,0.1475,0.1429,0.1400,0.1338,0.1382,0.1300,0.1283};
  //double wide[14] = {0.1626,0.1622,0.1655,0.1698,0.1559,0.1570,0.1491,0.1498,0.1387,0.1357,0.1300,0.1249,0.1200,0.1062};
  
  //double massn[14] = {700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3500};
  //double massw[14] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3500};

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
  gr1->SetLineWidth(2);
  //gr1->SetMarkerSize(2);
  TF1 *f1 = new TF1("fun1","pol4",500,3500);
  //f1->SetParameters(0.125111,8.90778*pow(10,-5),-7.9384*pow(10,-8),2.46732*pow(10,-11),-2.7542*pow(10,-15));
  f1->SetLineColor(2);
  f1->SetLineStyle(2);
  gr1->Fit(f1,"F");
  f1->SetLineColor(2);
  f1->SetLineStyle(2);
  f1->SetLineWidth(2);

  TGraph *gr2 = new TGraph(5,&massw[0],&wide[0]);
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(43);
  gr2->SetLineWidth(2);
  //gr1->SetMarkerSize(2);
  TF1 *f2 = new TF1("fun2","pol4",500,3500);
  //f2->SetParameters(0.113995,8.94788*pow(10,-5),-7.06752*pow(10,-8),1.84053*pow(10,-11),-1.64397*pow(10,-15));
  f2->SetLineColor(4);
  f2->SetLineStyle(2);
  gr2->Fit(f2,"F");
  f2->SetLineColor(4);
  f2->SetLineStyle(2);
  f1->SetLineWidth(2);

  TLegend *leg = new TLegend(0.15,0.15,0.30,0.25);
  leg->AddEntry(gr1,"FullSIM narrow","lp");
  leg->AddEntry(gr1,"FullSIM wide","lp");
  TCanvas *c0 = new TCanvas("c0","",1200,900);
  c0->cd();
  c0->SetBottomMargin(0.11);
  c0->SetLeftMargin(0.13);
  gr1->SetTitle("Total Efficiency");
  TAxis *xaxis = gr1->GetXaxis();
  TAxis *yaxis = gr1->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  c0->SetLogy();
  xaxis->SetLimits(500,3500); 
  yaxis->SetRangeUser(0.001,1);
  yaxis->SetTitle("Total Efficiency");
  yaxis->SetTitleOffset(1.1);
  f1->Draw("APL");
  gr1->Draw("AP");
  f1->Draw("APLSAME");
  CMS_lumi(c0,iPeriod,iPos);
  leg->Draw();
  c0->Print("EffN17S1.png");
  c0->Print("EffN17S1.pdf");
  c0->Print("EffN17S1.root");
  c0->Print("EffN17S1.svg");
}
