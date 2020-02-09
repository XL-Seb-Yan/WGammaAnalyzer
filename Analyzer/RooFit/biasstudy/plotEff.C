void plotEff(){
  gROOT->SetBatch(1);

  //16 madgraph
  double narrow[13] = {0.1579,0.1562,0.1571,0.1537,0.1559,0.1530,0.1477,0.1343,0.1379,0.1299,0.1291,0.1250,0.1204};
  //double wide[15] = {0.1473,0.1488,0.1500,0.1514,0.1466,0.1448,0.1387,0.1389,0.1265,0.1265,0.1257,0.1183,0.1106,0.1101,0.1041};
  
  //17 pythia
  //double narrow[15] = {0.1365,0.1346,0.1385,0.1352,0.1312,0.1324,0.1296,0.1233,0.1236,0.1190,0.1145,0.1065,0.1085,0.1092,0.1058};
  //double wide[15] = {0.0963,0.0951,0.0889,0.0829,0.0658,0.0524,0.0408,0.0320,0.0244,0.0177,0.0140,0.0102,0.0077,0.0060,0.0029};

  //17 madrgaph
  //double narrow[14] = {0.1645,0.1680,0.1690,0.1717,0.1665,0.1635,0.1619,0.1475,0.1429,0.1400,0.1338,0.1382,0.1300,0.1283};
  //double wide[14] = {0.1626,0.1622,0.1655,0.1698,0.1559,0.1570,0.1491,0.1498,0.1387,0.1357,0.1300,0.1249,0.1200,0.1062};
  
  double mass[13] = {700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000};

  TGraph *gr1 = new TGraph(13,&mass[0],&narrow[0]);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(21);
  TF1 *f1 = new TF1("fun","pol4",500,3500);
  f1->SetParameters(0.125111,8.90778*pow(10,-5),-7.9384*pow(10,-8),2.46732*pow(10,-11),-2.7542*pow(10,-15));
  //f1->SetParameters(0.87978,807.827,1801.19);
  f1->SetLineColor(2);
  gr1->Fit(f1,"F");
    //gr1->Fit("landau","F");
    //TF1 *f1 = gr1->GetFunction("landau");
  f1->SetLineColor(2);
  
  TCanvas *c0 = new TCanvas("c0","",1200,900);
  c0->cd();
  gr1->SetTitle("Total Efficiency");
  TAxis *xaxis = gr1->GetXaxis();
  TAxis *yaxis = gr1->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  c0->SetLogy();
  xaxis->SetLimits(500,3500); 
  yaxis->SetRangeUser(0.001,1);
  yaxis->SetTitle("Total Efficiency");
  gr1->Draw("AP");
  f1->Draw("SAME");
  c0->Print("EffN.png");
  c0->Print("EffN.pdf");
}
