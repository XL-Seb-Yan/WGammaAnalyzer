void plotEff(){
  gROOT->SetBatch(1);
  
  double narrow[15] = {0.1365,0.1346,0.1385,0.1352,0.1312,0.1324,0.1296,0.1233,0.1236,0.1190,0.1145,0.1065,0.1085,0.1092,0.1058};
  //double wide[13] = {0.0992,0.0997,0.0951,0.0893,0.0745,0.0625,0.0534,0.0409,0.0330,0.0279,0.0242,0.0219,0.0155};
  
  double mass[15] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500};

  TGraph *gr1 = new TGraph(15,&mass[0],&narrow[0]);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(21);
  gr1->Fit("landau","F");
  TF1 *f1 = gr1->GetFunction("landau");
  f1->SetLineColor(2);
  
  TCanvas *c0 = new TCanvas("c0","",1200,900);
  c0->cd();
  gr1->SetTitle("Total Efficiency (Narrow)");
  TAxis *xaxis = gr1->GetXaxis();
  TAxis *yaxis = gr1->GetYaxis();
  xaxis->SetRangeUser(500,3500);
  xaxis->SetTitle("m_{X} (GeV)");
  c0->SetLogy();
  yaxis->SetRangeUser(0.001,1);
  yaxis->SetTitle("Total Efficiency");
  gr1->Draw("AP");
  f1->Draw("SAME");
  c0->Print("EffNarrow.png");
}
