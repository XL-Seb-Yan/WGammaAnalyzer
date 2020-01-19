void plotEff(){
  gROOT->SetBatch(1);
  
  //double narrow[15] = {0.1365,0.1346,0.1385,0.1352,0.1312,0.1324,0.1296,0.1233,0.1236,0.1190,0.1145,0.1065,0.1085,0.1092,0.1058};
  double wide[15] = {0.0963,0.0951,0.0889,0.0829,0.0658,0.0524,0.0408,0.0320,0.0244,0.0177,0.0140,0.0102,0.0077,0.0060,0.0029};
  
  double mass[15] = {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500};

  TGraph *gr1 = new TGraph(15,&mass[0],&wide[0]);
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(21);
  TF1 *f1 = new TF1("fun","landau(0)-[3]*((x-[4])/13000)^2",500,3500);
  f1->SetParameters(0.53,810,258,0,10000);
  gr1->Fit(f1,"F");
    //gr1->Fit("landau","F");
    //TF1 *f1 = gr1->GetFunction("landau");
  f1->SetLineColor(4);
  
  TCanvas *c0 = new TCanvas("c0","",1200,900);
  c0->cd();
  gr1->SetTitle("Total Efficiency (Wide)");
  TAxis *xaxis = gr1->GetXaxis();
  TAxis *yaxis = gr1->GetYaxis();
  xaxis->SetRangeUser(500,3500);
  xaxis->SetTitle("m_{X} (GeV)");
  c0->SetLogy();
  yaxis->SetRangeUser(0.001,1);
  yaxis->SetTitle("Total Efficiency");
  gr1->Draw("AP");
  f1->Draw("SAME");
  c0->Print("EffWide.png");
}
