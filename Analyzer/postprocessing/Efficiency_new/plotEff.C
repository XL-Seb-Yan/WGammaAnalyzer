#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
void plotEff(){
  gROOT->SetBatch(1);
  lumi_13TeV = "";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.45;
  cmsTextSize = 0.55;
  int iPeriod = 5;
  int iPos = 0;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetHistLineWidth(3);
  gStyle->SetEndErrorSize(5);
  
  std::vector<float> massn17;
  std::vector<float> massw17;
  std::vector<float> effn17;
  std::vector<float> effw17;
  std::vector<float> effn_err17;
  std::vector<float> effw_err17;
  
  std::vector<float> massn17s1;
  std::vector<float> massw17s1;
  std::vector<float> effn17s1;
  std::vector<float> effw17s1;
  std::vector<float> effn_err17s1;
  std::vector<float> effw_err17s1;
  
  ifstream file17n("EffS0/log17N_eff_new.txt");
  ifstream file17w("EffS0/log17W_eff_new.txt");
  ifstream file17ns1("EffS1/log17NS1_newmpoints.txt");
  ifstream file17ws1("EffS1/log17WS1_newmpoints.txt");
  string str; 
  
  //These efficiencies are correct to reshow a similar behavior of reconstruction efficiencies
  float MCmassS0[20]={700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500,4000,5000,6000,7000,8000};
  float EffN[20]={0.1077,0.1076,0.1149,0.1175,0.1216,0.1218,0.1236,0.1249,0.1239,0.1197,0.1212,0.1185,0.1198,0.1182,0.1168,0.1140,0.1085,0.1069,0.0959,0.0907};
  for(int i=0; i<20; i++){
      massn17.push_back(MCmassS0[i]);
      effn17.push_back(EffN[i]);
      effn_err17.push_back(sqrt(EffN[i]*(1-EffN[i])/20000));
  }
  float EffW[20]={0.0989,0.1135,0.1110,0.1150,0.1198,0.1214,0.1240,0.1256,0.1206,0.1209,0.1231,0.1200,0.1181,0.1147,0.1063,0.1037,0.0993,0.0794,0.0719,0.0624};
  for(int i=0; i<20; i++){
      massw17.push_back(MCmassS0[i]);
      effw17.push_back(EffW[i]);
      effw_err17.push_back(sqrt(EffW[i]*(1-EffW[i])/20000));
  }
  
  float MCmassS1[8]={700,1500,2000,2500,3000,4000,6000,8000};
  float EffNS1[8]={0.136001, 0.164015, 0.163318, 0.159733, 0.152397, 0.144864, 0.131334, 0.107218};
  for(int i=0; i<8; i++){
      massn17s1.push_back(MCmassS1[i]);
      effn17s1.push_back(EffNS1[i]);
      effn_err17s1.push_back(sqrt(EffNS1[i]*(1-EffNS1[i])/20000));
  }
  float EffWS1[8]={0.112669, 0.152843, 0.157272, 0.151754, 0.144713, 0.137764, 0.126129, 0.100312};
  for(int i=0; i<8; i++){
      massw17s1.push_back(MCmassS1[i]);
      effw17s1.push_back(EffWS1[i]);
      effw_err17s1.push_back(sqrt(EffWS1[i]*(1-EffWS1[i])/20000));
  }
  
  cout<<"=========================Fitting 17==========================="<<endl;
  TGraphErrors *gr3 = new TGraphErrors(massn17.size(),&massn17[0],&effn17[0],0,&effn_err17[0]);
  gr3->SetMarkerColor(kPink+9);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(2);
  gr3->SetLineWidth(2);
  gr3->SetLineColor(kPink+9);
  TF1 *f3 = new TF1("fun3","pol4",500,8200);
  f3->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),0,0);
  f3->SetLineColor(kPink+9);
  gr3->Fit(f3,"R");
  gr3->Fit(f3,"R");
  f3->SetLineColor(kPink+9);
  f3->SetLineWidth(2);
  cout<<f3->GetParameter(0)<<","<<f3->GetParameter(1)<<","<<f3->GetParameter(2)<<","<<f3->GetParameter(3)<<","<<f3->GetParameter(4)<<endl;
  TGraphErrors *gr3_unc = new TGraphErrors(365);
  for(int i=0; i<385; i++){
    gr3_unc->SetPoint(i, 500+i*20, 0);
  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr3_unc);
  gr3_unc->SetFillStyle(3002);
  gr3_unc->SetFillColorAlpha(6,0.5);
  
  TGraphErrors *gr4 = new TGraphErrors(massw17.size(),&massw17[0],&effw17[0],0,&effw_err17[0]);
  gr4->SetMarkerColor(kBlue);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(2);
  gr4->SetLineWidth(2);
  gr4->SetLineColor(kBlue);
  TF1 *f4 = new TF1("fun3","pol4",500,8200);
  f4->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),0,0);
  f4->SetLineColor(kBlue);
  gr4->Fit(f4,"R");
  gr4->Fit(f4,"R");
  f4->SetLineColor(kBlue);
  f4->SetLineWidth(2);
  cout<<f4->GetParameter(0)<<","<<f4->GetParameter(1)<<","<<f4->GetParameter(2)<<","<<f4->GetParameter(3)<<","<<f4->GetParameter(4)<<endl;
  TGraphErrors *gr4_unc = new TGraphErrors(365);
  for(int i=0; i<385; i++){
    gr4_unc->SetPoint(i, 500+i*20, 0);
  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr4_unc);
  gr4_unc->SetFillStyle(3002);
  gr4_unc->SetFillColorAlpha(kBlue,0.5);
  
  cout<<endl;
  cout<<endl;
  cout<<"=========================Fitting 17 S1==========================="<<endl;
  TGraphErrors *gr5 = new TGraphErrors(8,&massn17s1[0],&effn17s1[0],0,&effn_err17s1[0]);
  gr5->SetMarkerColor(kPink+9);
  gr5->SetMarkerStyle(21);
  gr5->SetMarkerSize(2);
  gr5->SetLineWidth(2);
  gr5->SetLineColor(kPink+9);
  TF1 *f5 = new TF1("fun5","pol4",500,8200);
  f5->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11));
  f5->SetLineColor(kPink+9);
  f5->SetLineStyle(7);
  gr5->Fit(f5,"R");
  gr5->Fit(f5,"R");
  f5->SetLineColor(kPink+9);
  f5->SetLineStyle(2);
  f5->SetLineWidth(2);
  cout<<f5->GetParameter(0)<<","<<f5->GetParameter(1)<<","<<f5->GetParameter(2)<<","<<f5->GetParameter(3)<<","<<f5->GetParameter(4)<<endl;
  TGraphErrors *gr5_unc = new TGraphErrors(365);
  for(int i=0; i<385; i++){
    gr5_unc->SetPoint(i, 500+i*20, 0);
  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr5_unc);
  gr5_unc->SetFillStyle(3002);
  gr5_unc->SetFillColorAlpha(6,0.5);

  TGraphErrors *gr6 = new TGraphErrors(8,&massw17s1[0],&effw17s1[0],0,&effw_err17s1[0]);
  gr6->SetMarkerColor(kBlue);
  gr6->SetMarkerStyle(21);
  gr6->SetMarkerSize(2);
  gr6->SetLineColor(kBlue);
  gr6->SetLineWidth(2);
  TF1 *f6 = new TF1("fun6","pol4",500,8200);
  f6->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f6->SetLineColor(kBlue);
  f6->SetLineStyle(2);
  gr6->Fit(f6,"R");
  gr6->Fit(f6,"R");
  f6->SetLineColor(kBlue);
  f6->SetLineStyle(2);
  f6->SetLineWidth(2);
  cout<<f6->GetParameter(0)<<","<<f6->GetParameter(1)<<","<<f6->GetParameter(2)<<","<<f6->GetParameter(3)<<","<<f6->GetParameter(4)<<endl;
  TGraphErrors *gr6_unc = new TGraphErrors(365);
  for(int i=0; i<385; i++){
    gr6_unc->SetPoint(i, 500+i*20, 0);
  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr6_unc);
  gr6_unc->SetFillStyle(3002);
  gr6_unc->SetFillColorAlpha(kBlue,0.5);

  TLegend *leg = new TLegend(0.3,0.65,0.85,0.75);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;
  
  leg = new TLegend(0.3,0.65,0.85,0.75);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  TCanvas *c1 = new TCanvas("c1","",2400,1800);
  c1->cd();
  c1->SetBottomMargin(0.11);
  c1->SetLeftMargin(0.13);
  gr3->SetTitle("A #varepsilon");
  xaxis = gr3->GetXaxis();
  yaxis = gr3->GetYaxis();
  xaxis->SetTitle("m_{X} [GeV]");
  xaxis->SetLimits(500,8200); 
  yaxis->SetRangeUser(0,0.3);
  yaxis->SetTitle("A #varepsilon");
  yaxis->SetTitleOffset(1.2);
  
  f3->Draw("APL");
  f4->Draw("APLSAME");
  f5->Draw("APLSAME");
  f6->Draw("APLSAME");
  gr3->Draw("AP");
  gr3_unc->Draw("E3SAME");
  gr4->Draw("SAMEP");
  gr4_unc->Draw("E3SAME");
  gr5->Draw("SAMEP");
  gr5_unc->Draw("E3SAME");
  gr6->Draw("SAMEP");
  gr6_unc->Draw("E3SAME");
  f3->Draw("APLSAME");
  f4->Draw("APLSAME");
  f5->Draw("APLSAME");
  f6->Draw("APLSAME");
  lumi_13TeV = "";
  CMS_lumi(c1,5,iPos);
  leg->Clear();
  leg->AddEntry(gr3,"Spin-0, #Gamma_{X} / m_{X} = 0.01%","lep");
  leg->AddEntry(gr4,"Spin-0, #Gamma_{X} / m_{X} = 5%","lep");
  leg->AddEntry(gr5,"Spin-1, #Gamma_{X} / m_{X} = 0.01%","lep");
  leg->AddEntry(gr6,"Spin-1, #Gamma_{X} / m_{X} = 5%","lep");
  leg->Draw();
  //c1->SetGrid();
  //c1->SetLogx();
  c1->Print("EffAcc17.pdf");
  c1->Print("EffAcc17.svg");
  
 
  
  //17 S-0 Acceptance (new scalar model)
  float AccN[20]={0.2984,0.3112,0.3186,0.3301,0.333,0.34,0.3368,0.3385,0.3393,0.3374,0.3362,0.3325,0.3325,0.3301,0.327,0.318,0.312,0.296,0.287,0.28};
  float AccN_err[20];
  for(int i=0; i<20; i++){
      AccN_err[i]=sqrt(AccN[i]*(1-AccN[i])/20000);
  }
  float AccW[20]={0.2764,0.2938,0.3039,0.3155,0.326,0.3247,0.3382,0.3364,0.335,0.3317,0.3281,0.3198,0.3209,0.3139,0.297,0.2886,0.2707,0.2418,0.2093,0.1865};
  float AccW_err[20];
  for(int i=0; i<20; i++){
      AccW_err[i]=sqrt(AccW[i]*(1-AccW[i])/20000);
  }
  
  //17 S-1 Acceptance
  float MCmassNS1[8]={700,1500,2000,2500,3000,4000,6000,8000};
  float AccNS1[8]={0.3924,0.4332,0.4531,0.4473,0.4302,0.421,0.3838,0.3638};
  float AccNS1_err[8];
  for(int i=0; i<8; i++){
      int Ngen = 20000;
      if(abs(MCmassNS1[i]-1500) < 1) //M-4000 only has 7000 events
        Ngen = 13000;
      if(abs(MCmassNS1[i]-2500) < 1) //M-4000 only has 7000 events
        Ngen = 19000;
      if(abs(MCmassNS1[i]-6000) < 1) //M-4000 only has 7000 events
        Ngen = 19000;
      AccNS1_err[i]=sqrt(AccNS1[i]*(1-AccNS1[i])/Ngen);
  }
  float MCmassWS1[8]={700,1500,2000,2500,3000,4000,6000,8000};
  float AccWS1[8]={0.3501,0.4159,0.440,0.4309,0.41245,0.4085,0.3777,0.3355};
  float AccWS1_err[8];
  for(int i=0; i<8; i++){
      int Ngen = 20000;
      if(abs(MCmassWS1[i]-1500) < 1) //M-4000 only has 7000 events
        Ngen = 13000;
      if(abs(MCmassWS1[i]-2500) < 1) //M-4000 only has 7000 events
        Ngen = 19000;
      if(abs(MCmassWS1[i]-6000) < 1) //M-4000 only has 7000 events
        Ngen = 19000;
      AccWS1_err[i]=sqrt(AccWS1[i]*(1-AccWS1[i])/Ngen);
  }
  
  cout<<"=========================Fitting 17 Acc ==========================="<<endl;
  TGraphErrors *gr7 = new TGraphErrors(20,MCmassS0,AccN,0,AccN_err);
  gr7->SetMarkerColor(kPink+9);
  gr7->SetMarkerStyle(20);
  gr7->SetMarkerSize(2);
  gr7->SetLineColor(kPink+9);
  gr7->SetLineWidth(2);
  TF1 *f7 = new TF1("fun7","pol4",500,8200);
  f7->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f7->SetLineColor(kPink+9);
  gr7->Fit(f7,"R");
  gr7->Fit(f7,"R");
  f7->SetLineColor(kPink+9);
  f7->SetLineWidth(2);
  TGraphErrors *gr7_unc = new TGraphErrors(365);
  for(int i=0; i<385; i++){
    gr7_unc->SetPoint(i, 500+i*20, 0);
  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr7_unc);
  gr7_unc->SetFillStyle(3002);
  gr7_unc->SetFillColorAlpha(6,0.5);
  cout<<f7->GetParameter(0)<<","<<f7->GetParameter(1)<<","<<f7->GetParameter(2)<<","<<f7->GetParameter(3)<<","<<f7->GetParameter(4)<<endl;
  
  TGraphErrors *gr8 = new TGraphErrors(20,MCmassS0,AccW,0,AccW_err);
  gr8->SetMarkerColor(kBlue);
  gr8->SetMarkerStyle(20);
  gr8->SetMarkerSize(2);
  gr8->SetLineColor(kBlue);
  gr8->SetLineWidth(2);
  TF1 *f8 = new TF1("fun8","pol4",500,8200);
  f8->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f8->SetLineColor(kBlue);
  gr8->Fit(f8,"R");
  gr8->Fit(f8,"R");
  f8->SetLineColor(kBlue);
  f8->SetLineWidth(2);
  TGraphErrors *gr8_unc = new TGraphErrors(365);
  for(int i=0; i<385; i++){
    gr8_unc->SetPoint(i, 500+i*20, 0);
  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr8_unc);
  gr8_unc->SetFillStyle(3002);
  gr8_unc->SetFillColorAlpha(kBlue,0.5);
  cout<<f8->GetParameter(0)<<","<<f8->GetParameter(1)<<","<<f8->GetParameter(2)<<","<<f8->GetParameter(3)<<","<<f8->GetParameter(4)<<endl;
  
  cout<<endl;
  cout<<endl;
  cout<<"=========================Fitting 17 S1 Acc ==========================="<<endl;
  TGraphErrors *gr11 = new TGraphErrors(8,MCmassNS1,AccNS1,0,AccNS1_err);
  gr11->SetMarkerColor(kPink+9);
  gr11->SetMarkerStyle(21);
  gr11->SetMarkerSize(2);
  gr11->SetLineColor(kPink+9);
  gr11->SetLineWidth(2);
  TF1 *f11 = new TF1("fun11","pol4",500,8200);
  f11->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f11->SetLineColor(kPink+9);
  f11->SetLineStyle(2);
  gr11->Fit(f11,"R");
  gr11->Fit(f11,"R");
  f11->SetLineColor(kPink+9);
  f11->SetLineStyle(2);
  f11->SetLineWidth(2);
  TGraphErrors *gr11_unc = new TGraphErrors(365);
  for(int i=0; i<385; i++){
    gr11_unc->SetPoint(i, 500+i*20, 0);
  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr11_unc);
  gr11_unc->SetFillStyle(3002);
  gr11_unc->SetFillColorAlpha(6,0.5);
  cout<<f11->GetParameter(0)<<","<<f11->GetParameter(1)<<","<<f11->GetParameter(2)<<","<<f11->GetParameter(3)<<","<<f11->GetParameter(4)<<endl;
  
  TGraphErrors *gr12 = new TGraphErrors(8,MCmassWS1,AccWS1,0,AccWS1_err);
  gr12->SetMarkerColor(kBlue);
  gr12->SetMarkerStyle(21);
  gr12->SetMarkerSize(2);
  gr12->SetLineColor(kBlue);
  gr12->SetLineWidth(2);
  TF1 *f12 = new TF1("fun12","pol4",500,8200);
  f12->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f12->SetLineColor(kBlue);
  f12->SetLineStyle(2);
  gr12->Fit(f12,"R");
  gr12->Fit(f12,"R");
  f12->SetLineColor(kBlue);
  f12->SetLineStyle(2);
  f12->SetLineWidth(2);
  TGraphErrors *gr12_unc = new TGraphErrors(365);
  for(int i=0; i<385; i++){
    gr12_unc->SetPoint(i, 500+i*20, 0);
  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gr12_unc);
  gr12_unc->SetFillStyle(3002);
  gr12_unc->SetFillColorAlpha(kBlue,0.5);
  cout<<f12->GetParameter(0)<<","<<f12->GetParameter(1)<<","<<f12->GetParameter(2)<<","<<f12->GetParameter(3)<<","<<f12->GetParameter(4)<<endl;
  
  
  TCanvas *c2 = new TCanvas("c2","",2400,1800);
  c2->cd();
  c2->SetBottomMargin(0.11);
  c2->SetLeftMargin(0.13);
  gr7->SetTitle("A");
  xaxis = gr7->GetXaxis();
  yaxis = gr7->GetYaxis();
  xaxis->SetTitle("m_{X} [GeV]");
  xaxis->SetLimits(500,8200); 
  yaxis->SetRangeUser(0,0.7);
  yaxis->SetTitle("A");
  yaxis->SetTitleOffset(1.2);
  f7->Draw("APL");
  f8->Draw("APLSAME");
  f11->Draw("APLSAME");
  f12->Draw("APLSAME");
  gr7->Draw("AP");
  gr7_unc->Draw("E3SAME");
  gr8->Draw("SAMEP");
  gr8_unc->Draw("E3SAME");
  gr11->Draw("SAMEP");
  gr11_unc->Draw("E3SAME");
  gr12->Draw("SAMEP");
  gr12_unc->Draw("E3SAME");
  f7->Draw("APLSAME");
  f8->Draw("APLSAME");
  f11->Draw("APLSAME");
  f12->Draw("APLSAME");
  lumi_13TeV = "";
  CMS_lumi(c2,5,iPos);
  leg->Clear();
  leg->AddEntry(gr7,"Spin-0, #Gamma_{X} / m_{X} = 0.01%","lep");
  leg->AddEntry(gr8,"Spin-0, #Gamma_{X} / m_{X} = 5%","lep");
  leg->AddEntry(gr11,"Spin-1, #Gamma_{X} / m_{X} = 0.01%","lep");
  leg->AddEntry(gr12,"Spin-1, #Gamma_{X} / m_{X} = 5%","lep");
  leg->Draw();
  //c2->SetGrid();
  //c2->SetLogx();
  c2->Print("Acc17.pdf");
  c2->Print("Acc17.svg");
  
  TGraph *gr9 = new TGraph();
  TGraph *gr10 = new TGraph();
  TGraph *gr13 = new TGraph();
  TGraph *gr14 = new TGraph();
  
  int ipoints = 0;
  for(int i=700; i<8001; i+=10){
      gr9->SetPoint(ipoints,i,f3->Eval(i)/f7->Eval(i));
      gr10->SetPoint(ipoints,i,f4->Eval(i)/f8->Eval(i));
      gr13->SetPoint(ipoints,i,f5->Eval(i)/f11->Eval(i));
      gr14->SetPoint(ipoints,i,f6->Eval(i)/f12->Eval(i));
      ipoints++;
  }
  
  TCanvas *c3 = new TCanvas("c3","",2400,1800);
  c3->cd();
  c3->SetBottomMargin(0.11);
  c3->SetLeftMargin(0.13);
  xaxis = gr9->GetXaxis();
  yaxis = gr9->GetYaxis();
  xaxis->SetTitle("m_{X} [GeV]");
  xaxis->SetLimits(500,8200); 
  yaxis->SetRangeUser(0.0,0.6);
  yaxis->SetTitle("#varepsilon");
  yaxis->SetTitleOffset(1.2);
  gr9->SetLineColor(kPink+9);
  gr10->SetLineColor(kBlue);
  gr13->SetLineColor(kPink+9);
  gr14->SetLineColor(kBlue);
  gr9->SetLineWidth(2);
  gr10->SetLineWidth(2);
  gr13->SetLineWidth(2);
  gr14->SetLineWidth(2);
  gr13->SetLineStyle(2);
  gr14->SetLineStyle(2);
  gr9->Draw("AL");
  gr10->Draw("SAMEL");
  gr13->Draw("SAMEL");
  gr14->Draw("SAMEL");
  lumi_13TeV = "";
  CMS_lumi(c3,5,iPos);
  leg->Clear();
  leg->AddEntry(gr9,"Narrow S-0","l");
  leg->AddEntry(gr10,"Broad S-0","l");
  leg->AddEntry(gr13,"Narrow S-1","l");
  leg->AddEntry(gr14,"Broad S-1","l");
  leg->Draw();
  //c3->SetGrid();
  //c3->SetLogx();
  c3->Print("Eff17.pdf");
  c3->Print("Eff17.svg");
  
  /*
  // Yield calculation
  cout<<"S-0 Acc*Eff: ========================================="<<endl;
  for(int i=700; i<8001; i+=10){
    double mass = i;
    double yield = 137.19*1*f3->Eval(mass);
    cout<<"N"<<mass<<" "<<yield<<" ;"<<endl;
    // cout<<mass<<" "<<yield<<" "<<limit*yield*3<<" "<<limit*yield*10<<endl;
  }
  // Yield calculation
  for(int i=700; i<8001; i+=10){
    double mass = i;
    double yield = 137.19*1*f4->Eval(mass);
    cout<<"W"<<mass<<" "<<yield<<" ;"<<endl;
    // cout<<mass<<" "<<yield<<" "<<limit*yield*3<<" "<<limit*yield*10<<endl;
  }
  */
  /*
  cout<<"S-1 Acc*Eff: ========================================="<<endl;
  // Yield calculation
  for(int i=700; i<8001; i+=10){
    double mass = i;
    double yield = 137.19*1*f5->Eval(mass);
    cout<<"N"<<mass<<" "<<yield<<" ;"<<endl;
    // cout<<mass<<" "<<yield<<" "<<limit*yield*3<<" "<<limit*yield*10<<endl;
  }
  // Yield calculation
  for(int i=700; i<8001; i+=10){
    double mass = i;
    double yield = 137.19*1*f6->Eval(mass);
    cout<<"W"<<mass<<" "<<yield<<" ;"<<endl;
    // cout<<mass<<" "<<yield<<" "<<limit*yield*3<<" "<<limit*yield*10<<endl;
  }
  */
  /*
  cout<<"S-0 Eff: ========================================="<<endl;
  // Yield calculation
  for(int i=1500; i<8001; i+=10){
    double mass = i;
    double yield = 137.19*1*f3->Eval(mass)/f7->Eval(mass);
    cout<<"N"<<mass<<" "<<yield<<" "<<endl;
    // cout<<mass<<" "<<yield<<" "<<limit*yield*3<<" "<<limit*yield*10<<endl;
  }
  // Yield calculation
  for(int i=1500; i<8001; i+=10){
    double mass = i;
    double yield = 137.19*1*f4->Eval(mass)/f8->Eval(mass);
    cout<<"W"<<mass<<" "<<yield<<" "<<endl;
    // cout<<mass<<" "<<yield<<" "<<limit*yield*3<<" "<<limit*yield*10<<endl;
  }
  */
  
  /*
  cout<<"S-1 Eff: ========================================="<<endl;
  // Yield calculation
  for(int i=1500; i<8001; i+=20){
    double mass = i;
    double yield = 137.19*1*f5->Eval(mass)/f11->Eval(mass);
    cout<<"N"<<mass<<" "<<yield<<" "<<endl;
    // cout<<mass<<" "<<yield<<" "<<limit*yield*3<<" "<<limit*yield*10<<endl;
  }
  // Yield calculation
  for(int i=1500; i<8001; i+=20){
    double mass = i;
    double yield = 137.19*1*f6->Eval(mass)/f12->Eval(mass);
    cout<<"W"<<mass<<" "<<yield<<" "<<endl;
    // cout<<mass<<" "<<yield<<" "<<limit*yield*3<<" "<<limit*yield*10<<endl;
  }
  */
  
    for(int i=0; i<385; i++){
      int mass = i*20+500;
      cout<<mass<<" "<<f3->Eval(mass)<<" "<<gr3_unc->GetErrorYhigh(i)<<" "<<gr3_unc->GetErrorYlow(i)<<" "<<f4->Eval(mass)<<" "<<gr4_unc->GetErrorYhigh(i)<<" "<<gr4_unc->GetErrorYlow(i)<<" "<<f5->Eval(mass)<<" "<<gr5_unc->GetErrorYhigh(i)<<" "<<gr5_unc->GetErrorYlow(i)<<" "<<f6->Eval(mass)<<" "<<gr6_unc->GetErrorYhigh(i)<<" "<<gr6_unc->GetErrorYlow(i)<<endl;
  }
  
    for(int i=0; i<385; i++){
      int mass = i*20+500;
      cout<<mass<<" "<<f7->Eval(mass)<<" "<<gr7_unc->GetErrorYhigh(i)<<" "<<gr7_unc->GetErrorYlow(i)<<" "<<f8->Eval(mass)<<" "<<gr8_unc->GetErrorYhigh(i)<<" "<<gr8_unc->GetErrorYlow(i)<<" "<<f11->Eval(mass)<<" "<<gr11_unc->GetErrorYhigh(i)<<" "<<gr11_unc->GetErrorYlow(i)<<" "<<f12->Eval(mass)<<" "<<gr12_unc->GetErrorYhigh(i)<<" "<<gr12_unc->GetErrorYlow(i)<<endl;
  }
  
}