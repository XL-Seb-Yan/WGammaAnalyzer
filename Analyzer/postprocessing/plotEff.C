#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
void plotEff(){
  gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
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
  gStyle->SetLegendTextSize(0.022);
  gStyle->SetHistLineWidth(3);
  
  std::vector<float> massn16;
  std::vector<float> massw16;
  std::vector<float> massn17;
  std::vector<float> massw17;
  std::vector<float> effn16;
  std::vector<float> effw16;
  std::vector<float> effn17;
  std::vector<float> effw17;
  std::vector<float> effn_err16;
  std::vector<float> effw_err16;
  std::vector<float> effn_err17;
  std::vector<float> effw_err17;
  
  std::vector<float> massn17s1;
  std::vector<float> massw17s1;
  std::vector<float> effn17s1;
  std::vector<float> effw17s1;
  std::vector<float> effn_err17s1;
  std::vector<float> effw_err17s1;
  
  ifstream file16n("log16N_eff.txt");
  ifstream file17n("log17N_eff_new.txt");
  ifstream file16w("log16W_eff.txt");
  ifstream file17w("log17W_eff_new.txt");
  ifstream file17ns1("log17NS1_v4.txt");
  ifstream file17ws1("log17WS1_v4.txt");
  string str; 
  while (getline(file16n,str)) {
   if(str.find("++++") != std::string::npos){
	   std::stringstream ss(str);
	   int index = -1;
	   while(ss.good()){
          string substr;
          getline(ss,substr,',');
		  index++;
          if(index == 0) continue;
		  if(index == 1) massn16.push_back(strtof((substr).c_str(),0));
		  if(index == 2){
			  effn16.push_back(strtof((substr).c_str(),0) / 20000);
			  effn_err16.push_back(sqrt(strtof((substr).c_str(),0)) / 20000);
		  }
       }
   }
  }
  while (getline(file16w,str)) {
   if(str.find("++++") != std::string::npos){
	   std::stringstream ss(str);
	   int index = -1;
	   while(ss.good()){
          string substr;
          getline(ss,substr,',');
		  index++;
          if(index == 0) continue;
		  if(index == 1) massw16.push_back(strtof((substr).c_str(),0));
		  if(index == 2){
			  effw16.push_back(strtof((substr).c_str(),0) / 20000);
			  effw_err16.push_back(sqrt(strtof((substr).c_str(),0)) / 20000);
		  }
       }
   }
  }
  while (getline(file17n,str)) {
   if(str.find("++++") != std::string::npos){
	   std::stringstream ss(str);
	   int index = -1;
	   while(ss.good()){
          string substr;
          getline(ss,substr,',');
		  index++;
          if(index == 0) continue;
		  if(index == 1) massn17.push_back(strtof((substr).c_str(),0));
		  if(index == 2){
			  int Ngen = 20000;
			  if(abs(massn17.back()-4000) < 1) //M-4000 only has 19000 events
				  Ngen = 14000;
			  if(abs(massn17.back()-5000) < 1) //M-4000 only has 19000 events
				  Ngen = 12500;
			  effn17.push_back(strtof((substr).c_str(),0) / Ngen);
			  effn_err17.push_back(sqrt(strtof((substr).c_str(),0)) / Ngen);
		  }
       }
   }
  }
  while (getline(file17w,str)) {
   if(str.find("++++") != std::string::npos){
	   std::stringstream ss(str);
	   int index = -1;
	   while(ss.good()){
          string substr;
          getline(ss,substr,',');
		  index++;
          if(index == 0) continue;
		  if(index == 1) massw17.push_back(strtof((substr).c_str(),0));
		  if(index == 2){
			  int Ngen = 20000;
			  if(abs(massw17.back()-4000) < 1) //M-4000 only has 7000 events
				  Ngen = 19500;
              if(abs(massw17.back()-5000) < 1) //M-4000 only has 7000 events
                Ngen = 19500;
			  effw17.push_back(strtof((substr).c_str(),0) / Ngen);
			  effw_err17.push_back(sqrt(strtof((substr).c_str(),0)) / Ngen);
		  }
       }
   }
  }
  
  while (getline(file17ns1,str)) {
   if(str.find("++++") != std::string::npos){
	   std::stringstream ss(str);
	   int index = -1;
	   while(ss.good()){
          string substr;
          getline(ss,substr,',');
		  index++;
          if(index == 0) continue;
		  if(index == 1) massn17s1.push_back(strtof((substr).c_str(),0));
		  if(index == 2){
			  int Ngen = 10000;
			  effn17s1.push_back(strtof((substr).c_str(),0) / Ngen);
			  effn_err17s1.push_back(sqrt(strtof((substr).c_str(),0)) / Ngen);
		  }
       }
   }
  }
  
  
   while (getline(file17ws1,str)) {
   if(str.find("++++") != std::string::npos){
	   std::stringstream ss(str);
	   int index = -1;
	   while(ss.good()){
          string substr;
          getline(ss,substr,',');
		  index++;
          if(index == 0) continue;
		  if(index == 1) massw17s1.push_back(strtof((substr).c_str(),0));
		  if(index == 2){
			  int Ngen = 7000;
              if(abs(massw17s1.back()-700) < 1) //M-4000 only has 7000 events
                Ngen = 6500;
              if(abs(massw17s1.back()-3500) < 1) //M-4000 only has 7000 events
                Ngen = 6500;
              if(abs(massw17s1.back()-6000) < 1) //M-4000 only has 7000 events
                Ngen = 6000;
			  effw17s1.push_back(strtof((substr).c_str(),0) / Ngen);
			  effw_err17s1.push_back(sqrt(strtof((substr).c_str(),0)) / Ngen);
		  }
       }
   }
  }
  cout<<"=========================Fitting 16==========================="<<endl;
  TGraphErrors *gr1 = new TGraphErrors(massn16.size(),&massn16[0],&effn16[0],0,&effn_err16[0]);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(2.8);
  gr1->SetLineWidth(2);
  gr1->SetLineColor(2);
  TF1 *f1 = new TF1("fun1","pol4",700,3500);
  f1->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  f1->SetLineColor(2);
  gr1->Fit(f1,"R");
  gr1->Fit(f1,"R");
  f1->SetLineColor(2);
  f1->SetLineWidth(2);
  
  TGraphErrors *gr2 = new TGraphErrors(massw16.size(),&massw16[0],&effw16[0],0,&effw_err16[0]);
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(2.8);
  gr2->SetLineWidth(2);
  gr2->SetLineColor(4);
  TF1 *f2 = new TF1("fun2","pol4",700,3500);
  f2->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  f2->SetLineColor(4);
  gr2->Fit(f2,"R");
  gr2->Fit(f2,"R");
  f2->SetLineColor(4);
  f2->SetLineWidth(2);
  
  cout<<"=========================Fitting 17==========================="<<endl;
  TGraphErrors *gr3 = new TGraphErrors(massn17.size(),&massn17[0],&effn17[0],0,&effn_err17[0]);
  gr3->SetMarkerColor(2);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(2.8);
  gr3->SetLineWidth(2);
  gr3->SetLineColor(2);
  // TF1 *f3 = new TF1("fun3","pol4",700,6000);
  // f3->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  TF1 *f3 = new TF1("fun3","pol4",700,5000);
  f3->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),0,0);
  f3->SetLineColor(2);
  gr3->Fit(f3,"R");
  gr3->Fit(f3,"R");
  f3->SetLineColor(2);
  f3->SetLineWidth(2);
  
  TGraphErrors *gr4 = new TGraphErrors(massw17.size(),&massw17[0],&effw17[0],0,&effw_err17[0]);
  gr4->SetMarkerColor(4);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(2.8);
  gr4->SetLineWidth(2);
  gr4->SetLineColor(4);
  // TF1 *f4 = new TF1("fun4","pol4",700,6000);
  // f4->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  TF1 *f4 = new TF1("fun3","pol4",700,5000);
  f4->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),0,0);
  f4->SetLineColor(4);
  gr4->Fit(f4,"R");
  gr4->Fit(f4,"R");
  f4->SetLineColor(4);
  f4->SetLineWidth(2);
  
  cout<<"=========================Fitting 17 S1==========================="<<endl;
  TGraphErrors *gr5 = new TGraphErrors(7,&massn17s1[0],&effn17s1[0],0,&effn_err17s1[0]);
  gr5->SetMarkerColor(2);
  gr5->SetMarkerStyle(21);
  //gr1->SetMarkerStyle(20);
  gr5->SetMarkerSize(2.8);
  gr5->SetLineWidth(2);
  gr5->SetLineColor(2);
  gr5->SetLineStyle(2);
  TF1 *f5 = new TF1("fun5","pol4",700,5000);
  f5->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11));
  f5->SetLineColor(2);
  f5->SetLineStyle(2);
  gr5->Fit(f5,"R");
  gr5->Fit(f5,"R");
  f5->SetLineColor(2);
  f5->SetLineStyle(2);
  f5->SetLineWidth(2);

  TGraphErrors *gr6 = new TGraphErrors(7,&massw17s1[0],&effw17s1[0],0,&effw_err17s1[0]);
  gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(21);
  //gr2->SetMarkerStyle(20);
  gr6->SetMarkerSize(2.8);
  gr6->SetLineColor(4);
  gr6->SetLineStyle(2);
  gr6->SetLineWidth(2);
  TF1 *f6 = new TF1("fun6","pol4",700,5000);
  f6->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f6->SetLineColor(4);
  f6->SetLineStyle(2);
  gr6->Fit(f6,"R");
  gr6->Fit(f6,"R");
  f6->SetLineColor(4);
  f6->SetLineStyle(2);
  f6->SetLineWidth(2);

  TLegend *leg = new TLegend(0.6,0.8,0.9,0.9);
  
  TCanvas *c0 = new TCanvas("c0","",2400,1800);
  c0->cd();
  c0->SetBottomMargin(0.11);
  c0->SetLeftMargin(0.13);
  gr1->SetTitle("Total Efficiency");
  TAxis *xaxis = gr1->GetXaxis();
  TAxis *yaxis = gr1->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  xaxis->SetLimits(700,5000); 
  yaxis->SetRangeUser(0,0.3);
  yaxis->SetTitle("Acc. #times Eff.");
  yaxis->SetTitleOffset(1.2);
  f1->Draw("APL");
  f1->Draw("APL");
  gr1->Draw("AP");
  gr2->Draw("SAMEP");
  f1->Draw("PLSAME");
  f2->Draw("PLSAME");
  lumi_13TeV = "35.92 fb^{-1}";
  CMS_lumi(c0,4,iPos);
  leg->AddEntry(gr1,"2016 FullSIM Narrow S-0","lep");
  leg->AddEntry(gr2,"2016 FullSIM Wide S-0","lep");
  leg->Draw();
  c0->SetGrid();
  c0->Print("EffAcc16.png");
  c0->Print("EffAcc16.pdf");
  c0->Print("EffAcc16.svg");
  
  leg = new TLegend(0.6,0.7,0.9,0.9);
  TCanvas *c1 = new TCanvas("c1","",2400,1800);
  c1->cd();
  c1->SetBottomMargin(0.11);
  c1->SetLeftMargin(0.13);
  gr3->SetTitle("Total Efficiency");
  xaxis = gr3->GetXaxis();
  yaxis = gr3->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  xaxis->SetLimits(700,5000); 
  yaxis->SetRangeUser(0,0.3);
  yaxis->SetTitle("Acc. #times Eff.");
  yaxis->SetTitleOffset(1.2);
  f3->Draw("APL");
  f4->Draw("APLSAME");
  f5->Draw("APLSAME");
  f6->Draw("APLSAME");
  gr3->Draw("AP");
  gr4->Draw("SAMEP");
  gr5->Draw("SAMEP");
  gr6->Draw("SAMEP");
  f3->Draw("APLSAME");
  f4->Draw("APLSAME");
  f5->Draw("APLSAME");
  f6->Draw("APLSAME");
  lumi_13TeV = "41.53 fb^{-1}";
  CMS_lumi(c1,5,iPos);
  leg->Clear();
  leg->AddEntry(gr3,"2017 FullSIM Narrow S-0","lep");
  leg->AddEntry(gr4,"2017 FullSIM Wide S-0","lep");
  leg->AddEntry(gr5,"2017 FullSIM narrow S-1","lep");
  leg->AddEntry(gr6,"2017 FullSIM wide S-1","lep");
  leg->Draw();
  c1->SetGrid();
  c1->Print("EffAcc17.png");
  c1->Print("EffAcc17.pdf");
  c1->Print("EffAcc17.svg");
  
  for(int i=0; i<massn17.size(); i++){
	  cout<<massn17.at(i)<<" "<<effn17.at(i)<<endl;
  }
  cout<<endl;
  for(int i=0; i<massw17.size(); i++){
	  cout<<massw17.at(i)<<" "<<effw17.at(i)<<endl;
  }
  cout<<endl;
  for(int i=0; i<massn17s1.size(); i++){
	  cout<<massn17s1.at(i)<<" "<<effn17s1.at(i)<<endl;
  }
  cout<<endl;
  for(int i=0; i<massw17s1.size(); i++){
	  cout<<massw17s1.at(i)<<" "<<effw17s1.at(i)<<endl;
  }
  
  //17 S-0 Acceptance and Efficiency
  float MCmassN[16]={700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3500,4000,5000};
  float AccN[16]={0.3521,0.35145,0.355,0.3514,0.3488,0.348,0.346,0.33425,0.3298,0.3312,0.3278,0.33145,0.32615,0.3229,0.304357,0.30488};
  float MCmassW[16]={700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3500,4000,5000};
  float AccW[16]={0.3367,0.3473,0.3427,0.3438,0.33825,0.33215,0.33145,0.3261,0.31745,0.3159,0.31725,0.3113,0.29985,0.2863,0.279128,0.255897};
  
  TGraphErrors *gr7 = new TGraphErrors(16,MCmassN,AccN,0,0);
  gr7->SetMarkerColor(2);
  gr7->SetMarkerStyle(20);
  gr7->SetMarkerSize(2.8);
  gr7->SetLineColor(2);
  gr7->SetLineStyle(2);
  gr7->SetLineWidth(2);
  TGraphErrors *gr8 = new TGraphErrors(16,MCmassW,AccW,0,0);
  gr8->SetMarkerColor(4);
  gr8->SetMarkerStyle(21);
  gr8->SetMarkerSize(2.8);
  gr8->SetLineColor(4);
  gr8->SetLineStyle(2);
  gr8->SetLineWidth(2);
  
  TCanvas *c2 = new TCanvas("c2","",2400,1800);
  c2->cd();
  c2->SetBottomMargin(0.11);
  c2->SetLeftMargin(0.13);
  gr7->SetTitle("Acceptance");
  xaxis = gr7->GetXaxis();
  yaxis = gr7->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  xaxis->SetLimits(700,5000); 
  yaxis->SetRangeUser(0,0.5);
  yaxis->SetTitle("Acceptance");
  yaxis->SetTitleOffset(1.2);
  gr7->Draw("AP");
  gr8->Draw("SAMEP");
  lumi_13TeV = "";
  CMS_lumi(c2,5,iPos);
  leg->Clear();
  leg->AddEntry(gr7,"2017 FullSIM Narrow S-0","lep");
  leg->AddEntry(gr8,"2017 FullSIM Wide S-0","lep");
  leg->Draw();
  c2->SetGrid();
  c2->Print("Acc17.png");
  c2->Print("Acc17.pdf");
  c2->Print("Acc17.root");
  c2->Print("Acc17.svg");
  
  std::vector<float> Eff_N, Eff_W;
  for(int i=0; i<16; i++){
	  Eff_N.push_back(effn17[i]/AccN[i]);
	  Eff_W.push_back(effw17[i]/AccW[i]);
  }
  
  TGraphErrors *gr9 = new TGraphErrors(16,MCmassN,&Eff_N[0],0,0);
  gr9->SetMarkerColor(2);
  gr9->SetMarkerStyle(20);
  gr9->SetMarkerSize(2.8);
  gr9->SetLineColor(2);
  gr9->SetLineStyle(2);
  gr9->SetLineWidth(2);
  TF1 *f9 = new TF1("fun9","pol4",700,5000);
  f9->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f9->SetLineColor(2);
  f9->SetLineStyle(2);
  gr9->Fit(f9,"R");
  gr9->Fit(f9,"R");
  f9->SetLineColor(2);
  f9->SetLineStyle(2);
  f9->SetLineWidth(2);
  TGraphErrors *gr10 = new TGraphErrors(16,MCmassW,&Eff_W[0],0,0);
  gr10->SetMarkerColor(4);
  gr10->SetMarkerStyle(20);
  gr10->SetMarkerSize(2.8);
  gr10->SetLineColor(4);
  gr10->SetLineStyle(2);
  gr10->SetLineWidth(2);
  TF1 *f10 = new TF1("fun10","pol4",700,5000);
  f10->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f10->SetLineColor(4);
  f10->SetLineStyle(2);
  gr10->Fit(f10,"R");
  gr10->Fit(f10,"R");
  f10->SetLineColor(4);
  f10->SetLineStyle(2);
  f10->SetLineWidth(2);

  
  TCanvas *c3 = new TCanvas("c3","",2400,1800);
  c3->cd();
  c3->SetBottomMargin(0.11);
  c3->SetLeftMargin(0.13);
  gr9->SetTitle("Acceptance");
  xaxis = gr9->GetXaxis();
  yaxis = gr9->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  xaxis->SetLimits(700,5000); 
  yaxis->SetRangeUser(0,0.5);
  yaxis->SetTitle("Efficiency");
  yaxis->SetTitleOffset(1.2);
  f9->Draw("APL");
  f10->Draw("APLSAME");
  gr9->Draw("AP");
  gr10->Draw("SAMEP");
  f9->Draw("APLSAME");
  f10->Draw("APLSAME");
  lumi_13TeV = "";
  CMS_lumi(c3,5,iPos);
  leg->Clear();
  leg->AddEntry(gr7,"2017 FullSIM Narrow S-0","lep");
  leg->AddEntry(gr8,"2017 FullSIM Wide S-0","lep");
  leg->Draw();
  c3->SetGrid();
  c3->Print("Eff17.png");
  c3->Print("Eff17.pdf");
  c3->Print("Eff17.svg");
  
  //17 S-1 Acceptance and Efficiency
  float MCmassNS1[7]={700,1200,2000,2800,3500,4000,5000};
  float AccNS1[7]={0.4075,0.4466,0.4518,0.4339,0.4276,0.421,0.408};
  float MCmassWS1[7]={700,1200,2000,2800,3500,4000,5000};
  float AccWS1[7]={0.3501,0.4159,0.4413,0.4419,0.4144,0.4033,0.4049};
  
  TGraphErrors *gr11 = new TGraphErrors(7,MCmassNS1,AccNS1,0,0);
  gr11->SetMarkerColor(2);
  gr11->SetMarkerStyle(21);
  gr11->SetMarkerSize(2.8);
  gr11->SetLineColor(2);
  gr11->SetLineStyle(2);
  gr11->SetLineWidth(2);
  TGraphErrors *gr12 = new TGraphErrors(7,MCmassWS1,AccWS1,0,0);
  gr12->SetMarkerColor(4);
  gr12->SetMarkerStyle(21);
  gr12->SetMarkerSize(2.8);
  gr12->SetLineColor(4);
  gr12->SetLineStyle(2);
  gr12->SetLineWidth(2);
  
  TCanvas *c4 = new TCanvas("c4","",2400,1800);
  c4->cd();
  c4->SetBottomMargin(0.11);
  c4->SetLeftMargin(0.13);
  gr11->SetTitle("Acceptance");
  xaxis = gr11->GetXaxis();
  yaxis = gr11->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  xaxis->SetLimits(700,5000); 
  yaxis->SetRangeUser(0,0.7);
  yaxis->SetTitle("Acceptance");
  yaxis->SetTitleOffset(1.2);
  gr11->Draw("AP");
  gr12->Draw("SAMEP");
  lumi_13TeV = "";
  CMS_lumi(c4,5,iPos);
  leg->Clear();
  leg->AddEntry(gr11,"2017 FullSIM Narrow S-1","lep");
  leg->AddEntry(gr12,"2017 FullSIM Wide S-1","lep");
  leg->Draw();
  c4->SetGrid();
  c4->Print("Acc17S1.png");
  c4->Print("Acc17S1.pdf");
  c4->Print("Acc17S1.root");
  c4->Print("Acc17S1.svg");
  
  std::vector<float> Eff_NS1, Eff_WS1;
  for(int i=0; i<7; i++){
	  Eff_NS1.push_back(effn17s1[i]/AccNS1[i]);
	  Eff_WS1.push_back(effw17s1[i]/AccWS1[i]);
  }
  
  TGraphErrors *gr13 = new TGraphErrors(7,MCmassNS1,&Eff_NS1[0],0,0);
  gr13->SetMarkerColor(2);
  gr13->SetMarkerSize(2.8);
  gr13->SetLineColor(2);
  gr13->SetLineStyle(2);
  gr13->SetLineWidth(2);
  TGraphErrors *gr14 = new TGraphErrors(7,MCmassWS1,&Eff_WS1[0],0,0);
  gr14->SetMarkerColor(4);
  gr14->SetMarkerStyle(21);
  gr14->SetMarkerSize(2.8);
  gr14->SetLineColor(4);
  gr14->SetLineStyle(2);
  gr14->SetLineWidth(2);
  TCanvas *c5 = new TCanvas("c5","",2400,1800);
  c5->cd();
  c5->SetBottomMargin(0.11);
  c5->SetLeftMargin(0.13);
  gr14->SetTitle("Acceptance");
  xaxis = gr14->GetXaxis();
  yaxis = gr14->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  xaxis->SetLimits(700,5000); 
  yaxis->SetRangeUser(0,0.5);
  yaxis->SetTitle("Efficiency");
  yaxis->SetTitleOffset(1.2);
  gr14->Draw("AP");
  gr14->Draw("SAMEP");
  lumi_13TeV = "";
  CMS_lumi(c5,5,iPos);
  leg->Clear();
  leg->AddEntry(gr11,"2017 FullSIM Narrow S-1","lep");
  leg->AddEntry(gr14,"2017 FullSIM Wide S-1","lep");
  leg->Draw();
  c5->SetGrid();
  c5->Print("Eff17S1.png");
  c5->Print("Eff17S1.pdf");
  c5->Print("Eff17S1.svg");
}