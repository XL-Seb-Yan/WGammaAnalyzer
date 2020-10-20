#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
void plotEff(){
  gROOT->SetBatch(1);
  lumi_13TeV = "";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Simulation";
  lumiTextSize = 0.45;
  cmsTextSize = 0.6;
  int iPeriod = 5;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetHistLineWidth(3);
  gStyle->SetEndErrorSize(5);
  
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
  
  ifstream file16n("EffS0/log16N_eff.txt");
  ifstream file17n("EffS0/log17N_eff_new.txt");
  ifstream file16w("EffS0/log16W_eff.txt");
  ifstream file17w("EffS0/log17W_eff_new.txt");
  ifstream file17ns1("EffS1/log17NS1_newmpoints.txt");
  ifstream file17ws1("EffS1/log17WS1_newmpoints.txt");
  string str; 
  
  
  // while (getline(file16n,str)) {
   // if(str.find("++++") != std::string::npos){
	   // std::stringstream ss(str);
	   // int index = -1;
	   // while(ss.good()){
          // string substr;
          // getline(ss,substr,',');
		  // index++;
          // if(index == 0) continue;
		  // if(index == 1) massn16.push_back(strtof((substr).c_str(),0));
		  // if(index == 2){
			  // effn16.push_back(strtof((substr).c_str(),0) / 20000);
			  // effn_err16.push_back(sqrt(strtof((substr).c_str(),0)) / 20000);
		  // }
       // }
   // }
  // }
  // while (getline(file16w,str)) {
   // if(str.find("++++") != std::string::npos){
	   // std::stringstream ss(str);
	   // int index = -1;
	   // while(ss.good()){
          // string substr;
          // getline(ss,substr,',');
		  // index++;
          // if(index == 0) continue;
		  // if(index == 1) massw16.push_back(strtof((substr).c_str(),0));
		  // if(index == 2){
			  // effw16.push_back(strtof((substr).c_str(),0) / 20000);
			  // effw_err16.push_back(sqrt(strtof((substr).c_str(),0)) / 20000);
		  // }
       // }
   // }
  // }
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
			  effn_err17.push_back(sqrt(strtof((substr).c_str(),0)*(1-strtof((substr).c_str(),0) / Ngen)) / Ngen);
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
			  effw_err17.push_back(sqrt(strtof((substr).c_str(),0)*(1-strtof((substr).c_str(),0) / Ngen))/ Ngen);
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
			  int Ngen = 20000;
              if(abs(massn17s1.back()-1500) < 1) //M-4000 only has 7000 events
                Ngen = 13000;
              if(abs(massn17s1.back()-2500) < 1) //M-4000 only has 7000 events
                Ngen = 19000;
              if(abs(massn17s1.back()-6000) < 1) //M-4000 only has 7000 events
                Ngen = 19000;
			  effn17s1.push_back(strtof((substr).c_str(),0) / Ngen);
			  //effn_err17s1.push_back(sqrt(strtof((substr).c_str(),0)) / Ngen);
              effn_err17s1.push_back(sqrt(strtof((substr).c_str(),0)*(1-strtof((substr).c_str(),0) / Ngen))/ Ngen);
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
			  int Ngen = 20000;
              if(abs(massw17s1.back()-700) < 1) //M-4000 only has 7000 events
                Ngen = 19000;
              if(abs(massw17s1.back()-1500) < 1) //M-4000 only has 7000 events
                Ngen = 15000;
              if(abs(massw17s1.back()-2000) < 1) //M-4000 only has 7000 events
                Ngen = 16000;
              if(abs(massw17s1.back()-4000) < 1) //M-4000 only has 7000 events
                Ngen = 18000;
              if(abs(massw17s1.back()-6000) < 1) //M-4000 only has 7000 events
                Ngen = 17000;
              if(abs(massw17s1.back()-8000) < 1) //M-4000 only has 7000 events
                Ngen = 19000;
			  effw17s1.push_back(strtof((substr).c_str(),0) / Ngen);
			  effw_err17s1.push_back(sqrt(strtof((substr).c_str(),0)*(1-strtof((substr).c_str(),0) / Ngen))/ Ngen);
		  }
       }
   }
  }
  // cout<<"=========================Fitting 16==========================="<<endl;
  // TGraphErrors *gr1 = new TGraphErrors(massn16.size(),&massn16[0],&effn16[0],0,&effn_err16[0]);
  // gr1->SetMarkerColor(2);
  // gr1->SetMarkerStyle(20);
  // gr1->SetMarkerSize(2.8);
  // gr1->SetLineWidth(2);
  // gr1->SetLineColor(2);
  // TF1 *f1 = new TF1("fun1","pol4",700,3500);
  // f1->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  // f1->SetLineColor(2);
  // gr1->Fit(f1,"R");
  // gr1->Fit(f1,"R");
  // f1->SetLineColor(2);
  // f1->SetLineWidth(2);
  
  // TGraphErrors *gr2 = new TGraphErrors(massw16.size(),&massw16[0],&effw16[0],0,&effw_err16[0]);
  // gr2->SetMarkerColor(4);
  // gr2->SetMarkerStyle(20);
  // gr2->SetMarkerSize(2.8);
  // gr2->SetLineWidth(2);
  // gr2->SetLineColor(4);
  // TF1 *f2 = new TF1("fun2","pol4",700,3500);
  // f2->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  // f2->SetLineColor(4);
  // gr2->Fit(f2,"R");
  // gr2->Fit(f2,"R");
  // f2->SetLineColor(4);
  // f2->SetLineWidth(2);
  
  cout<<"=========================Fitting 17==========================="<<endl;
  TGraphErrors *gr3 = new TGraphErrors(massn17.size(),&massn17[0],&effn17[0],0,&effn_err17[0]);
  gr3->SetMarkerColor(2);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(2);
  gr3->SetLineWidth(2);
  gr3->SetLineColor(2);
  TF1 *f3 = new TF1("fun3","pol4",700,8000);
  f3->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),0,0);
  f3->SetLineColor(2);
  gr3->Fit(f3,"R");
  gr3->Fit(f3,"R");
  f3->SetLineColor(2);
  f3->SetLineWidth(2);
  cout<<f3->GetParameter(0)<<","<<f3->GetParameter(1)<<","<<f3->GetParameter(2)<<","<<f3->GetParameter(3)<<","<<f3->GetParameter(4)<<endl;
  
  TGraphErrors *gr4 = new TGraphErrors(massw17.size(),&massw17[0],&effw17[0],0,&effw_err17[0]);
  gr4->SetMarkerColor(4);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(2);
  gr4->SetLineWidth(2);
  gr4->SetLineColor(4);
  TF1 *f4 = new TF1("fun3","pol4",700,8000);
  f4->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),0,0);
  f4->SetLineColor(4);
  gr4->Fit(f4,"R");
  gr4->Fit(f4,"R");
  f4->SetLineColor(4);
  f4->SetLineWidth(2);
  cout<<f4->GetParameter(0)<<","<<f4->GetParameter(1)<<","<<f4->GetParameter(2)<<","<<f4->GetParameter(3)<<","<<f4->GetParameter(4)<<endl;
  
  cout<<endl;
  cout<<endl;
  cout<<"=========================Fitting 17 S1==========================="<<endl;
  TGraphErrors *gr5 = new TGraphErrors(8,&massn17s1[0],&effn17s1[0],0,&effn_err17s1[0]);
  gr5->SetMarkerColor(2);
  gr5->SetMarkerStyle(21);
  gr5->SetMarkerSize(2);
  gr5->SetLineWidth(2);
  gr5->SetLineColor(2);
  TF1 *f5 = new TF1("fun5","pol4",700,8000);
  f5->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11));
  f5->SetLineColor(2);
  f5->SetLineStyle(7);
  gr5->Fit(f5,"R");
  gr5->Fit(f5,"R");
  f5->SetLineColor(2);
  f5->SetLineStyle(2);
  f5->SetLineWidth(2);
  cout<<f5->GetParameter(0)<<","<<f5->GetParameter(1)<<","<<f5->GetParameter(2)<<","<<f5->GetParameter(3)<<","<<f5->GetParameter(4)<<endl;

  TGraphErrors *gr6 = new TGraphErrors(8,&massw17s1[0],&effw17s1[0],0,&effw_err17s1[0]);
  gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(21);
  gr6->SetMarkerSize(2);
  gr6->SetLineColor(4);
  gr6->SetLineWidth(2);
  TF1 *f6 = new TF1("fun6","pol4",700,8000);
  f6->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f6->SetLineColor(4);
  f6->SetLineStyle(2);
  gr6->Fit(f6,"R");
  gr6->Fit(f6,"R");
  f6->SetLineColor(4);
  f6->SetLineStyle(2);
  f6->SetLineWidth(2);
  cout<<f6->GetParameter(0)<<","<<f6->GetParameter(1)<<","<<f6->GetParameter(2)<<","<<f6->GetParameter(3)<<","<<f6->GetParameter(4)<<endl;

  TLegend *leg = new TLegend(0.6,0.8,0.9,0.9);
  TAxis *xaxis = NULL;
  TAxis *yaxis = NULL;
  
  // TCanvas *c0 = new TCanvas("c0","",2400,1800);
  // c0->cd();
  // c0->SetBottomMargin(0.11);
  // c0->SetLeftMargin(0.13);
  // gr1->SetTitle("Total Efficiency");
  // TAxis *xaxis = gr1->GetXaxis();
  // TAxis *yaxis = gr1->GetYaxis();
  // xaxis->SetTitle("m_{X} [GeV]");
  // xaxis->SetLimits(700,5000); 
  // yaxis->SetRangeUser(0,0.3);
  // yaxis->SetTitle("Acc. #times Eff.");
  // yaxis->SetTitleOffset(1.2);
  // f1->Draw("APL");
  // f1->Draw("APL");
  // gr1->Draw("AP");
  // gr2->Draw("SAMEP");
  // f1->Draw("PLSAME");
  // f2->Draw("PLSAME");
  // lumi_13TeV = "35.92 fb^{-1}";
  // CMS_lumi(c0,4,iPos);
  // leg->AddEntry(gr1,"2016 FullSIM Narrow S-0","lep");
  // leg->AddEntry(gr2,"2016 FullSIM Broad S-0","lep");
  // leg->Draw();
  // c0->SetGrid();
  // c0->Print("EffAcc16.png");
  // c0->Print("EffAcc16.pdf");
  // c0->Print("EffAcc16.svg");
  
  leg = new TLegend(0.5,0.65,0.85,0.75);
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
  xaxis->SetLimits(700,8000); 
  yaxis->SetRangeUser(0,0.3);
  yaxis->SetTitle("A #varepsilon");
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
  lumi_13TeV = "";
  CMS_lumi(c1,5,iPos);
  leg->Clear();
  leg->AddEntry(gr3,"Narrow S-0","lep");
  leg->AddEntry(gr4,"Broad S-0","lep");
  leg->AddEntry(gr5,"Narrow S-1","lep");
  leg->AddEntry(gr6,"Broad S-1","lep");
  leg->Draw();
  c1->SetGrid();
  c1->Print("EffAcc17.pdf");
  c1->Print("EffAcc17.png");
  c1->Print("EffAcc17.svg");
  
  // for(int i=0; i<massn17.size(); i++){
	  // cout<<massn17.at(i)<<" "<<effn17.at(i)<<endl;
  // }
  // cout<<endl;
  // for(int i=0; i<massw17.size(); i++){
	  // cout<<massw17.at(i)<<" "<<effw17.at(i)<<endl;
  // }
  // cout<<endl;
  // for(int i=0; i<massn17s1.size(); i++){
	  // cout<<massn17s1.at(i)<<" "<<effn17s1.at(i)<<endl;
  // }
  // cout<<endl;
  // for(int i=0; i<massw17s1.size(); i++){
	  // cout<<massw17s1.at(i)<<" "<<effw17s1.at(i)<<endl;
  // }
  
  //17 S-0 Acceptance
  float MCmassN[19]={700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3500,4000,5000,6000,7000,8000};
  float AccN[19]={0.3521,0.35145,0.355,0.3514,0.3488,0.348,0.346,0.33425,0.3298,0.3312,0.3278,0.33145,0.32615,0.3229,0.314357,0.30488,0.3065,0.2919,0.288};
  float AccN_err[19];
  for(int i=0; i<19; i++){
      int Ngen = 20000;
      if(abs(MCmassN[i]-4000) < 1) //M-4000 only has 19000 events
          Ngen = 14000;
      if(abs(MCmassN[i]-5000) < 1) //M-4000 only has 19000 events
          Ngen = 12500;
      AccN_err[i]=sqrt(AccN[i]*(1-AccN[i])/Ngen);
  }
  float MCmassW[19]={700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3500,4000,5000,6000,7000,8000};
  float AccW[19]={0.3367,0.3473,0.3427,0.3438,0.33825,0.33215,0.33145,0.3261,0.31745,0.3159,0.31725,0.3113,0.29985,0.2863,0.279128,0.255897,0.229,0.1959,0.1427};
  float AccW_err[19];
  for(int i=0; i<19; i++){
      int Ngen = 20000;
      if(abs(MCmassW[i]-4000) < 1) //M-4000 only has 19000 events
          Ngen = 19500;
      if(abs(MCmassW[i]-5000) < 1) //M-4000 only has 19000 events
          Ngen = 19500;
      AccW_err[i]=sqrt(AccW[i]*(1-AccW[i])/Ngen);
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
  TGraphErrors *gr7 = new TGraphErrors(19,MCmassN,AccN,0,AccN_err);
  gr7->SetMarkerColor(2);
  gr7->SetMarkerStyle(20);
  gr7->SetMarkerSize(2);
  gr7->SetLineColor(2);
  gr7->SetLineWidth(2);
  TF1 *f7 = new TF1("fun7","pol4",700,8000);
  f7->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f7->SetLineColor(2);
  gr7->Fit(f7,"R");
  gr7->Fit(f7,"R");
  f7->SetLineColor(2);
  f7->SetLineWidth(2);
  TGraphErrors *gr8 = new TGraphErrors(19,MCmassW,AccW,0,AccW_err);
  gr8->SetMarkerColor(4);
  gr8->SetMarkerStyle(20);
  gr8->SetMarkerSize(2);
  gr8->SetLineColor(4);
  gr8->SetLineWidth(2);
  TF1 *f8 = new TF1("fun8","pol4",700,8000);
  f8->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f8->SetLineColor(4);
  gr8->Fit(f8,"R");
  gr8->Fit(f8,"R");
  f8->SetLineColor(4);
  f8->SetLineWidth(2);
  
  cout<<endl;
  cout<<endl;
  cout<<"=========================Fitting 17 S1 Acc ==========================="<<endl;
  TGraphErrors *gr11 = new TGraphErrors(8,MCmassNS1,AccNS1,0,AccNS1_err);
  gr11->SetMarkerColor(2);
  gr11->SetMarkerStyle(21);
  gr11->SetMarkerSize(2);
  gr11->SetLineColor(2);
  gr11->SetLineWidth(2);
  TF1 *f11 = new TF1("fun11","pol4",700,8000);
  f11->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f11->SetLineColor(2);
  f11->SetLineStyle(2);
  gr11->Fit(f11,"R");
  gr11->Fit(f11,"R");
  f11->SetLineColor(2);
  f11->SetLineStyle(2);
  f11->SetLineWidth(2);
  TGraphErrors *gr12 = new TGraphErrors(8,MCmassWS1,AccWS1,0,AccWS1_err);
  gr12->SetMarkerColor(4);
  gr12->SetMarkerStyle(21);
  gr12->SetMarkerSize(2);
  gr12->SetLineColor(4);
  gr12->SetLineWidth(2);
  TF1 *f12 = new TF1("fun12","pol4",700,8000);
  f12->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11));
  f12->SetLineColor(4);
  f12->SetLineStyle(2);
  gr12->Fit(f12,"R");
  gr12->Fit(f12,"R");
  f12->SetLineColor(4);
  f12->SetLineStyle(2);
  f12->SetLineWidth(2);
  
  TCanvas *c2 = new TCanvas("c2","",2400,1800);
  c2->cd();
  c2->SetBottomMargin(0.11);
  c2->SetLeftMargin(0.13);
  gr7->SetTitle("A");
  xaxis = gr7->GetXaxis();
  yaxis = gr7->GetYaxis();
  xaxis->SetTitle("m_{X} [GeV]");
  xaxis->SetLimits(700,8000); 
  yaxis->SetRangeUser(0,0.7);
  yaxis->SetTitle("A");
  yaxis->SetTitleOffset(1.2);
  f7->Draw("APL");
  f8->Draw("APLSAME");
  f11->Draw("APLSAME");
  f12->Draw("APLSAME");
  gr7->Draw("AP");
  gr8->Draw("SAMEP");
  gr11->Draw("SAMEP");
  gr12->Draw("SAMEP");
  f7->Draw("APLSAME");
  f8->Draw("APLSAME");
  f11->Draw("APLSAME");
  f12->Draw("APLSAME");
  lumi_13TeV = "";
  CMS_lumi(c2,5,iPos);
  leg->Clear();
  leg->AddEntry(gr7,"Narrow S-0","lep");
  leg->AddEntry(gr8,"Broad S-0","lep");
  leg->AddEntry(gr11,"Narrow S-1","lep");
  leg->AddEntry(gr12,"Broad S-1","lep");
  leg->Draw();
  c2->SetGrid();
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
  xaxis->SetLimits(700,8000); 
  yaxis->SetRangeUser(0.0,0.6);
  yaxis->SetTitle("#varepsilon");
  yaxis->SetTitleOffset(1.2);
  gr9->SetLineColor(2);
  gr10->SetLineColor(4);
  gr13->SetLineColor(2);
  gr14->SetLineColor(4);
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
  c3->SetGrid();
  c3->Print("Eff17.pdf");
  c3->Print("Eff17.svg");
  
  /*
  std::vector<float> Eff_N, Eff_W;
  for(int i=0; i<effn17.size(); i++){
	  Eff_N.push_back(effn17[i]/AccN[i]);
	  Eff_W.push_back(effw17[i]/AccW[i]);
  }
  
  cout<<endl;
  cout<<endl;
  cout<<"=========================Fitting 17 n Eff ==========================="<<endl;
  TGraphErrors *gr9 = new TGraphErrors(19,MCmassN,&Eff_N[0],0,0);
  gr9->SetMarkerColor(2);
  gr9->SetMarkerStyle(20);
  gr9->SetMarkerSize(2);
  gr9->SetLineColor(2);
  gr9->SetLineWidth(2);
  TF1 *f9 = new TF1("fun9","[0]*(1.-TMath::Erf((x-[1])/[2]))+[3]",700,8000);
  f9->SetParameters(0.05,3000,100,0.36);
  f9->SetLineColor(2);
  gr9->Fit(f9,"R");
  gr9->Fit(f9,"R");
  f9->SetLineColor(2);
  f9->SetLineWidth(2);
  cout<<"=========================Fitting 17 w Eff ==========================="<<endl;
  TGraphErrors *gr10 = new TGraphErrors(19,MCmassW,&Eff_W[0],0,0);
  gr10->SetMarkerColor(4);
  gr10->SetMarkerStyle(20);
  gr10->SetMarkerSize(2);
  gr10->SetLineColor(4);
  gr10->SetLineWidth(2);
  TF1 *f10 = new TF1("fun10","[0]*(1.-TMath::Erf((x-[1])/[2]))+[3]",700,8000);
  f10->SetParameters(0.02,2620,1200,0.32);
  f10->SetParLimits(2,600,1000);
  f10->SetLineColor(4);
  gr10->Fit(f10,"R");
  gr10->Fit(f10,"R");
  f10->SetLineColor(4);
  f10->SetLineWidth(2);
  
  std::vector<float> Eff_NS1, Eff_WS1;
  for(int i=0; i<8; i++){
	  Eff_NS1.push_back(effn17s1[i]/AccNS1[i]);
	  Eff_WS1.push_back(effw17s1[i]/AccWS1[i]);
  }
  
  TGraphErrors *gr13 = new TGraphErrors(8,MCmassNS1,&Eff_NS1[0],0,0);
  gr13->SetMarkerColor(2);
  gr13->SetMarkerStyle(21);
  gr13->SetMarkerSize(2);
  gr13->SetLineColor(2);
  gr13->SetLineWidth(2);
  TF1 *f13 = new TF1("fun13","[0]*(1.-TMath::Erf((x-[1])/[2]))+[3]",700,8000);
  f13->SetParameters(0.02,2620,1200,0.32);
  f13->SetParLimits(2,600,1000);
  f13->SetLineColor(2);
  f13->SetLineStyle(2);
  gr13->Fit(f13,"R");
  gr13->Fit(f13,"R");
  f13->SetLineColor(2);
  f13->SetLineStyle(2);
  f13->SetLineWidth(2);

  TGraphErrors *gr14 = new TGraphErrors(8,MCmassWS1,&Eff_WS1[0],0,0);
  gr14->SetMarkerColor(4);
  gr14->SetMarkerStyle(21);
  gr14->SetMarkerSize(2);
  gr14->SetLineColor(4);
  gr14->SetLineWidth(2);
  TF1 *f14 = new TF1("fun14","[0]*(1.-TMath::Erf((x-[1])/[2]))+[3]",700,8000);
  f14->SetParameters(0.02,2620,1200,0.32);
  f14->SetParLimits(2,600,1000);
  f14->SetLineColor(4);
  f14->SetLineStyle(2);
  gr14->Fit(f14,"R");
  gr14->Fit(f14,"R");
  f14->SetLineColor(4);
  f14->SetLineStyle(2);
  f14->SetLineWidth(2);

  
  TCanvas *c3 = new TCanvas("c3","",2400,1800);
  c3->cd();
  c3->SetBottomMargin(0.11);
  c3->SetLeftMargin(0.13);
  gr9->SetTitle("Acceptance");
  xaxis = gr9->GetXaxis();
  yaxis = gr9->GetYaxis();
  xaxis->SetTitle("m_{X} [GeV]");
  xaxis->SetLimits(700,8000); 
  yaxis->SetRangeUser(0.2,0.5);
  yaxis->SetTitle("Efficiency");
  yaxis->SetTitleOffset(1.2);
  f9->Draw("APL");
  f10->Draw("APLSAME");
  gr9->Draw("AP");
  gr10->Draw("SAMEP");
  f9->Draw("APLSAME");
  f10->Draw("APLSAME");
  f13->Draw("APLSAME");
  f14->Draw("APLSAME");
  gr13->Draw("SAMEP");
  gr14->Draw("SAMEP");
  f13->Draw("APLSAME");
  f14->Draw("APLSAME");
  lumi_13TeV = "";
  CMS_lumi(c3,5,iPos);
  leg->Clear();
  leg->AddEntry(gr9,"FullSIM Narrow S-0","lep");
  leg->AddEntry(gr10,"FullSIM Broad S-0","lep");
  leg->AddEntry(gr13,"FullSIM Narrow S-1","lep");
  leg->AddEntry(gr14,"FullSIM Broad S-1","lep");
  leg->Draw();
  c3->SetGrid();
  c3->Print("Eff17.png");
  c3->Print("Eff17.pdf");
  c3->Print("Eff17.svg");
  */
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
  
}