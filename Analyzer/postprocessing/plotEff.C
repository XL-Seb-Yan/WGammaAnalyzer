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
  
  ifstream file16n("log16N.txt");
  ifstream file17n("log17N.txt");
  ifstream file16w("log16W.txt");
  ifstream file17w("log17W.txt");
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
			  effn17.push_back(strtof((substr).c_str(),0) / 20000);
			  effn_err17.push_back(sqrt(strtof((substr).c_str(),0)) / 20000);
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
			  effw17.push_back(strtof((substr).c_str(),0) / 20000);
			  effw_err17.push_back(sqrt(strtof((substr).c_str(),0)) / 20000);
		  }
       }
   }
  }
  
  // READ SIGNAL NORMALIZATION FOR GOF TEST
    // while (getline(file17n,str)) {
   // if(str.find("Normalization:") != std::string::npos){
	   // std::stringstream ss(str);
	   // int index = -1;
	   // while(ss.good()){
          // string substr;
          // getline(ss,substr,' ');
		  // index++;
          // if(index == 0) continue;
		  // if(index == 2) massn17.push_back(strtof((substr).c_str(),0));
		  // if(index == 1){
			  // effn17.push_back(strtof((substr).c_str(),0) / 20000);
			  // effn_err17.push_back(sqrt(strtof((substr).c_str(),0)) / 20000);
		  // }
       // }
   // }
  // }
  // while (getline(file17w,str)) {
   // if(str.find("Normalization:") != std::string::npos){
	   // std::stringstream ss(str);
	   // int index = -1;
	   // while(ss.good()){
          // string substr;
          // getline(ss,substr,' ');
		  // index++;
          // if(index == 0) continue;
		  // if(index == 2) massw17.push_back(strtof((substr).c_str(),0));
		  // if(index == 1){
			  // effw17.push_back(strtof((substr).c_str(),0) / 20000);
			  // effw_err17.push_back(sqrt(strtof((substr).c_str(),0)) / 20000);
		  // }
       // }
   // }
  // }

  //17 spin1
  double narrow[5] = {0.160777,0.175728,0.174048,0.1568783,0.149526316};
  double narrow_err[5] = {0.004014262,0.004202095,0.004184388,0.003978605,0.003967318};
  double wide[5] = {0.146991,0.165444,0.161994,0.150717,0.144221};
  double wide_err[5] = {0.00404133,0.004287507,0.0040248478,0.0038822287,0.003797644};
  double massn[5] = {700,1200,2000,2800,3500};
  double massw[5] = {700,1200,2000,2800,3500};
  
  
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
  TF1 *f3 = new TF1("fun3","pol4",700,3500);
  f3->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
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
  TF1 *f4 = new TF1("fun4","pol4",700,3500);
  f4->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  f4->SetLineColor(4);
  gr4->Fit(f4,"R");
  gr4->Fit(f4,"R");
  f4->SetLineColor(4);
  f4->SetLineWidth(2);
  
  cout<<"=========================Fitting 17 S1==========================="<<endl;
  TGraphErrors *gr5 = new TGraphErrors(5,&massn[0],&narrow[0],0,&narrow_err[0]);
  gr5->SetMarkerColor(2);
  gr5->SetMarkerStyle(21);
  //gr1->SetMarkerStyle(20);
  gr5->SetMarkerSize(2.8);
  gr5->SetLineWidth(2);
  gr5->SetLineColor(2);
  gr5->SetLineStyle(2);
  TF1 *f5 = new TF1("fun5","pol4",700,3500);
  f5->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),-6.1749*pow(10,-15));
  f5->SetLineColor(2);
  f5->SetLineStyle(2);
  gr5->Fit(f5,"R");
  gr5->Fit(f5,"R");
  f5->SetLineColor(2);
  f5->SetLineStyle(2);
  f5->SetLineWidth(2);

  TGraphErrors *gr6 = new TGraphErrors(5,&massw[0],&wide[0],0,&wide_err[0]);
  gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(21);
  //gr2->SetMarkerStyle(20);
  gr6->SetMarkerSize(2.8);
  gr6->SetLineColor(4);
  gr6->SetLineStyle(2);
  gr6->SetLineWidth(2);
  TF1 *f6 = new TF1("fun6","pol4",700,3500);
  f6->SetParameters(0.122574,0.000114036,-9.69702*pow(10,-8),2.80535*pow(10,-11),-2.86354*pow(10,-15));
  f6->SetLineColor(4);
  f6->SetLineStyle(2);
  gr6->Fit(f6,"R");
  gr6->Fit(f6,"R");
  f6->SetLineColor(4);
  f6->SetLineStyle(2);
  f6->SetLineWidth(2);

  
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

  TLegend *leg = new TLegend(0.6,0.8,0.9,0.9);
  
  TCanvas *c0 = new TCanvas("c0","",2400,1800);
  c0->cd();
  c0->SetBottomMargin(0.11);
  c0->SetLeftMargin(0.13);
  gr1->SetTitle("Total Efficiency");
  TAxis *xaxis = gr1->GetXaxis();
  TAxis *yaxis = gr1->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  xaxis->SetLimits(700,3500); 
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
  c0->Print("Eff16.png");
  c0->Print("Eff16.pdf");
  c0->Print("Eff16.root");
  c0->Print("Eff16.svg");
  
  leg = new TLegend(0.6,0.7,0.9,0.9);
  TCanvas *c1 = new TCanvas("c1","",2400,1800);
  c1->cd();
  c1->SetBottomMargin(0.11);
  c1->SetLeftMargin(0.13);
  gr3->SetTitle("Total Efficiency");
  xaxis = gr3->GetXaxis();
  yaxis = gr3->GetYaxis();
  xaxis->SetTitle("m_{X} (GeV)");
  xaxis->SetLimits(700,3500); 
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
  c1->Print("Eff17.png");
  c1->Print("Eff17.pdf");
  c1->Print("Eff17.root");
  c1->Print("Eff17.svg");
  
  cout<<massn17.size()<<endl;
  for(int i=0; i<massn17.size(); i++){
	  cout<<massn17.at(i)<<" "<<effn17.at(i)<<endl;
  }
  for(int i=0; i<massw17.size(); i++){
	  cout<<massw17.at(i)<<" "<<effw17.at(i)<<endl;
  }
}
