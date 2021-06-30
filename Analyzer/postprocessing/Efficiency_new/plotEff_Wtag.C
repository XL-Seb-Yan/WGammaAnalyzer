#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
void plotEff_Wtag(){
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
  
  //These efficiencies are correct to reshow a similar behavior of reconstruction efficiencies
  float MCmassS0[20]={700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500,4000,5000,6000,7000,8000};
  float EffN[20]={0.6325,0.5927,0.5674,0.5764,0.5617,0.5595,0.5329,0.5202,0.5170,0.4916,0.4875,0.4872,0.4910,0.4888,0.4762,0.4736,0.4627,0.4728,0.4385,0.3987};
  for(int i=0; i<20; i++){
      massn17.push_back(MCmassS0[i]);
      effn17.push_back(EffN[i]);
      effn_err17.push_back(sqrt(EffN[i]*(1-EffN[i])/20000));
  }
  float EffW[20]={0.6275,0.6173,0.5820,0.5785,0.5643,0.5503,0.5321,0.5354,0.5060,0.4961,0.4939,0.5009,0.4981,0.4905,0.4772,0.4767,0.4823,0.4451,0.4239,0.3988};
  for(int i=0; i<20; i++){
      massw17.push_back(MCmassS0[i]);
      effw17.push_back(EffW[i]);
      effn_err17.push_back(sqrt(EffW[i]*(1-EffW[i])/20000));
  }
  
  //These efficiencies are correct to reshow a similar behavior of reconstruction efficiencies
  float MCmassS1[8]={700,1500,2000,2500,3000,4000,6000,8000};
  float EffNS1[8]={0.6209,0.5496,0.5269,0.4940,0.4701,0.4858,0.4601,0.4015};
  for(int i=0; i<8; i++){
      massn17s1.push_back(MCmassS1[i]);
      effn17s1.push_back(EffNS1[i]);
      effn_err17s1.push_back(sqrt(EffNS1[i]*(1-EffNS1[i])/10000));
  }
  float EffWS1[8]={0.6002,0.5460,0.5158,0.4862,0.4758,0.4826,0.4693,0.4136};
  for(int i=0; i<8; i++){
      massw17s1.push_back(MCmassS1[i]);
      effw17s1.push_back(EffWS1[i]);
      effw_err17s1.push_back(sqrt(EffWS1[i]*(1-EffWS1[i])/10000));
  }
  
  cout<<"=========================Fitting 17==========================="<<endl;
  TGraphErrors *gr3 = new TGraphErrors(massn17.size(),&massn17[0],&effn17[0],0,&effn_err17[0]);
  gr3->SetMarkerColor(2);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(2);
  gr3->SetLineWidth(2);
  gr3->SetLineColor(kRed-4);
  TF1 *f3 = new TF1("fun3","pol4",700,8000);
  f3->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11),0,0);
  f3->SetLineColor(kRed-4);
  gr3->Fit(f3,"R");
  gr3->Fit(f3,"R");
  f3->SetLineColor(kRed-4);
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
  gr5->SetLineColor(kRed-4);
  TF1 *f5 = new TF1("fun5","pol4",700,8000);
  f5->SetParameters(0.0796303,0.000227561,-1.8753*pow(10,-7),5.7758*pow(10,-11));
  f5->SetLineColor(kRed-4);
  f5->SetLineStyle(7);
  gr5->Fit(f5,"R");
  gr5->Fit(f5,"R");
  f5->SetLineColor(kRed-4);
  f5->SetLineStyle(2);
  f5->SetLineWidth(2);
  cout<<f5->GetParameter(0)<<","<<f5->GetParameter(1)<<","<<f5->GetParameter(2)<<","<<f5->GetParameter(3)<<","<<f5->GetParameter(4)<<endl;

  for(int i=0; i<8; i++){
    cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
    cout<<effw17s1[i]<<endl;
  }
  
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
  
  TF1 effacc_S0N("effacc","pol4",500,8200);
  effacc_S0N.SetParameters(0.0743178,5.97851e-05,-2.35105e-08,3.22656e-12,-1.55813e-16);
  
  TF1 acc_S0N("acc","pol4",500,8200);
  acc_S0N.SetParameters(0.203957,0.000159745,-6.38119e-08,9.10934e-12,-4.58625e-16);
  
  TGraphErrors *gr1 = new TGraphErrors();
  int ipoints = 0;
	for(int i=700; i<8001; i+=10){
    TH1F temp("","",1000,0,1);
    temp.Fill(f3->Eval(i));
	  temp.Fill(f4->Eval(i));
	  temp.Fill(f5->Eval(i));
	  temp.Fill(f6->Eval(i));
	  gr1->SetPoint(ipoints,i,((f3->Eval(i) + f4->Eval(i) + f5->Eval(i) + f6->Eval(i))/4));
	  gr1->SetPointError(ipoints,0,temp.GetStdDev());
	  ipoints++;
    //cout<<i<<" "<<137.19 * effacc_S0N.Eval(i) / acc_S0N.Eval(i) / ((f3->Eval(i) + f4->Eval(i) + f5->Eval(i) + f6->Eval(i))/4) <<" :"<<endl;
    cout<<i<<" "<<((f3->Eval(i) + f4->Eval(i) + f5->Eval(i) + f6->Eval(i))/4)<<" "<<temp.GetStdDev()<<endl;
	}
  
  gr1->SetLineColor(6);
	gr1->SetFillColor(6);
	gr1->SetFillStyle(3002);
  gr1->SetLineWidth(3);
  gr1->SetFillColorAlpha(6,0.5);

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
  xaxis = gr1->GetXaxis();
  yaxis = gr1->GetYaxis();
  xaxis->SetTitle("m_{X} [GeV]");
  xaxis->SetLimits(500,8200); 
  yaxis->SetRangeUser(0.35,0.70);
  yaxis->SetTitle("#varepsilon_{W-tag}");
  yaxis->SetTitleOffset(1.2);
  gr1->Draw("A4 L SAME");
  // f3->Draw("SAME");
  // f4->Draw("SAME");
  // f5->Draw("SAME");
  // f6->Draw("SAME");
  lumi_13TeV = "";
  CMS_lumi(c1,5,iPos);
  leg->Clear();
  leg->AddEntry(gr1,"W-tagging efficiency average","fl");
  // leg->AddEntry(f3,"Spin-0, #Gamma_{X} / m_{X} = 0.01%","lep");
  // leg->AddEntry(f4,"Spin-0, #Gamma_{X} / m_{X} = 5%","lep");
  // leg->AddEntry(f5,"Spin-1, #Gamma_{X} / m_{X} = 0.01%","lep");
  // leg->AddEntry(f6,"Spin-1, #Gamma_{X} / m_{X} = 5%","lep");
  leg->Draw();
  c1->Print("EffWtag_17.pdf");
  c1->Print("EffWtag_17.svg");
  
}