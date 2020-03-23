//Plot histograms with scale factors applied directly on historgams
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/tdrstyle.C"
#include "/afs/cern.ch/work/x/xuyan/work5/PROD17/AN/AN-19-280/utils/general/CMS_lumi.C"
void Histoplot()
{

  gROOT->SetBatch(1);
  lumi_13TeV = "41.53 fb^{-1}";
  writeExtraText = 1;
  lumiTextOffset = 0.15;
  bool plot_CMS = true;
  extraText = "Preliminary";
  lumiTextSize = 0.35;
  cmsTextSize = 0.45;
  int iPeriod = 4;
  int iPos = 11;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetBarWidth(2);
  gStyle->SetHistLineWidth(2);

  // Write histos to root file
  TFile *file1 = TFile::Open("Histogram_M2800W.root");
  TFile *file2 = TFile::Open("Histogram_GJets_tau21.root");
  TFile *file3 = TFile::Open("Histogram_QCD_tau21.root");
  TH1* hist1_1 = (TH1*)file1->Get("1"); //p_pt
  TH1* hist1_2 = (TH1*)file1->Get("2"); //p_eta
  TH1* hist1_3 = (TH1*)file1->Get("3"); //j_pt
  TH1* hist1_4 = (TH1*)file1->Get("4"); //j_eta
  TH1* hist1_5 = (TH1*)file1->Get("5"); //j_e
  TH1* hist1_6 = (TH1*)file1->Get("6"); //j_mass
  TH1* hist1_7 = (TH1*)file1->Get("7"); //j_tau21
  TH1* hist1_8 = (TH1*)file1->Get("8"); //s_cos
  TH1* hist1_9 = (TH1*)file1->Get("9"); //s_ptm
  TH1* hist1_10 = (TH1*)file1->Get("10"); //s_invmass
  TH1* hist2_1 = (TH1*)file2->Get("1");
  TH1* hist2_2 = (TH1*)file2->Get("2");
  TH1* hist2_3 = (TH1*)file2->Get("3");
  TH1* hist2_4 = (TH1*)file2->Get("4");
  TH1* hist2_5 = (TH1*)file2->Get("5");
  TH1* hist2_6 = (TH1*)file2->Get("6");
  TH1* hist2_7 = (TH1*)file2->Get("7");
  TH1* hist2_8 = (TH1*)file2->Get("8");
  TH1* hist2_9 = (TH1*)file2->Get("9");
  TH1* hist2_10 = (TH1*)file2->Get("10");
  TH1* hist3_1 = (TH1*)file3->Get("1");
  TH1* hist3_2 = (TH1*)file3->Get("2");
  TH1* hist3_3 = (TH1*)file3->Get("3");
  TH1* hist3_4 = (TH1*)file3->Get("4");
  TH1* hist3_5 = (TH1*)file3->Get("5");
  TH1* hist3_6 = (TH1*)file3->Get("6");
  TH1* hist3_7 = (TH1*)file3->Get("7");
  TH1* hist3_8 = (TH1*)file3->Get("8");
  TH1* hist3_9 = (TH1*)file3->Get("9");
  TH1* hist3_10 = (TH1*)file3->Get("10");

  hist2_1->Scale(1.39);
  hist2_2->Scale(1.39);
  hist2_3->Scale(1.39);
  hist2_4->Scale(1.39);
  hist2_5->Scale(1.39);
  hist2_6->Scale(1.39);
  hist2_7->Scale(1.39);
  hist2_8->Scale(1.39);
  hist2_9->Scale(1.39);
  hist2_10->Scale(1.39);
  hist3_1->Scale(0.5);
  hist3_2->Scale(0.5);
  hist3_3->Scale(0.5);
  hist3_4->Scale(0.5);
  hist3_5->Scale(0.5);
  hist3_6->Scale(0.5);
  hist3_7->Scale(0.5);
  hist3_8->Scale(0.5);
  hist3_9->Scale(0.5);
  hist3_10->Scale(0.5);

  int color1 = 2;
  int color2 = 8;
  int color3 = 7;

  hist1_1->SetLineColor(1);
  hist1_2->SetLineColor(1);
  hist1_3->SetLineColor(1);
  hist1_4->SetLineColor(1);
  hist1_5->SetLineColor(1);
  hist1_6->SetLineColor(1);
  hist1_7->SetLineColor(1);
  hist1_8->SetLineColor(1);
  hist1_9->SetLineColor(1);
  hist1_10->SetLineColor(1);
  hist2_1->SetLineColor(1);
  hist2_2->SetLineColor(1);
  hist2_3->SetLineColor(1);
  hist2_4->SetLineColor(1);
  hist2_5->SetLineColor(1);
  hist2_6->SetLineColor(1);
  hist2_7->SetLineColor(1);
  hist2_8->SetLineColor(1);
  hist2_9->SetLineColor(1);
  hist2_10->SetLineColor(1);
  hist3_1->SetLineColor(1);
  hist3_2->SetLineColor(1);
  hist3_3->SetLineColor(1);
  hist3_4->SetLineColor(1);
  hist3_5->SetLineColor(1);
  hist3_6->SetLineColor(1);
  hist3_7->SetLineColor(1);
  hist3_8->SetLineColor(1);
  hist3_9->SetLineColor(1);
  hist3_10->SetLineColor(1);

  hist2_1->SetFillColor(color2);
  hist2_2->SetFillColor(color2);
  hist2_3->SetFillColor(color2);
  hist2_4->SetFillColor(color2);
  hist2_5->SetFillColor(color2);
  hist2_6->SetFillColor(color2);
  hist2_7->SetFillColor(color2);
  hist2_8->SetFillColor(color2);
  hist2_9->SetFillColor(color2);
  hist2_10->SetFillColor(color2);
  hist3_1->SetFillColor(color3);
  hist3_2->SetFillColor(color3);
  hist3_3->SetFillColor(color3);
  hist3_4->SetFillColor(color3);
  hist3_5->SetFillColor(color3);
  hist3_6->SetFillColor(color3);
  hist3_7->SetFillColor(color3);
  hist3_8->SetFillColor(color3);
  hist3_9->SetFillColor(color3);
  hist3_10->SetFillColor(color3);
  
  hist1_1->SetMarkerStyle(8);
  hist1_2->SetMarkerStyle(8);
  hist1_3->SetMarkerStyle(8);
  hist1_4->SetMarkerStyle(8);
  hist1_5->SetMarkerStyle(8);
  hist1_6->SetMarkerStyle(8);
  hist1_7->SetMarkerStyle(8);
  hist1_8->SetMarkerStyle(8);
  hist1_9->SetMarkerStyle(8);
  hist1_10->SetMarkerStyle(8);
  hist1_1->SetMarkerSize(1.5);
  hist1_2->SetMarkerSize(1.5);
  hist1_3->SetMarkerSize(1.5);
  hist1_4->SetMarkerSize(1.5);
  hist1_5->SetMarkerSize(1.5);
  hist1_6->SetMarkerSize(1.5);
  hist1_7->SetMarkerSize(1.5);
  hist1_8->SetMarkerSize(1.5);
  hist1_9->SetMarkerSize(1.5);
  hist1_10->SetMarkerSize(1.5);

  cout<<"OK";

  //Stacked plot

  THStack *stack1 = new THStack("stack1","pt_{p}");
  THStack *stack2 = new THStack("stack2","eta_{p}");
  THStack *stack3 = new THStack("stack3","pt_{j}");
  THStack *stack4 = new THStack("stack4","eta_{j}");
  THStack *stack5 = new THStack("stack5","e_{j}");
  THStack *stack6 = new THStack("stack6","Jet softdrop mass");
  THStack *stack7 = new THStack("stack7","Jet tau21");
  THStack *stack8 = new THStack("stack8","cos(#theta*)");
  THStack *stack9 = new THStack("stack9","pt_{j} / M");
  THStack *stack10 = new THStack("stack10","Invariant mass");
  stack1->Add(hist3_1); stack1->Add(hist2_1);
  stack2->Add(hist3_2); stack2->Add(hist2_2);
  stack3->Add(hist3_3); stack3->Add(hist2_3);
  stack4->Add(hist3_4); stack4->Add(hist2_4);
  stack5->Add(hist3_5); stack5->Add(hist2_5);
  stack6->Add(hist3_6); stack6->Add(hist2_6);
  stack7->Add(hist3_7); stack7->Add(hist2_7);
  stack8->Add(hist3_8); stack8->Add(hist2_8);
  stack9->Add(hist3_9); stack9->Add(hist2_9);
  stack10->Add(hist3_10); stack10->Add(hist2_10);

  TLegend *legend = new TLegend(0.54,0.78,0.9,0.9);
  // Residual plot
  TH1 *data = NULL;
  TH1 *bkg = NULL;
  TH1 *pull = NULL;
  TAxis *xaxis1 = NULL;
  TAxis *yaxis1 = NULL;
  TAxis *xaxis2 = NULL;
  TAxis *yaxis2 = NULL;

  //===========================================================
  //gStyle->SetHistMinimumZero();
  TCanvas *c01 = new TCanvas("c01","",2400,2200);
  c01->cd();
  TPad *p01a = new TPad("p01a","p01a",0.1,0.30,0.9,1.0);
  TPad *p01b = new TPad("p01b","p01b",0.1,0.1,0.9,0.315);
  p01a->Draw();
  p01b->Draw();
  p01a->cd();
  p01a->SetBottomMargin(0.11);
  p01a->SetLogy();
  xaxis1 = hist1_1->GetXaxis();
  yaxis1 = hist1_1->GetYaxis();
  xaxis1->SetTitle("pt_{#gamma} (GeV)");
  yaxis1->SetTitle("Entries / 60 GeV");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_1->Draw("E1");
  stack1->Draw("SAMEHIST");
  hist1_1->Draw("E1SAME");
  hist1_1->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_1,"M-2800 Wide","lep");
  legend->AddEntry(hist2_1,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_1,"2017 MC, QCD(weighted)","f");
  legend->Draw();
  
  p01b->cd();
  p01b->SetTopMargin(0.037);
  p01b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("pt_{#gamma} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p01b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c01->Print("p_pt.png");
  c01->Print("p_pt.pdf");
  c01->Print("p_pt.svg");
  c01->Print("p_pt.root");
  //==========================================================

  //===========================================================
  TCanvas *c02 = new TCanvas("c02","",2400,2200);
  c02->cd();
  TPad *p02a = new TPad("p02a","p02a",0.1,0.30,0.9,1.0);
  TPad *p02b = new TPad("p02b","p02b",0.1,0.1,0.9,0.315);
  p02a->Draw();
  p02b->Draw();
  p02a->cd();
  p02a->SetBottomMargin(0.11);
  p02a->SetLogy();
  xaxis1 = hist1_2->GetXaxis();
  yaxis1 = hist1_2->GetYaxis();
  xaxis1->SetTitle("#eta_{#gamma}");
  yaxis1->SetTitle("Entries / 0.08");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_2->Draw("E1");
  stack2->Draw("SAMEHIST");
  hist1_2->Draw("E1SAME");
  hist1_2->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_2,"M-2800 Wide","lep");
  legend->AddEntry(hist2_2,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_2,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p02b->cd();
  p02b->SetTopMargin(0.037);
  p02b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_2->Clone();
  bkg = (TH1*)hist2_2->Clone();
  bkg->Add((TH1*)hist3_2->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("#eta_{#gamma}");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p02b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c02->Print("p_eta.png");
  c02->Print("p_eta.pdf");
  c02->Print("p_eta.svg");
  c02->Print("p_eta.root");
  //==========================================================

  //===========================================================
  TCanvas *c03 = new TCanvas("c03","",2400,2200);
  c03->cd();
  TPad *p03a = new TPad("p03a","p03a",0.1,0.30,0.9,1.0);
  TPad *p03b = new TPad("p03b","p03b",0.1,0.1,0.9,0.315);
  p03a->Draw();
  p03b->Draw();
  p03a->cd();
  p03a->SetBottomMargin(0.11);
  p03a->SetLogy();
  xaxis1 = hist1_3->GetXaxis();
  yaxis1 = hist1_3->GetYaxis();
  xaxis1->SetTitle("pt_{AK8jet} (GeV)");
  yaxis1->SetTitle("Entries / 60 (GeV)");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_3->Draw("E1");
  stack3->Draw("SAMEHIST");
  hist1_3->Draw("E1SAME");
  hist1_3->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_3,"M-2800 Wide","lep");
  legend->AddEntry(hist2_3,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_3,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p03b->cd();
  p03b->SetTopMargin(0.037);
  p03b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_3->Clone();
  bkg = (TH1*)hist2_3->Clone();
  bkg->Add((TH1*)hist3_3->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("pt_{AK8jet} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p03b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c03->Print("j_pt.png");
  c03->Print("j_pt.pdf");
  c03->Print("j_pt.svg");
  c03->Print("j_pt.root");
  //==========================================================

  //===========================================================
  TCanvas *c04 = new TCanvas("c04","",2400,2200);
  c04->cd();
  TPad *p04a = new TPad("p04a","p04a",0.1,0.30,0.9,1.0);
  TPad *p04b = new TPad("p04b","p04b",0.1,0.1,0.9,0.315);
  p04a->Draw();
  p04b->Draw();
  p04a->cd();
  p04a->SetBottomMargin(0.11);
  p04a->SetLogy();
  xaxis1 = hist1_4->GetXaxis();
  yaxis1 = hist1_4->GetYaxis();
  xaxis1->SetTitle("#eta_{AK8jet}");
  yaxis1->SetTitle("Entries / 0.08");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_4->Draw("E1");
  stack4->Draw("SAMEHIST");
  hist1_4->Draw("E1SAME");
  hist1_4->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_4,"M-2800 Wide","lep");
  legend->AddEntry(hist2_4,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_4,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p04b->cd();
  p04b->SetTopMargin(0.037);
  p04b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_4->Clone();
  bkg = (TH1*)hist2_4->Clone();
  bkg->Add((TH1*)hist3_4->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("#eta_{AK8jet}");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p04b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c04->Print("j_eta.png");
  c04->Print("j_eta.pdf");
  c04->Print("j_eta.svg");
  c04->Print("j_eta.root");
  //==========================================================

  //===========================================================
  TCanvas *c05 = new TCanvas("c05","",2400,2200);
  c05->cd();
  TPad *p05a = new TPad("p05a","p05a",0.1,0.30,0.9,1.0);
  TPad *p05b = new TPad("p05b","p05b",0.1,0.1,0.9,0.315);
  p05a->Draw();
  p05b->Draw();
  p05a->cd();
  p05a->SetBottomMargin(0.11);
  p05a->SetLogy();
  xaxis1 = hist1_5->GetXaxis();
  yaxis1 = hist1_5->GetYaxis();
  xaxis1->SetTitle("E_{AK8jet} (GeV)");
  yaxis1->SetTitle("Entries / 60 (GeV)");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_5->Draw("E1");
  stack5->Draw("SAMEHIST");
  hist1_5->Draw("E1SAME");
  hist1_5->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_5,"M-2800 Wide","lep");
  legend->AddEntry(hist2_5,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_5,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p05b->cd();
  p05b->SetTopMargin(0.037);
  p05b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_5->Clone();
  bkg = (TH1*)hist2_5->Clone();
  bkg->Add((TH1*)hist3_5->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("E_{AK8jet} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p05b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c05->Print("j_e.png");
  c05->Print("j_e.pdf");
  c05->Print("j_e.svg");
  c05->Print("j_e.root");
  //==========================================================

  //===========================================================
  TCanvas *c06 = new TCanvas("c06","",2400,2200);
  c06->cd();
  TPad *p06a = new TPad("p06a","p06a",0.1,0.30,0.9,1.0);
  TPad *p06b = new TPad("p06b","p06b",0.1,0.1,0.9,0.315);
  p06a->Draw();
  p06b->Draw();
  p06a->cd();
  p06a->SetBottomMargin(0.11);
  p06a->SetLogy();
  xaxis1 = hist1_6->GetXaxis();
  yaxis1 = hist1_6->GetYaxis();
  xaxis1->SetTitle("Softdrop m_{AK8jet} (GeV)");
  yaxis1->SetTitle("Entries / 2 (GeV)");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_6->Draw("E1");
  stack6->Draw("SAMEHIST");
  hist1_6->Draw("E1SAME");
  hist1_6->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_6,"M-2800 Wide","lep");
  legend->AddEntry(hist2_6,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_6,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p06b->cd();
  p06b->SetTopMargin(0.037);
  p06b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_6->Clone();
  bkg = (TH1*)hist2_6->Clone();
  bkg->Add((TH1*)hist3_6->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("Softdrop m_{AK8jet} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p06b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c06->Print("j_msd.png");
  c06->Print("j_msd.pdf");
  c06->Print("j_msd.svg");
  c06->Print("j_msd.root");
  //==========================================================

  //===========================================================
  TCanvas *c07 = new TCanvas("c07","",2400,2200);
  c07->cd();
  TPad *p07a = new TPad("p07a","p07a",0.1,0.30,0.9,1.0);
  TPad *p07b = new TPad("p07b","p07b",0.1,0.1,0.9,0.315);
  p07a->Draw();
  p07b->Draw();
  p07a->cd();
  p07a->SetBottomMargin(0.11);
  p07a->SetLogy();
  xaxis1 = hist1_7->GetXaxis();
  yaxis1 = hist1_7->GetYaxis();
  xaxis1->SetTitle("AK8 jet #tau_{21}");
  yaxis1->SetTitle("Entries / 0.02");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_7->Draw("E1");
  stack7->Draw("SAMEHIST");
  hist1_7->Draw("E1SAME");
  hist1_7->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_7,"M-2800 Wide","lep");
  legend->AddEntry(hist2_7,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_7,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p07b->cd();
  p07b->SetTopMargin(0.037);
  p07b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_7->Clone();
  bkg = (TH1*)hist2_7->Clone();
  bkg->Add((TH1*)hist3_7->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("AK8 jet #tau_{21}");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p07b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c07->Print("j_tau21.png");
  c07->Print("j_tau21.pdf");
  c07->Print("j_tau21.svg");
  c07->Print("j_tau21.root");
  //==========================================================

   //===========================================================
  TCanvas *c08 = new TCanvas("c08","",2400,2200);
  c08->cd();
  TPad *p08a = new TPad("p08a","p08a",0.1,0.30,0.9,1.0);
  TPad *p08b = new TPad("p08b","p08b",0.1,0.1,0.9,0.315);
  p08a->Draw();
  p08b->Draw();
  p08a->cd();
  p08a->SetBottomMargin(0.11);
  p08a->SetLogy();
  xaxis1 = hist1_8->GetXaxis();
  yaxis1 = hist1_8->GetYaxis();
  xaxis1->SetTitle("cos(#theta*)");
  yaxis1->SetTitle("Entries / 0.02");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_8->Draw("E1");
  stack8->Draw("SAMEHIST");
  hist1_8->Draw("E1SAME");
  hist1_8->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_8,"M-2800 Wide","lep");
  legend->AddEntry(hist2_8,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_8,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p08b->cd();
  p08b->SetTopMargin(0.037);
  p08b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_8->Clone();
  bkg = (TH1*)hist2_8->Clone();
  bkg->Add((TH1*)hist3_8->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("cos(#theta*)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p08b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c08->Print("s_cos.png");
  c08->Print("s_cos.pdf");
  c08->Print("s_cos.svg");
  c08->Print("s_cos.root");
  //==========================================================

     //===========================================================
  TCanvas *c09 = new TCanvas("c09","",2400,2200);
  c09->cd();
  TPad *p09a = new TPad("p09a","p09a",0.1,0.30,0.9,1.0);
  TPad *p09b = new TPad("p09b","p09b",0.1,0.1,0.9,0.315);
  p09a->Draw();
  p09b->Draw();
  p09a->cd();
  p09a->SetBottomMargin(0.11);
  p09a->SetLogy();
  xaxis1 = hist1_9->GetXaxis();
  yaxis1 = hist1_9->GetYaxis();
  xaxis1->SetTitle("pt_{#gamma}/M_{j#gamma}");
  yaxis1->SetTitle("Entries / 0.04");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_9->Draw("E1");
  stack9->Draw("SAMEHIST");
  hist1_9->Draw("E1SAME");
  hist1_9->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_9,"M-2800 Wide","lep");
  legend->AddEntry(hist2_9,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_9,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p09b->cd();
  p09b->SetTopMargin(0.037);
  p09b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_9->Clone();
  bkg = (TH1*)hist2_9->Clone();
  bkg->Add((TH1*)hist3_9->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("pt_{#gamma}/M_{j#gamma}");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p09b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c09->Print("s_ptm.png");
  c09->Print("s_ptms.pdf");
  c09->Print("s_ptm.svg");
  c09->Print("s_ptm.root");
  //==========================================================

  //===========================================================
  TCanvas *c10 = new TCanvas("c10","",2400,2200);
  c10->cd();
  TPad *p10a = new TPad("p10a","p10a",0.1,0.30,0.9,1.0);
  TPad *p10b = new TPad("p10b","p10b",0.1,0.1,0.9,0.315);
  p10a->Draw();
  p10b->Draw();
  p10a->cd();
  p10a->SetBottomMargin(0.11);
  p10a->SetLogy();
  xaxis1 = hist1_10->GetXaxis();
  yaxis1 = hist1_10->GetYaxis();
  xaxis1->SetTitle("M_{j#gamma} (GeV)");
  yaxis1->SetTitle("Entries / 80 GeV");
  xaxis1->SetTitleOffset(1.1);
  yaxis1->SetTitleOffset(1.3);
  yaxis1->SetRangeUser(0.01,1000000);
  hist1_10->Draw("E1");
  stack10->Draw("SAMEHIST");
  hist1_10->Draw("E1SAME");
  hist1_10->Draw("AXISSAME");
  legend->Clear();
  legend->AddEntry(hist1_10,"M-2800 Wide","lep");
  legend->AddEntry(hist2_10,"2017 MC, GJets(weighted)","f");
  legend->AddEntry(hist3_10,"2017 MC, QCD(weighted)","f");
  legend->Draw();

  p10b->cd();
  p10b->SetTopMargin(0.037);
  p10b->SetBottomMargin(0.5);
  /*
  pull = (TH1*)hist1_1->Clone();
  data = (TH1*)hist1_1->Clone();
  bkg = (TH1*)hist2_1->Clone();
  bkg->Add((TH1*)hist3_1->Clone());
  pull->Add(bkg,-1);
  pull->Divide(data);
  */
  pull = (TH1*)hist1_10->Clone();
  bkg = (TH1*)hist2_10->Clone();
  bkg->Add((TH1*)hist3_10->Clone());
  pull->Divide(bkg);
  xaxis2 = pull->GetXaxis();
  yaxis2 = pull->GetYaxis();
  xaxis2->SetTitle("M_{j#gamma} (GeV)");
  xaxis2->SetTitleOffset(1.25);
  yaxis2->SetTitle("data/MC");
  yaxis2->SetTitleOffset(0.5);
  yaxis2->SetRangeUser(0,2);
  xaxis2->SetLabelSize(0.1);
  xaxis2->SetTitleSize(0.1);
  yaxis2->SetLabelSize(0.1);
  yaxis2->SetNdivisions(5);
  yaxis2->SetTitleSize(0.08);
  p10b->SetGrid();
  p10b->SetGrid();
  pull->SetFillColor(kViolet);
  //pull->SetLineColor(kViolet);
  pull->SetFillColorAlpha(kViolet, 0.35);
  pull->Draw("PE1");
  //pull->Draw("BAR HIST");
  c10->Print("s_m.png");
  c10->Print("s_m.pdf");
  c10->Print("s_m.svg");
  c10->Print("s_m.root");
  //==========================================================
}
