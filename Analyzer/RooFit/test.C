/// \file
/// \ingroup tutorial_math
/// \notebook
/// Example of CrystalBall Function and its distribution (pdf and cdf)
///
/// \macro_image
/// \macro_code
///
/// \author Lorenzo Moneta

void test()  {
  gROOT->SetBatch(1);
  auto c1 = new TCanvas("1","1",1200,900);
   // crystal ball function
   c1->cd();

   //(alpha, n sigma, mu)
   auto f1 = new TF1("pdf","ROOT::Math::crystalball_pdf(x, 1, 1.01, 1, 0)",-5,5);
   f1->SetLineColor(kRed);
   f1->Draw();
   auto f2 = new TF1("pdf","ROOT::Math::crystalball_pdf(x, 1, 1.21, 1, 0)",-5,5);
   f2->SetLineColor(kGreen);
   f2->Draw("same");
   auto f3 = new TF1("pdf","ROOT::Math::crystalball_pdf(x, 1, 1.41, 1, 0)",-5,5);
   f3->SetLineColor(kBlue);
   f3->Draw("same");

   auto legend = new TLegend(0.5,0.7,0.6,0.8);
   legend->AddEntry(f1,"a=1","L");
   legend->AddEntry(f2,"a=2","L");
   legend->AddEntry(f3,"a=3","L");
   legend->Draw();
   TAxis *axis = f1->GetYaxis();
   axis->SetRangeUser(0,0.3);

   c1->Print("test2.png");
}
