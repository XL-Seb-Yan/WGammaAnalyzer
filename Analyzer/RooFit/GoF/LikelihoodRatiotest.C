#include <TMath.h>
#include <TLegend.h>
void LikelihoodRatiotest()
{
  gROOT->SetBatch(1);
  double nnl1 = 6386.1933;
  double nnl2 = 6385.0444;
  int p1=1;
  int p2=2;
  double LR = 2*(nnl1-nnl2);
  cout<<"likelihood ratio is: "<<LR<<endl;
  double p_value = ROOT::Math::chisquared_cdf_c(LR,p2-p1);
  cout<<"p value is: "<<p_value<<endl;
}
