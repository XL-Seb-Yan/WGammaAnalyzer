#include <TMath.h>
#include <TLegend.h>
void LikelihoodRatiotest()
{
  gROOT->SetBatch(1);
  double nnl1 = 2968.725;
  double nnl2 = 2968.479;
  int p1=1;
  int p2=2;
  double LR = 2*(nnl1-nnl2);
  cout<<"likelihood ratio is: "<<LR<<endl;
  double p_value = ROOT::Math::chisquared_cdf_c(LR,p2-p1);
  cout<<"p value is: "<<p_value<<endl;
}
