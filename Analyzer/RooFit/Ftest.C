#include <TMath.h>
#include <TLegend.h>
void Ftest()
{
  gROOT->SetBatch(1);
  double SSR1 = 520.72;
  double SSR2 = 481.09;
  int p1 = 2;
  int p2 = 3;
  int n = 120;
  double f_statistic = ((SSR1-SSR2)/(p2-p1))/((SSR2)/(n-p2));
  cout<<"f statistic is: "<<f_statistic<<endl;
  double p_value = ROOT::Math::fdistribution_cdf_c(f_statistic,p2-p1,n-p2);
  cout<<"p value is: "<<p_value<<endl;
}
