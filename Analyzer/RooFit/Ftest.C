#include <TMath.h>
#include <TLegend.h>
void Ftest()
{
  gROOT->SetBatch(1);
  double SSR1 = 56481.336353;
  double SSR2 = 56102.489294;
  int p1 = 3;
  int p2 = 4;
  int n = 95;
  double f_statistic = ((SSR1-SSR2)/(p2-p1))/((SSR2)/(n-p2));
  cout<<"f statistic is: "<<f_statistic<<endl;
  double p_value = ROOT::Math::fdistribution_cdf_c(f_statistic,p2-p1,n-p2);
  cout<<"p value is: "<<p_value<<endl;
}
