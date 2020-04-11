#!/bin/bash

#for i in {"GJets_HT-100To200","GJets_HT-200To400","GJets_HT-400To600","GJets_HT-600ToInf","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"}
#for i in {"SinglePhoton2016","123"}
for i in {"EGamma2018","123"}
do
  root -l -q postproc.C+\(\"${i}\",1\)
done
rm *.so *.d *.pcm