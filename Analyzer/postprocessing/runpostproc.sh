#!/bin/bash
# for i in {"GJets_HT-100To200","GJets_HT-200To400","GJets_HT-400To600","GJets_HT-600ToInf","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"}
# do
	# root -l -q postproc.C+\(\"${i}\",0,17\)
# done

for i in {"SinglePhoton2016",""}
do
		root -l -q postproc.C+\(\"${i}\",1,16\)
done


# for i in {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500}
 # for i in {700,1200,2000,2800,3500}
# do
	# root -l -q postproc.C+\(${i},0,17\)
# done