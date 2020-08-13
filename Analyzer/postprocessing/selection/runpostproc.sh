#!/bin/bash
# for i in {"QCD","GJets"}
# do
	# root -l -q postproc.C+\(\"${i}\",0,17\)
# done

#for i in {"SinglePhoton2016",""}
#for i in {"SinglePhoton2017",""}
 # for i in {"EGamma2018",""}
# do
		# root -l -q postproc.C+\(\"${i}\",1,18\)
# done


#for i in {700,2000}
for i in {7000,8000}
do
	root -l -q postproc.C+\(${i},0,17\)
done