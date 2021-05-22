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

# for i in {700,1500,2000,2500,3000,4000,6000,8000}
for i in {6000,6500}
do
	root -l -q postproc.C+\(${i},0,17\)
done
rm *.d
rm *.so
rm *.pcm