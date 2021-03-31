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
for i in {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500,4000,5000,6000,7000,8000}
do
	root -l -q postproc.C+\(${i},0,17\)
done
rm *.d
rm *.so
rm *.pcm