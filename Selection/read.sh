name1="/afs/cern.ch/work/x/xuyan/work5/CMSSW_9_4_0/src/VgammaTuplizer/Ntuplizer/crab_jobs_2017B_Apr19/crab_Wgamma94XSingleMuonTuples_Apr19_2017B/results/"
name2=" 0 NONE"
files=$(ls /afs/cern.ch/work/x/xuyan/work5/CMSSW_9_4_0/src/VgammaTuplizer/Ntuplizer/crab_jobs_2017B_Apr19/crab_Wgamma94XSingleMuonTuples_Apr19_2017B/results)
for line in $files
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> muon.conf
done
