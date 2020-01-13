name1="/isilon/hadoop/store/user/xuyan/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/Wgamma949_QCD_HT300to500_Aug20/190904_224519/0000/"
name2=" 0 NONE"
name0="$ QCD_HT300to500"
files=$(ls /isilon/hadoop/store/user/xuyan/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/Wgamma949_QCD_HT300to500_Aug20/190904_224519/0000)
echo ${name0} >> QCD_HT300to500.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT300to500.conf
done
'
name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Aug07_2017E/190807_143226/0001/"
files=$(ls /isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Aug07_2017E/190807_143226/0001)
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton_2017E_Aug07.conf
done
'

