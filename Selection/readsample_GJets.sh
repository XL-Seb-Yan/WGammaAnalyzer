name1="/isilon/hadoop/store/user/xuyan/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/Wgamma949_GJets_HT-100To200_Jan10/200110_203037/0000/"
name2=" 0 NONE"
name0="$ QCD_HT-100To200"
files=$(ls $name1)
echo ${name0} >> QCD_HT-100To200.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT-100To200.conf
done

name1="/isilon/hadoop/store/user/xuyan/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/Wgamma949_GJets_HT-200To400_Jan10/200110_203215/0000/"
name2=" 0 NONE"
name0="$ QCD_HT-200To400"
files=$(ls $name1)
echo ${name0} >> QCD_HT-200To400.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT-200To400.conf
done

name1="/isilon/hadoop/store/user/xuyan/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/Wgamma949_GJets_HT-400To600_Jan10/200110_203401/0000/"
name2=" 0 NONE"
name0="$ QCD_HT700to1000"
files=$(ls $name1)
echo ${name0} >> QCD_HT-400To600.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT-400To600.conf
done

name1="/isilon/hadoop/store/user/xuyan/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/Wgamma949_GJets_HT-600ToInf_Jan10/200110_203541/0000/"
name2=" 0 NONE"
name0="$ QCD_HT1000to1500"
files=$(ls $name1)
echo ${name0} >> QCD_HT-600ToInf.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT-600ToInf.conf
done
