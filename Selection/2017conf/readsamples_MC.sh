name1="/isilon/hadoop/store/user/xuyan/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/Wgamma949_GJets_HT-100To200_Jan10/200403_082614/0000/"
name2=" 0 NONE"
name0="$ GJets_HT-100To200"
files=$(ls $name1)
echo ${name0} >> GJets_HT-100To200.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> GJets_HT-100To200.conf
done


name1="/isilon/hadoop/store/user/xuyan/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/Wgamma949_GJets_HT-200To400_Jan10/200403_082702/0000/"
name2=" 0 NONE"
name0="$ GJets_HT-200To400"
files=$(ls $name1)
echo ${name0} >> GJets_HT-200To400.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> GJets_HT-200To400.conf
done

name1="/isilon/hadoop/store/user/xuyan/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/Wgamma949_GJets_HT-400To600_Jan10/200403_082752/0000/"
name2=" 0 NONE"
name0="$ GJets_HT-400To600"
files=$(ls $name1)
echo ${name0} >> GJets_HT-400To600.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> GJets_HT-400To600.conf
done

name1="//isilon/hadoop/store/user/xuyan/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/Wgamma949_GJets_HT-600ToInf_Jan10/200403_082842/0000/"
name2=" 0 NONE"
name0="$ GJets_HT-600ToInf"
files=$(ls $name1)
echo ${name0} >> GJets_HT-600ToInf.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> GJets_HT-600ToInf.conf
done

name1="/isilon/hadoop/store/user/xuyan/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/Wgamma949_QCD_HT300to500_Jan10/200403_082115/0000/"
name2=" 0 NONE"
name0="$ QCD_HT300to500"
files=$(ls $name1)
echo ${name0} >> QCD_HT300to500.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT300to500.conf
done

name1="/isilon/hadoop/store/user/xuyan/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/Wgamma949_QCD_HT500to700_Jan10/200403_082205/0000/"
name2=" 0 NONE"
name0="$ QCD_HT500to700"
files=$(ls $name1)
echo ${name0} >> QCD_HT500to700.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT500to700.conf
done

name1="/isilon/hadoop/store/user/xuyan/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/Wgamma949_QCD_HT700to1000_Jan10/200403_082255/0000/"
name2=" 0 NONE"
name0="$ QCD_HT700to1000"
files=$(ls $name1)
echo ${name0} >> QCD_HT700to1000.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT700to1000.conf
done

name1="/isilon/hadoop/store/user/xuyan/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/Wgamma949_QCD_HT1000to1500_Jan10/200403_082345/0000/"
name2=" 0 NONE"
name0="$ QCD_HT1000to1500"
files=$(ls $name1)
echo ${name0} >> QCD_HT1000to1500.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT1000to1500.conf
done

name1="/isilon/hadoop/store/user/xuyan/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/Wgamma949_QCD_HT1500to2000_Jan10/200403_082435/0000/"
name2=" 0 NONE"
name0="$ QCD_HT1500to2000"
files=$(ls $name1)
echo ${name0} >> QCD_HT1500to2000.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT1500to2000.conf
done

name1="/isilon/hadoop/store/user/xuyan/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/Wgamma949_QCD_HT2000toInf_Jan10/200403_082525/0000/"
name2=" 0 NONE"
name0="$ QCD_HT2000toInf"
files=$(ls $name1)
echo ${name0} >> QCD_HT2000toInf.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> QCD_HT2000toInf.conf
done
