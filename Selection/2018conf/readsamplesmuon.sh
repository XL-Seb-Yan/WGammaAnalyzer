name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma102XSingleMuon_Jan16_2018A/200125_115106/0000/"
name2=" 0 NONE"
name0="$ 2018SingleMuonA1"
files=$(ls $name1)
echo ${name0} >> 2018SingleMuonA1.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2018SingleMuonA1.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma102XSingleMuon_Jan16_2018A/200125_115106/0001/"
name2=" 0 NONE"
name0="$ 2018SingleMuonA2"
files=$(ls $name1)
echo ${name0} >> 2018SingleMuonA2.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2018SingleMuonA2.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma102XSingleMuon_Jan16_2018B/200125_115248/0000/"
name2=" 0 NONE"
name0="$ 2018SingleMuonB"
files=$(ls $name1)
echo ${name0} >> 2018SingleMuonB.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2018SingleMuonB.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma102XSingleMuon_Jan16_2018C/200125_115431/0000/"
name2=" 0 NONE"
name0="$ 2018SingleMuonC"
files=$(ls $name1)
echo ${name0} >> 2018SingleMuonC.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2018SingleMuonC.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma102XSingleMuon_Jan16_2018D/200125_114144/0000/"
name2=" 0 NONE"
name0="$ 2018SingleMuonD1"
files=$(ls $name1)
echo ${name0} >> 2018SingleMuonD1.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2018SingleMuonD1.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma102XSingleMuon_Jan16_2018D/200125_114144/0001/"
name2=" 0 NONE"
name0="$ 2018SingleMuonD2"
files=$(ls $name1)
echo ${name0} >> 2018SingleMuonD2.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2018SingleMuonD2.conf
done
