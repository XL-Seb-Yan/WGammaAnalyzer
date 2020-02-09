name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan16_2016B/200117_005809/0000/"
name2=" 0 NONE"
name0="$ 2016SingleMuonB"
files=$(ls $name1)
echo ${name0} >> 2016SingleMuonB.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2016SingleMuonB.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan16_2016C/200117_010010/0000/"
name2=" 0 NONE"
name0="$ 2016SingleMuonC"
files=$(ls $name1)
echo ${name0} >> 2016SingleMuonC.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2016SingleMuonC.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan16_2016D/200117_010214/0000/"
name2=" 0 NONE"
name0="$ 2016SingleMuonD"
files=$(ls $name1)
echo ${name0} >> 2016SingleMuonD.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2016SingleMuonD.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan16_2016E/200117_010424/0000/"
name2=" 0 NONE"
name0="$ 2016SingleMuonE"
files=$(ls $name1)
echo ${name0} >> 2016SingleMuonE.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2016SingleMuonE.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan16_2016F/200117_010633/0000/"
name2=" 0 NONE"
name0="$ 2016SingleMuonF"
files=$(ls $name1)
echo ${name0} >> 2016SingleMuonF.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2016SingleMuonF.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan16_2016G/200117_010837/0000/"
name2=" 0 NONE"
name0="$ 2016SingleMuonG"
files=$(ls $name1)
echo ${name0} >> 2016SingleMuonG.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2016SingleMuonG.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan16_2016H/200117_011040/0000/"
name2=" 0 NONE"
name0="$ 2016SingleMuonH"
files=$(ls $name1)
echo ${name0} >> 2016SingleMuonH.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2016SingleMuonH.conf
done
