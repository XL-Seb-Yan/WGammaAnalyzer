name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan10_2017B/200110_232359/0000/"
name2=" 0 NONE"
name0="$ 2017SingleMuonB"
files=$(ls $name1)
echo ${name0} >> 2017SingleMuonB.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2017SingleMuonB.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan10_2017C/200110_232802/0000/"
name2=" 0 NONE"
name0="$ 2017SingleMuonC"
files=$(ls $name1)
echo ${name0} >> 2017SingleMuonC.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2017SingleMuonC.conf
done
name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan10_2017C/200110_232802/0001/"
files=$(ls $name1)
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2017SingleMuonC.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan10_2017D/200110_233213/0000/"
name2=" 0 NONE"
name0="$ 2017SingleMuonD"
files=$(ls $name1)
echo ${name0} >> 2017SingleMuonD.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2017SingleMuonD.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan10_2017E/200110_233708/0000/"
name2=" 0 NONE"
name0="$ 2017SingleMuonE"
files=$(ls $name1)
echo ${name0} >> 2017SingleMuonE.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2017SingleMuonE.conf
done
name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan10_2017E/200110_233708/0001/"
files=$(ls $name1)
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2017SingleMuonE.conf
done

name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan10_2017F/200110_234115/0000/"
name2=" 0 NONE"
name0="$ 2017SingleMuonF"
files=$(ls $name1)
echo ${name0} >> 2017SingleMuonF.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2017SingleMuonF.conf
done
name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuon_Jan10_2017F/200110_234115/0001/"
files=$(ls $name1)
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> 2017SingleMuonF.conf
done
