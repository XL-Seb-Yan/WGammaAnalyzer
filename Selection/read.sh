name1="/isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuonTuples_Apr19_2017B/190419_053455/0000/"
name2=" 0 NONE"
name0="$ SingleMuon2017B"
files=$(ls /isilon/hadoop/store/user/xuyan/SingleMuon/Wgamma94XSingleMuonTuples_Apr19_2017B/190419_053455/0000)
echo ${name0} >> SingleMuon_2017B.conf
for line in $files
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SingleMuon_2017B.conf
done
