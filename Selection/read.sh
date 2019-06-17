name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhotonTuplesPhotonrich_Jun16_2017B/190616_083801/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017B"
files=$(ls /isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhotonTuplesPhotonrich_Jun16_2017B/190616_083801/0000)
echo ${name0} >> SinglePhoton_2017B_noHLTFilter.conf
for line in $files
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton_2017B_noHLTFilter.conf
done
