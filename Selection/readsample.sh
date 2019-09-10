name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Aug20_2017C/190819_204741/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017C"
files=$(ls /isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Aug20_2017C/190819_204741/0000)
echo ${name0} >> SinglePhoton_2017C_Aug20.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton_2017C_Aug20.conf
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

