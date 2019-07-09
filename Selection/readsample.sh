name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jul06_2017D/190706_205102/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017D"
files=$(ls /isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jul06_2017D/190706_205102/0000)
echo ${name0} >> SinglePhoton_2017D_Jul06.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton_2017D_Jul06.conf
done
