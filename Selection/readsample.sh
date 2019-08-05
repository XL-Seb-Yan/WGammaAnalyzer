name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Aug02_2017F/190802_223739/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017F"
files=$(ls /isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Aug02_2017F/190802_223739/0000)
echo ${name0} >> SinglePhoton_2017F_Aug02.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton_2017F_Aug02.conf
done

name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Aug02_2017F/190802_223739/0001/"
files=$(ls /isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Aug02_2017F/190802_223739/0001)
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton_2017F_Aug02.conf
done

