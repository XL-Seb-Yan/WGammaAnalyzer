name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017B/200121_090937/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017B"
files=$(ls $name1)
echo ${name0} >> SinglePhoton2017B.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton2017B.conf
done

name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017C/200121_091116/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017C"
files=$(ls $name1)
echo ${name0} >> SinglePhoton2017C.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton2017C.conf
done

name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017D/200121_091257/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017D"
files=$(ls $name1)
echo ${name0} >> SinglePhoton2017D.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton2017D.conf
done

name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017E/200121_091431/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017E"
files=$(ls $name1)
echo ${name0} >> SinglePhoton2017E.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton2017E.conf
done

name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017F/200121_091605/0000/"
name2=" 0 NONE"
name0="$ SinglePhoton2017F"
files=$(ls $name1)
echo ${name0} >> SinglePhoton2017F.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> SinglePhoton2017F.conf
done

