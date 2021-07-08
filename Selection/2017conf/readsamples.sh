name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017E/200406_090735/0000/"
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
'
name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017B/200406_090501/0000/"
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

name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017C/200406_090552/0000/"
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

name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017D/200406_090645/0000/"
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

name1="/isilon/hadoop/store/user/xuyan/SinglePhoton/Wgamma94XSinglePhoton_Jan10_2017F/200406_090830/0000/"
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
'
