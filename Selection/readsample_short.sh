name1="/isilon/hadoop/store/user/xuyan/EGamma/Wgamma102XEGamma_Jan10v2_2018A/200114_121945/0000/"
name2=" 0 NONE"
name0="$ EGamma2018A"
files=$(ls $name1)
echo ${name0} >> EGamma2018A.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> EGamma2018A.conf
done

name1="/isilon/hadoop/store/user/xuyan/EGamma/Wgamma102XEGamma_Jan10v2_2018B/200114_122121/0000/"
name2=" 0 NONE"
name0="$ EGamma2018B"
files=$(ls $name1)
echo ${name0} >> EGamma2018B.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> EGamma2018B.conf
done

name1="/isilon/hadoop/store/user/xuyan/EGamma/Wgamma102XEGamma_Jan10v2_2018C/200114_122254/0000/"
name2=" 0 NONE"
name0="$ EGamma2018C"
files=$(ls $name1)
echo ${name0} >> EGamma2018C.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> EGamma2018C.conf
done

name1="/isilon/hadoop/store/user/xuyan/EGamma/Wgamma102XEGamma_Jan10v2_2018D/200114_122428/0000/"
name2=" 0 NONE"
name0="$ EGamma2018D"
files=$(ls $name1)
echo ${name0} >> EGamma2018D.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> EGamma2018D.conf
done

