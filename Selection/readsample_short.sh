name1="/isilon/hadoop/store/user/xuyan/EGamma/Wgamma102XEGamma_Jan16_2018D/200116_173300/0002/"
name2=" 0 NONE"
name0="$ EGamma2018D5"
files=$(ls $name1)
echo ${name0} >> EGamma2018D5.conf
for line in $files 
do
echo $line
name3=${name1}${line}${name2}
echo ${name3} >> EGamma2018D5.conf
done

