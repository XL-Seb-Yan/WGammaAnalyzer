name1="/home/xyan13/WGProj/CMSSW_9_4_0/src/WGammaAnalyzer/Selection/data/prescale/"
files=$(ls /home/xyan13/WGProj/CMSSW_9_4_0/src/WGammaAnalyzer/Selection/data/prescale)
for line in $files
do
echo $line
name3=${name1}${line}
echo ${name3} >> Photon175_json.txt
done
