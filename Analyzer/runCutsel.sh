#!/bin/bash

for i in {700,800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500}
do
    root -l -q cutselection.C+\(${i}\)
    #root -l -q effcal.C+\(${i}\)
done
rm *.so *.d *.pcm
