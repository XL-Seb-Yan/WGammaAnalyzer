#!/bin/bash

#for i in {750,850,950,1050,1100,1150,1250,1300,1350,1450,1500,1550,1600,1650,1700,1750,1850,1900,1950,2050,2100,2150,2250,2300,2350,2450,2500,2550,2650,2700,2750,2850,2900,2950,3050,3100,3150,3200,3250,3300,3350,3400,3450}
#for i in {3050,3100,3150,3200,3250,3300,3350,3400,3450}
#for i in {700,800,900}
#for i in {1000,1200,1400,1600,1800}
#for i in {2200,2400,2600}
for i in {2800,3000,3500}
do
    root -l -q make_signal_narrow_shapes_Bukin.C+\(${i},300\)
done
rm *.so *.d *.pcm
