#!/bin/bash

for i in {800,900,1000}
#for i in {1200,1400,1600,1800,2000}
#for i in {2200,2400,2600,2800,3000}
#for i in {3500,4000}
do
    root -l -q make_signal_v2_shapes.C+\(${i}\)
done
rm *.so *.d *.pcm
