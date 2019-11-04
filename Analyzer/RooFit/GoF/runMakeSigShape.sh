#!/bin/bash

for i in {800,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000}
do
    root -l make_signal_shapes.C+\(${i}\)
done
rm *.so *.d *.pcm
