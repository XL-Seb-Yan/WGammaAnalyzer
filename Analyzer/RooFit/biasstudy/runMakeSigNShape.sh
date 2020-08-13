#!/bin/bash

for i in {1000,1200,1400,1600,1800}
do
    #root -l -q GetNormalization.C+\(${i},500\)
	root -l -q make_signal_narrow_shapes_CBGaus.C+\(${i},600\)
done
rm *.so *.d *.pcm
