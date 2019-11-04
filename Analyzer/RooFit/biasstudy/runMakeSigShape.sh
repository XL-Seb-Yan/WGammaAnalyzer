#!/bin/bash

for i in {3500,4000}
do
    root -l -q make_signal_v2_shapes.C+\(${i}\)
done
rm *.so *.d *.pcm
