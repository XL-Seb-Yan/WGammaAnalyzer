#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

# integrated luminosity for data
LUMI=2215

root -l -q select.C+\(\"samples.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d *.pcm
