#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

# integrated luminosity for data
LUMI=2215

root -l -q select_trigger.C+\(\"SinglePhoton_2017B_Jun19.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d *.pcm
