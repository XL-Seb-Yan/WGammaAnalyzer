#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

# integrated luminosity for data
LUMI=2215

root -l -q select_trig_debug.C+\(\"SingleMuon_2017B.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d *.pcm
