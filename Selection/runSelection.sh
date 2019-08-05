#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

root -l -q select_TMVA.C+\(\"samples.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d *.pcm
