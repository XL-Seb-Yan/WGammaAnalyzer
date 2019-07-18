#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

root -l -q select_TMVA.C+\(\"SinglePhoton_2017C_Jul06.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d *.pcm
