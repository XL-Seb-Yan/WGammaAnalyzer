#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

root -l -q select.C+\(\"sampleswide.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d *.pcm
