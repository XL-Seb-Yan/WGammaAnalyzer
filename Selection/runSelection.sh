#!/bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/work/x/xuyan/work5/CMSSW_9_4_0/src/VgammaTuplizer/Analyzer/flat

# integrated luminosity for data
LUMI=2215

root -l -q select.C+\(\"muon.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d *.pcm
