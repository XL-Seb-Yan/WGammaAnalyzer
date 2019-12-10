#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

#weight:
#GJets100-200: 34.6195
#GJets200-400: 4.9300
#GJets400-600: 2.1225
#GJets600-Inf: 1.0607
#QCD300-500: 214.8054
#QCD500-700: 2.4840
#QCD700-1000: 5.1894
#QCD1000-1500: 2.5156
#QCD1500-2000: 0.3635
#QCD2000-Inf: 0.1436
root -l -q select.C+\(\"sampleswide_temp.conf\",\"${NTUPDIR}\",1\)

rm *.so *.d *.pcm
