#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

#weight:
#GJets100-200: 36.03
#GJets200-400: 4.89
#GJets400-600: 2.32
#GJets600-Inf: 1.08
#QCD300-500: 222.10
#QCD500-700: 22.15
#QCD700-1000: 5.514
#QCD1000-1500: 2.72
#QCD1500-2000: 0.35
#QCD2000-Inf: 0.14
root -l -q ../select_trigger1718.C+\(\"sample_sig_SF.conf\",\"${NTUPDIR}\",\"Signal\"\)
