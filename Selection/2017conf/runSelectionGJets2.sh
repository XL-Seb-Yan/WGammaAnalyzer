#!/bin/bash

# output ntuple directory
NTUPDIR=SelOutPut

#weight:
#GJets100-200: 36.02895
#GJets200-400: 4.895264
#GJets400-600: 2.322734
#GJets600-Inf: 1.076998
#QCD300-500: 222.6718
#QCD500-700: 22.26943
#QCD700-1000: 5.490946
#QCD1000-1500: 2.737698
#QCD1500-2000: 0.354744
#QCD2000-Inf: 0.141199
root -l -q ../select201817_Final.C+\(\"GJets_HT-400To600.conf\",\"${NTUPDIR}\",2.322734\)
root -l -q ../select201817_Final.C+\(\"GJets_HT-600ToInf.conf\",\"${NTUPDIR}\",1.076998\)
