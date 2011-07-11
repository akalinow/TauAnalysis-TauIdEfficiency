#!/bin/csh -f

set fitResults                        = '/data1/veelken/tmp/muonPtGt20/V4b/fitTauIdEff_wConstraints.root'
set compTauIdEffPreselNumbers_Ztautau = '/data1/veelken/tmp/muonPtGt20/V4b/compTauIdEffPreselNumbers_Ztautau_pythia_2011Jul06_mauroV4_noTauSel.root'
#set compTauIdEffPreselNumbers_Ztautau = '/data1/veelken/tmp/muonPtGt20/V4b/compTauIdEffPreselNumbers_Ztautau_2011Jun30v2.root'

set compTauIdEffFinalNumbers          = '/data1/veelken/tmp/muonPtGt20/V4b/compTauIdEffFinalNumbers_input.root'

#echo hadd $compTauIdEffFinalNumbers $fitResults $compTauIdEffPreselNumbers_Ztautau
hadd $compTauIdEffFinalNumbers $fitResults $compTauIdEffPreselNumbers_Ztautau
