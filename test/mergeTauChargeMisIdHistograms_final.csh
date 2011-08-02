#!/bin/csh -f

set fitResults                        = '/data1/veelken/tmp/muonPtGt20/V6/fitTauIdEff.root'
set compTauIdEffPreselNumbers_Ztautau = '/data1/veelken/tmp/muonPtGt20/V6/compTauIdEffPreselNumbers_Ztautau_powheg_2011Jul23V6_noTauSel_tauChargeMisIdRate.root'
set compTauIdEffFinalNumbers          = '/data1/veelken/tmp/muonPtGt20/V6/compTauChargeMisIdFinalNumbers_input.root'

#echo hadd -f $compTauIdEffFinalNumbers $fitResults $compTauIdEffPreselNumbers_Ztautau
hadd -f $compTauIdEffFinalNumbers $fitResults $compTauIdEffPreselNumbers_Ztautau
