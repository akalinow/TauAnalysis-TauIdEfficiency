import FWCore.ParameterSet.Config as cms

#tauid pset {name, taucollection, ditaucollection, discName (pat discriminator name)} get them @ line 102. discName is what gets written in the tautuple, redundant? should change name to be = discName
#-------------------------------------------------------------------------------------------
#                                  OLD TANC
#-------------------------------------------------------------------------------------------
oldTancOnePercent = cms.PSet(
    name = cms.string("tauDiscrTaNCfrOnePercent"),
    taucollection = cms.string("selectedPatPFTausShrinkingConePFRelIsoCumulative"),
    ditaucollection = cms.string("selectedMuPFTauShrinkingConePairsForTauIdEffCumulative"),
    discName = cms.string("byTaNCfrOnePercent")
    )

oldTancHalfPercent = oldTancOnePercent.clone(
    name = cms.string("tauDiscrTaNCfrHalfPercent"),
    discName = cms.string("byTaNCfrHalfPercent")
    )

oldTancQuarterPercent = oldTancOnePercent.clone(
    name = cms.string("tauDiscrTaNCfrQuarterPercent"),
    discName = cms.string("byTaNCfrQuarterPercent")
    )

#-------------------------------------------------------------------------------------------
#                                  OLD HPS
#-------------------------------------------------------------------------------------------
oldHPSLoose = cms.PSet(
    name = cms.string("tauDiscrIsolationLoose"),
    taucollection = cms.string("selectedPatPFTausHPSPFRelIsoCumulative"),
    ditaucollection = cms.string("selectedMuPFTauHPSpairsForTauIdEffCumulative"),
    discName = cms.string("byIsolationLoose")
    )

oldHPSMedium = oldHPSLoose.clone(
    name = cms.string("tauDiscrIsolationMedium"),
    discName = cms.string("byIsolationMedium")
    )

oldHPSTight = oldHPSLoose.clone(
    name = cms.string("tauDiscrIsolationTight"),
    discName = cms.string("byIsolationTight")
    )

#-------------------------------------------------------------------------------------------
#                                  COMBINED HPS
#-------------------------------------------------------------------------------------------
combinedHPSLoose = cms.PSet(
    name = cms.string("tauDiscrHPSloose"),
    taucollection = cms.string("selectedPatPFTausHPSpTaNCPFRelIsoCumulative"),
    ditaucollection = cms.string("selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative"),
    discName = cms.string("byIsolationLoose")
    )

combinedHPSMedium = combinedHPSLoose.clone(
    name = cms.string("tauDiscrHPSmedium"),
    discName = cms.string("byIsolationMedium")
    )

combinedHPSTight = combinedHPSLoose.clone(
    name = cms.string("tauDiscrHPStight"),
    discName = cms.string("byIsolationTight")
    )

#-------------------------------------------------------------------------------------------
#                                  COMBINED TANC
#-------------------------------------------------------------------------------------------
combinedTancOnePercent = combinedHPSLoose.clone(
    name = cms.string("combinedTaNCOnePercent"),
    discName = cms.string("byTaNCfrOnePercent")
    )

combinedTancHalfPercent = combinedTancOnePercent.clone(
    name = cms.string("combinedTaNCHalfPercent"),
    discName = cms.string("byTaNCfrHalfPercent")
    )

combinedTancQuarterPercent = combinedTancOnePercent.clone(
    name = cms.string("combinedTaNCQuarterPercent"),
    discName = cms.string("byTaNCfrQuarterPercent")
    )
