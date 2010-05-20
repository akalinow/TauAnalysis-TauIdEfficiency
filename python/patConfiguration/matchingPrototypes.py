import FWCore.ParameterSet.Config as cms

'''

Prototypes of the various and cleaning modules used to produce
the relevant collections of taus for different measurement
channels

Author: Evan K. Friis (UC Davis)

TODO: ppMuX matcher

'''

# Cleaner is used to match pat::Taus to different probe jets
dijetCleanerPrototype = cms.EDProducer(
    "PATTauCleaner",
    src = cms.InputTag("myPatTaus"),
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        # check overlaps on tag jet
        TagJet = cms.PSet(
            src = cms.InputTag("pfJetsTagAndProbes", "tagObject"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            deltaR = cms.double(0.3),
            # don't check if they share some AOD object ref
            checkRecoComponents = cms.bool(False), 
            pairCut = cms.string(""),
            # overlaps don't cause the electron to be discared
            requireNoOverlaps = cms.bool(False), 
        ),
        HighestPtProbeJet = cms.PSet(
            src = cms.InputTag("pfJetsTagAndProbes", "highestPtProbe"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            deltaR = cms.double(0.3),
            # don't check if they share some AOD object ref
            checkRecoComponents = cms.bool(False), 
            pairCut = cms.string(""),
            # overlaps don't cause the electron to be discared
            requireNoOverlaps = cms.bool(False), 
        ),
        SecondHighestPtProbe = cms.PSet(
            src = cms.InputTag("pfJetsTagAndProbes", "secondHighestPtProbe"),
            algorithm = cms.string("byDeltaR"),
            preselection = cms.string(""),
            deltaR = cms.double(0.3),
            # don't check if they share some AOD object ref
            checkRecoComponents = cms.bool(False), 
            pairCut = cms.string(""),
            # overlaps don't cause the electron to be discared
            requireNoOverlaps = cms.bool(False), 
        )
    ),
    finalCut = cms.string('')
)

