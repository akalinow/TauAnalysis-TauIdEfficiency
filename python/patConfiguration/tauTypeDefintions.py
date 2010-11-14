import FWCore.ParameterSet.Config as cms

'''
Dictionary matching PATTauProducer parameters to values
appropriate for the different tau collections

Author: Evan K. Friis (UC Davis)

'''

_shrinkingConeOptions = {
    # off for now so we dont' need to rerun PFT sequence
    'tauSource' : cms.InputTag('shrinkingConePFTauProducer'),
    'tauIDSources' : cms.PSet(
        leadingTrackFinding = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut"),
        leadingPionPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingPionPtCut"),
        trackIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation"),
        trackIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion"),
        ecalIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation"),
        ecalIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion"),
        byIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation"),
        byIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion"),
        againstElectron = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon"),
        byTaNC = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"),
        byTaNCfrOnePercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrOnePercent"),
        byTaNCfrHalfPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrHalfPercent"),
        byTaNCfrQuarterPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent"),
        byTaNCfrTenthPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrTenthPercent")
    )
}

_fixedConeOptions = {
    'tauSource' : cms.InputTag('fixedConePFTauProducer'),
    'tauIDSources' : cms.PSet(
        leadingTrackFinding = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackFinding"),
        leadingTrackPtCut = cms.InputTag("fixedConePFTauDiscriminationByLeadingTrackPtCut"),
        leadingPionPtCut = cms.InputTag("fixedConePFTauDiscriminationByLeadingPionPtCut"),
        trackIsolation = cms.InputTag("fixedConePFTauDiscriminationByTrackIsolation"),
        trackIsolationUsingLeadingPion = cms.InputTag("fixedConePFTauDiscriminationByTrackIsolationUsingLeadingPion"),
        ecalIsolation = cms.InputTag("fixedConePFTauDiscriminationByECALIsolation"),
        ecalIsolationUsingLeadingPion = cms.InputTag("fixedConePFTauDiscriminationByECALIsolationUsingLeadingPion"),
        byIsolation = cms.InputTag("fixedConePFTauDiscriminationByIsolation"),
        byIsolationUsingLeadingPion = cms.InputTag("fixedConePFTauDiscriminationByIsolationUsingLeadingPion"),
        againstElectron = cms.InputTag("fixedConePFTauDiscriminationAgainstElectron"),
        againstMuon = cms.InputTag("fixedConePFTauDiscriminationAgainstMuon"),
    )
}

patTauProducerOptions = {
    'shrinkingCone' : _shrinkingConeOptions,
    'fixedCone' : _fixedConeOptions
}




