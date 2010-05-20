import FWCore.ParameterSet.Config as cms

# Ntuple values to use for PAT taus built from PFTaus

variables = cms.PSet(
    # Signal occupancy
    nSignalCharged = cms.string('signalPFChargedHadrCands().size()'),
    nSignalGammas = cms.string('signalPFGammaCands().size()'),
    nSignalNeutralHadrons = cms.string('signalPFNeutrHadrCands().size()'),

    # Isolation occupancy
    nIsoCharged = cms.string('isolationPFChargedHadrCands().size()'),
    nIsoGammas = cms.string('isolationPFGammaCands().size()'),
    nIsoNeutralHadrons = cms.string('isolationPFNeutrHadrCands().size()'),

    # Quantitiies from underlying jet
    jetPt = cms.string('pfTauTagInfoRef().pfjetRef().pt()'),
    jetMass = cms.string('pfTauTagInfoRef().pfjetRef().mass()'),
    jetPhi = cms.string('pfTauTagInfoRef().pfjetRef().phi()'),
    jetEta = cms.string('pfTauTagInfoRef().pfjetRef().eta()'),
    jetEnergy = cms.string('pfTauTagInfoRef().pfjetRef().energy()'),

    # Jet width moments
    etaetaMoment = cms.string('etaetaMoment()'),
    phiphiMoment = cms.string('phiphiMoment()'),
    etaphiMoment = cms.string('etaphiMoment()'),

    # discriminators
    leadingTrackFinding = cms.string('tauID("leadingTrackFinding")'),
    leadingTrackPtCut = cms.string('tauID("leadingTrackPtCut")'),
    leadingPionPtCut = cms.string('tauID("leadingPionPtCut")'),
    trackIsolation = cms.string('tauID("trackIsolation")'),
    trackIsolationUsingLeadingPion = cms.string('tauID("trackIsolationUsingLeadingPion")'),
    ecalIsolation = cms.string('tauID("ecalIsolation")'),
    ecalIsolationUsingLeadingPion = cms.string('tauID("ecalIsolationUsingLeadingPion")'),
    byIsolation = cms.string('tauID("byIsolation")'),
    byIsolationUsingLeadingPion = cms.string('tauID("byIsolationUsingLeadingPion")'),
    againstElectron = cms.string('tauID("againstElectron")'),
    againstMuon = cms.string('tauID("againstMuon")'),
)
