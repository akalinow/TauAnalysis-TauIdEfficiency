import FWCore.ParameterSet.Config as cms

# Common values taken from the reco::Candidate data type

commonVariables = cms.PSet(
    pt = cms.string('pt()'),
    mass = cms.string('mass()'),
    phi = cms.string('phi()'),
    eta = cms.string('eta()'),
    absEta = cms.string('abs(eta())'),
    energy = cms.string('energy()'),
)

