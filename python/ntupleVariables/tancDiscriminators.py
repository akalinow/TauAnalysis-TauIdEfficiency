import FWCore.ParameterSet.Config as cms

variables = cms.PSet(
    byTaNC = cms.string('tauID("byTaNC")'),
    byTaNCfrOnePercent = cms.string('tauID("byTaNCfrOnePercent")'),
    byTaNCfrHalfPercent = cms.string('tauID("byTaNCfrHalfPercent")'),
    byTaNCfrQuarterPercent = cms.string('tauID("byTaNCfrQuarterPercent")'),
    byTaNCfrTenthPercent = cms.string('tauID("byTaNCfrTenthPercent")'),
    decayMode = cms.string('decayMode()'),
)
