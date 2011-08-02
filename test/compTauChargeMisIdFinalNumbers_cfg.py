import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring('/data1/veelken/tmp/muonPtGt20/V6/compTauChargeMisIdFinalNumbers_input.root')
)
    
process.fwliteOutput = cms.PSet(
    fileName = cms.string('/data1/veelken/tmp/muonPtGt20/V6/compTauChargeMisIdFinalNumbers_output.root')
)

process.compTauChargeMisIdFinalNumbers = cms.PSet(

    # CV: set to '' if determining tau id. efficiency for whole Tag & Probe sample,
    #     set to name of one individual bin in case you want to measure the tau id. efficiency as function of tauPt, tauEta,...
    #    (needs as many 'fitTauIdEff' jobs to be run in parallel as there are bins)
    directory = cms.string(''),

    passed_region = cms.string('C1p'),
    failed_region = cms.string('D1p'),
        
    tauIds = cms.vstring(
        'tauDiscrHPSloose', # "new" HPS implemented in HPS+TaNC combined algorithm
        'tauDiscrHPSlooseDBcorr',
        'tauDiscrHPScombLooseDBcorr',
        'tauDiscrHPSmedium',
        'tauDiscrHPSmediumDBcorr',
        'tauDiscrHPScombMediumDBcorr',
        'tauDiscrHPStight',
        'tauDiscrHPStightDBcorr',
        'tauDiscrHPScombTightDBcorr'
    ),

    fitVariables = cms.vstring(
        'diTauVisMass',
        #'diTauVisMassFromJet' # CV: diTauVisMass always computed from PFJet momenta if using PAT-tuple workflow
    )
)
