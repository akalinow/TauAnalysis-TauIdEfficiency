import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    #fileNames   = cms.vstring('fitTauIdEff_wConstraints_2011June30_matthew.root'),
    fileNames   = cms.vstring('/data1/veelken/tmp/analyzeTauIdEffHistograms_all_2011Jul01_mauroV2.root'),
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('/data1/veelken/tmp/fitTauIdEff.root')
)

process.fitTauIdEff = cms.PSet(

    # CV: set to '' if determining tau id. efficiency for whole Tag & Probe sample,
    #     set to name of one individual bin in case you want to measure the tau id. efficiency as function of tauPt, tauEta,...
    #    (needs as many 'fitTauIdEff' jobs to be run in parallel as there are bins)
    directory = cms.string(''),

    #runClosureTest = cms.bool(False),
    runClosureTest = cms.bool(True),

    takeQCDfromData = cms.bool(False),
    #takeQCDfromData = cms.bool(True),
    
    runSysUncertainties = cms.bool(False),
    #runSysUncertainties = cms.bool(True),

    numPseudoExperiments = cms.uint32(10000),

    regions = cms.vstring(
        'B1',  # needed to access QCD template obtained from Data
        'C1',
        'C1p',
        'C1f'
    ),
    
    tauIds = cms.vstring(
        'tauDiscrHPSloose', # "new" HPS implemented in HPS+TaNC combined algorithm
        'tauDiscrHPSmedium',
        'tauDiscrHPStight'
    ),

    fitVariables = cms.vstring(
        'diTauVisMass',
        #'diTauVisMassFromJet' # CV: diTauVisMass always computed from PFJet momenta if using PAT-tuple workflow
    ),

    sysUncertainties = cms.vstring(
        "CENTRAL_VALUE",
        #"SysTauJetEn", # needed for diTauVisMass/diTauVisMassFromJet
        #"SysJetEnUp"   # needed for diTauMt
    )
)