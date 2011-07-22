import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    #fileNames   = cms.vstring('fitTauIdEff_wConstraints_2011June30_matthew.root'),
    fileNames   = cms.vstring('/data1/veelken/tmp/muonPtGt20/V4d/analyzeTauIdEffHistograms_all_2011Jul06_mauroV4.root'),
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('/data1/veelken/tmp/muonPtGt20/V4b/fitTauIdEff.root')
)

process.fitTauIdEff = cms.PSet(

    # CV: set to '' if determining tau id. efficiency for whole Tag & Probe sample,
    #     set to name of one individual bin in case you want to measure the tau id. efficiency as function of tauPt, tauEta,...
    #    (needs as many 'fitTauIdEff' jobs to be run in parallel as there are bins)
    directory = cms.string(''),

    runClosureTest = cms.bool(False),
    #runClosureTest = cms.bool(True),

    regions = cms.vstring(
        'A1',  # needed to access QCD template obtained from Data
        'B1',  # needed to access QCD template obtained from Data
        'B1p',
        'B1f',
        'C1',
        'C1p',
        'C1f',
        'D1',  # needed to measure tau charge misidentification rate
        'D1p',
        'D1f'
    ),

    # regions (in Data) from which templates for QCD background are taken
    regionTakeQCDtemplateFromData_passed = cms.string('B1p'),
    regionTakeQCDtemplateFromData_failed = cms.string('B1f'),
    #takeQCDfromData = cms.bool(False),
    takeQCDfromData = cms.bool(True),

    # define "passed" and "failed" regions
    passed_region = cms.string('C1p'),
    failed_region = cms.string('C1f'), # use for tau id. efficiency measurement
    #failed_region = cms.string('D1p'), # use for measurement of tau charge misidentification efficiency
    
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
        #diTauSVfitMass
    ),

    allowTemplateMorphing = cms.bool(True),
    morphSysUncertainty = cms.string(
        "sysTauJetEn"
    ),
    sysVariedByNsigma = cms.double(3.0),
    
    loadSysUncertainties = cms.vstring(
        "sysTauJetEn"
    ),

    runPseudoExperiments = cms.bool(False),
    #runPseudoExperiments = cms.bool(True),
    numPseudoExperiments = cms.uint32(10000),
    varySysUncertainties = cms.vstring(
        "sysTauJetEn"
    ),

    #makeControlPlots = cms.bool(True)
    makeControlPlots = cms.bool(False)
)
