import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    #fileNames   = cms.vstring('fitTauIdEff_wConstraints_2011June30_matthew.root'),
    fileNames   = cms.vstring('/data1/veelken/tmp/muonPtGt20/V6/analyzeTauIdEffHistograms_all_2011Jul23V6.root'),
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('/data1/veelken/tmp/muonPtGt20/V6/fitTauIdEff_wConstraints.root')
)

process.fitTauIdEff_wConstraints = cms.PSet(

    # CV: set to '' if determining tau id. efficiency for whole Tag & Probe sample,
    #     set to name of one individual bin in case you want to measure the tau id. efficiency as function of tauPt, tauEta,...
    #    (needs as many 'fitTauIdEff' jobs to be run in parallel as there are bins)
    directory = cms.string(''),

    runClosureTest = cms.bool(False),
    #runClosureTest = cms.bool(True),

    # regions (in Data) from which templates for QCD background are taken
    regionTakeQCDtemplateFromData_all    = cms.string('B1'),
    regionTakeQCDtemplateFromData_passed = cms.string('B1p'),
    regionTakeQCDtemplateFromData_failed = cms.string('B1f'),
    #takeQCDfromData = cms.bool(False),
    takeQCDfromData = cms.bool(True),

    # CV: fitting fake-rates of background processes
    #     in C2f/C2p regions causes bias of fit result (2011/06/28)
    fitTauIdEffC2 = cms.bool(False),
    #fitTauIdEffC2 = cms.bool(True),

    regions = cms.vstring(
        'ABCD',
        'A',
        'A1',  # QCD enriched control region (OS, loose muon isolation, Mt && Pzeta cuts applied)
        'B',
        'B1',  # QCD enriched control region (SS, loose muon isolation, Mt && Pzeta cuts applied)
        'B1p',
        'B1f',
        'C',
        'C1',
        'C1p',
        'C1f',
        'C2',
        'C2p',
        'C2f',
        'D',   # generic background control region (SS, tight muon isolation)
        #'D1',
        #'D1p',
        #'D1f',
        #'D2',
        #'D2p',
        #'D2f'
    ),
    
    tauIds = cms.vstring(
        ##'tauDiscrHPSloose', # "new" HPS implemented in HPS+TaNC combined algorithm
        ##'tauDiscrHPSlooseDBcorr',
        'tauDiscrHPScombLooseDBcorr',
        ##'tauDiscrHPSmedium',
        ##'tauDiscrHPSmediumDBcorr',
        'tauDiscrHPScombMediumDBcorr',
        ##'tauDiscrHPStight',
        ##'tauDiscrHPStightDBcorr',
        'tauDiscrHPScombTightDBcorr'
    ),

    fitVariables = cms.vstring(
        'diTauVisMass',
        #'diTauVisMassFromJet' # CV: diTauVisMass always computed from PFJet momenta if using PAT-tuple workflow
        #'diTauSVfitMass'
    ),

    allowTemplateMorphing = cms.bool(False),
    morphSysUncertainty = cms.string(
        "sysTauJetEn"
    ),
    sysVariedByNsigma = cms.double(3.0),
    morphQCDinABD = cms.bool(True),
    
    loadSysUncertainties = cms.vstring(
        "sysTauJetEn",
        "sysJetEn",
        "sysAddPUsmearing"
    ),

    runPseudoExperiments = cms.bool(False),
    #runPseudoExperiments = cms.bool(True),
    numPseudoExperiments = cms.uint32(10000),
    varySysUncertainties = cms.vstring(
        "sysTauJetEn"
    ),

    makeControlPlots = cms.bool(True)
    #makeControlPlots = cms.bool(False)
)
