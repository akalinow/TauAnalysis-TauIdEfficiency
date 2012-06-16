
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('/data1/veelken/tmp/muonPtGt17/v1_4_fitEWKbgSum_v1_wFall11.bak/2012RunAplusB/tauIdEfficiency/analyzeTauIdEffHistograms_all_corrQCDtemplates_2012May12v1_4.root')
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('debug.root')
)

process.smoothTauIdEffTemplates = cms.PSet(
    directory = cms.string(''),
    
    histogramsToSmooth = cms.VPSet(
        cms.PSet(
            histogramName = cms.string('WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake'),
            fitFunctionType = cms.string("LG1")
        ),
        cms.PSet(
            histogramName = cms.string('WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake'),
            fitFunctionType = cms.string("CB1")
        ),
        cms.PSet(
            histogramName = cms.string('WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake'),
            fitFunctionType = cms.string("CB2")
        ),
        cms.PSet(
            histogramName = cms.string('WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake'),
            fitFunctionType = cms.string("CB2")
        ),
        cms.PSet(
            histogramName = cms.string('QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all'),
            fitFunctionType = cms.string("CB3")
        ),
        cms.PSet(
            histogramName = cms.string('QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all'),
            fitFunctionType = cms.string("CB3")
        )           
    ),

    makeControlPlots = cms.bool(True),
    controlPlotFilePath = cms.string('/data1/veelken/tmp/muonPtGt17/v1_4_fitEWKbgSum_v1_wFall11.bak/2012RunAplusB/tauIdEfficiency/plots'),
)
