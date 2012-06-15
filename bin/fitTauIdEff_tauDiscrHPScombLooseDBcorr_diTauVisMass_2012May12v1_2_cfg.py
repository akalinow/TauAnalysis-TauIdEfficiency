
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

#histogram from FWlite
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('/afs/cern.ch/work/c/calpas/CMSSW/FWliteHisto/analyzeTauIdEffHistograms_all_corrQCDtemplates_2012May12v1_2.root')
)

#output of the fit template histogram   
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('/afs/cern.ch/work/c/calpas/CMSSW/FWliteHisto/Fit/fitTauIdEff_tauDiscrHPScombLooseDBcorr_diTauVisMass_2012May12v1_2.root')
)

#pass all the cfg parameter 
process.fitTauIdEff = cms.PSet(

#histogram to smooth
    HistogramsToFit = cms.vstring(
'WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake', #LG1
'WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake' #CB1
'WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake', #CB2 
'WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake', #CB2 
'QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all', #CB3 ??
'QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all', #CB3 ??
),


#function to use for the smoothing // no necessary, the fit depends on the input histograme name
    FitFunctions = cms.vstring(
'LandauConvGauss_v1' #LG1
'CristalBall_v1' #CB1
'CristalBall_v2',#CB2
'CristalBall_v3' #CB3
),

#the name of the output histogram fitted // no necessary depend on the input histo name
    HistoFitted = cms.vstring(
'HistoFitted.root'
),


#the name of the control fit plot // no necessary depend on the input histo name
    ControlPlots = cms.vstring(
'controlPlots.png',
'controlPlots.pdf',
'controlPlots.eps'
),


    directory = cms.string(''),

    runClosureTest = cms.bool(False),
    

    regions = cms.vstring(
'A','B','D','C1f',
   ),


    # regions (Data or MC) from which templates for QCD background are taken
    regionQCDtemplate_passed = cms.string('C1p_qcd'),
    regionQCDtemplate_failed = cms.string('C1f_qcd'),
    regionQCDtemplate_D      = cms.string('D_qcd'),

    # define "passed" and "failed" regions
    #region_passed = cms.string('C1p'),
    #region_failed = cms.string('C1f'),
    #region_D = cms.string('D'),

    region_A = cms.string('A'),
    region_B = cms.string('B'),
    region_D = cms.string('D'),
    region_C1f = cms.string('C1f'),
   


    tauId = cms.vstring('tauDiscrHPScombLooseDBcorr'),

    fitVariable = cms.vstring(
'diTauVisMass','diTauMt'
),


    fitIndividualProcesses = cms.bool(False),


    templateMorphingMode = cms.string('vertical'), 

    sysUncertainties = cms.vstring(
'JetEn','TauJetEn','TauJetRes','UnclusteredEn'
    ),
    sysVariedByNsigma = cms.double(3.0),

    runPseudoExperiments = cms.bool(False),
    #runPseudoExperiments = cms.bool(True),
    numPseudoExperiments = cms.uint32(10000),

    intLumiData = cms.double(0.560700),

    makeControlPlots = cms.bool(True),
 
    controlPlotFilePath = cms.string('/afs/cern.ch/work/c/calpas/muonPtGt17/v1_2_fitEWKbgSum_v7/2012RunA/tauIdEfficiency/plots')

  )
