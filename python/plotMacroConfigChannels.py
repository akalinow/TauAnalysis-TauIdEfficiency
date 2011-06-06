import FWCore.ParameterSet.Config as cms

#------------------------------------------------------------------------------------------------------------------------------------
#            Collection of samples used by the fitting {name, totalEvevts,[Xsec:MC],,inputFileName,label,histoColor (fill), dataset}
#------------------------------------------------------------------------------------------------------------------------------------
import TauAnalysis.Configuration.dbsInterface as dbsInterface
import TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi as samples

data = cms.PSet( #data_SingleMu_Run2011A_PromptReco_v1,v2
    name = cms.string("data"), #keyword! Do not change
    totalEvevts = cms.double( dbsInterface.GetNumEvents(samples.RECO_SAMPLES['data_SingleMu_Run2011A_PromptReco_v2']['datasetpath']) + samples.RECO_SAMPLES['data_SingleMu_Run2011A_PromptReco_v1']['datasetpath']) ),
    inputFileName = cms.string(""),
    label =  cms.string("Data"),
    histoColor = cms.int32(1)
    )

Ztautau_powhegZ2 = cms.PSet( #data_SingleMu_Run2011A_PromptReco_v1,v2
    name = cms.string("Ztautau"), #keyword! Do not change
    totalEvevts = cms.double( dbsInterface.GetNumEvents(samples.RECO_SAMPLES['Ztautau_powhegZ2']['datasetpath']) ),
    inputFileName = cms.string(""),
    label =  cms.string("Z/#gamma^{*} #rightarrow #tau^{+} #tau^{-}"),
    histoColor = cms.int32(628)
    )

Zmumu_pythia = cms.PSet( #data_SingleMu_Run2011A_PromptReco_v1,v2
    name = cms.string("data"), #keyword! Do not change
    totalEvevts = cms.double( dbsInterface.GetNumEvents(samples.RECO_SAMPLES['Zmumu_pythia']['datasetpath']) ),
    inputFileName = cms.string(""),
    label =  cms.string("Z/#gamma^{*} #rightarrow #mu^{+} #mu^{-}"),
    histoColor = cms.int32(596)
    )

PPmuXptGt20Mu15 = cms.PSet( #data_SingleMu_Run2011A_PromptReco_v1,v2
    name = cms.string("data"), #keyword! Do not change
    totalEvevts = cms.double( dbsInterface.GetNumEvents(samples.RECO_SAMPLES['PPmuXptGt20Mu15']['datasetpath']) ),
    inputFileName = cms.string(""),
    label =  cms.string("QCD"),
    histoColor = cms.int32(797)
    )

WplusJets_madgraph = cms.PSet( #data_SingleMu_Run2011A_PromptReco_v1,v2
    name = cms.string("data"), #keyword! Do not change
    totalEvevts = cms.double( dbsInterface.GetNumEvents(samples.RECO_SAMPLES['WplusJets_madgraph']['datasetpath']) ),
    inputFileName = cms.string(""),
    label =  cms.string("W + jets"),
    histoColor = cms.int32(856)
    )

TTplusJets_madgraph = cms.PSet( #data_SingleMu_Run2011A_PromptReco_v1,v2
    name = cms.string("data"), #keyword! Do not change
    totalEvevts = cms.double( dbsInterface.GetNumEvents(samples.RECO_SAMPLES['TTplusJets_madgraph']['datasetpath']) ),
    inputFileName = cms.string(""),
    label =  cms.string("t#bar{t} + jets"),
    histoColor = cms.int32(618)
    )



