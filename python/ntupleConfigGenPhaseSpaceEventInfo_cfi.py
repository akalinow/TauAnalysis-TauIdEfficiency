import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables for Monte Carlo generator information
# (e.g. PtHat value in events produced by PYTHIA)
#--------------------------------------------------------------------------------

genPhaseSpaceEventInfo_template = cms.PSet(
    vector = cms.bool(False),
    
    pluginType = cms.string("GenPhaseSpaceEventInfoExtractor"),

    src = cms.InputTag('generator'),

    columns = cms.PSet(
        genPtHat = cms.string("ptHat"),
    )
)    
