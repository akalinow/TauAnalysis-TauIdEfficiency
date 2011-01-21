import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define run number, luminosity section and event number variables
#--------------------------------------------------------------------------------

runLumiSectionEventNumber_template = cms.PSet(
    vector = cms.bool(False),
    
    pluginType = cms.string("RunLumiSectionEventNumberExtractor"),

    src = cms.InputTag(''),

    columns = cms.PSet(
        run   = cms.string("run"),
        ls    = cms.string("ls"),
        event = cms.string("event")
    )
)    
