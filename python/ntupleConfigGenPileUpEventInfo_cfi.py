import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables for multiplicity of mixed-in pile-up interactions
# added on Monte Carlo generator level 
#--------------------------------------------------------------------------------

genPileUpEventInfo_template = cms.PSet(
    vector = cms.bool(False),
    
    pluginType = cms.string("GenPileUpEventInfoExtractor"),

    src = cms.InputTag('addPileupInfo'),

    columns = cms.PSet(
        genNumPileUp = cms.string("numPileUpInteractions"),
    )
)    
