import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables for HLT trigger paths
#--------------------------------------------------------------------------------

trigger_template = cms.PSet(
    vector = cms.bool(False),
    
    pluginType = cms.string("HLTInfoExtractor"),

    src = cms.InputTag("TriggerResults::HLT"),

    columns = cms.PSet(
        hltL1Jet6U         = cms.string("HLT_L1Jet6U"),
        hltJet15U          = cms.string("HLT_Jet15U"),
        hltJet30U          = cms.string("HLT_Jet30U"),
        hltJet50U          = cms.string("HLT_Jet50U"),
        hltMinBiasBSC      = cms.string("HLT_MinBiasBSC"),
        hltMinBiasBSCNoBPTX= cms.string("HLT_MinBiasBSC_NoBPTX")
    )
)    
