import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables for HLT trigger paths
#--------------------------------------------------------------------------------

trigger_template = cms.PSet(
    vector = cms.bool(False),
    
    pluginType = cms.string("PATTriggerInfoExtractor"),

    src = cms.InputTag('patTriggerEvent'),

    columns = cms.PSet(
        hltL1Jet6Ubit               = cms.string("HLT_L1Jet6U:bit"),
        hltL1Jet6Uprescale          = cms.string("HLT_L1Jet6U:prescale"),
        hltJet15Ubit                = cms.string("HLT_Jet15U:bit"),
        hltJet15Uprescale           = cms.string("HLT_Jet15U:prescale"),
        hltJet30Ubit                = cms.string("HLT_Jet30U:bit"),
        hltJet30Uprescale           = cms.string("HLT_Jet30U:prescale"),
        hltJet50Ubit                = cms.string("HLT_Jet50U:bit"),
        hltJet50Uprescale           = cms.string("HLT_Jet50U:prescale"),
        hltMinBiasBSCbit            = cms.string("HLT_MinBiasBSC:bit"),
        hltMinBiasBSCprescale       = cms.string("HLT_MinBiasBSC:prescale"),
        hltMinBiasBSCnoBPTXbit      = cms.string("HLT_MinBiasBSC_NoBPTX:bit"),
        hltMinBiasBSCnoBPTXprescale = cms.string("HLT_MinBiasBSC_NoBPTX:prescale")
    ),

    maxWarnings = cms.int32(1)
)    
