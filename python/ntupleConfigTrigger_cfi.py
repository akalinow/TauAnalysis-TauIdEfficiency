import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables for HLT trigger paths
#--------------------------------------------------------------------------------

trigger_template = cms.PSet(
    vector = cms.bool(False),
    
    pluginType = cms.string("PATTriggerInfoExtractor"),

    src = cms.InputTag('patTriggerEvent'),

    columns = cms.PSet(
        hltJet15Ubit       = cms.string("HLT_Jet15U:bit"),
        hltJet15Uprescale  = cms.string("HLT_Jet15U:prescale"),                                    
        hltJet30Ubit       = cms.string("HLT_Jet30U:bit"),
        hltJet30Uprescale  = cms.string("HLT_Jet30U:prescale"),                                    
        hltJet50Ubit       = cms.string("HLT_Jet50U:bit"),
        hltJet50Uprescale  = cms.string("HLT_Jet50U:prescale"),
        hltMu5bit          = cms.string("HLT_Mu5:bit"),
        hltMu5prescale     = cms.string("HLT_Mu5:prescale"),                                    
        hltMu9bit          = cms.string("HLT_Mu9:bit"),
        hltMu9prescale     = cms.string("HLT_Mu9:prescale"),
        hltIsoMu9bit       = cms.string("HLT_IsoMu9:bit"),
        hltIsoMu9prescale  = cms.string("HLT_IsoMu9:prescale"),                                    
        hltMu11bit         = cms.string("HLT_Mu11:bit"),
        hltMu11prescale    = cms.string("HLT_Mu11:prescale"),      
    ),

    # define L1 seeds of all HLT trigger paths used in (edm)Ntuple filling
    # NOTE: definition will become obsolete once associated between HLT paths and L1 seeds (algorithms)
    #       is implemented in pat::TriggerEvent
    l1Seeds = cms.PSet(
        HLT_Jet15U = cms.string("L1_SingleJet6U"),
        HLT_Jet30U = cms.string("L1_SingleJet20U"),
        HLT_Jet50U = cms.string("L1_SingleJet30U"),
        HLT_Mu5 = cms.string("L1_SingleMu3"),
        HLT_Mu9 = cms.string("L1_SingleMu7"),
        HLT_IsoMu9 = cms.string("L1_SingleMu7"),
        HLT_Mu11 = cms.string("L1_SingleMu7")
    ),

    maxWarnings = cms.int32(1)
)    
