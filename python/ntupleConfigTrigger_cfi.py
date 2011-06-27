import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables for HLT trigger paths
#--------------------------------------------------------------------------------

trigger_template = cms.PSet(
    vector = cms.bool(False),
    
    pluginType = cms.string("PATTriggerInfoExtractor"),

    src = cms.InputTag('patTriggerEvent'),

    columns = cms.PSet(
        hltJet30v1bit         = cms.string("HLT_Jet30_v1:bit"),
        hltJet30v1prescale    = cms.string("HLT_Jet30_v1:prescale"),                                    
        hltJet60v1bit         = cms.string("HLT_Jet60_v1:bit"),
        hltJet60v1prescale    = cms.string("HLT_Jet60_v1:prescale"),
        hltJet80v1bit         = cms.string("HLT_Jet80_v1:bit"),
        hltJet80v1prescale    = cms.string("HLT_Jet80_v1:prescale"),                                    
        hltJet110v1bit        = cms.string("HLT_Jet110_v1:bit"),
        hltJet110v1prescale   = cms.string("HLT_Jet110_v1:prescale"),
        hltJet150v1bit        = cms.string("HLT_Jet150_v1:bit"),
        hltJet150v1prescale   = cms.string("HLT_Jet150_v1:prescale"),
        hltJet190v1bit        = cms.string("HLT_Jet190_v1:bit"),
        hltJet190v1prescale   = cms.string("HLT_Jet190_v1:prescale"),
        hltJet240v1bit        = cms.string("HLT_Jet240_v1:bit"),
        hltJet240v1prescale   = cms.string("HLT_Jet240_v1:prescale"),
        hltMu8v1bit           = cms.string("HLT_Mu8_v1:bit"),
        hltMu8v1prescale      = cms.string("HLT_Mu8_v1:prescale"),
        hltMu12v1bit          = cms.string("HLT_Mu12_v1:bit"),
        hltMu12v1prescale     = cms.string("HLT_Mu12_v1:prescale"),
        hltIsoMu12v1bit       = cms.string("HLT_IsoMu12_v1:bit"),
        hltIsoMu12v1prescale  = cms.string("HLT_IsoMu12_v1:prescale"),
        hltMu15v1bit          = cms.string("HLT_Mu15_v1:bit"),
        hltMu15v1prescale     = cms.string("HLT_Mu15_v1:prescale"),
        hltMu15v2bit          = cms.string("HLT_Mu15_v2:bit"),
        hltMu15v2prescale     = cms.string("HLT_Mu15_v2:prescale"),
        hltIsoMu17v5bit       = cms.string("HLT_IsoMu17_v5:bit"),
        hltIsoMu17v5prescale  = cms.string("HLT_IsoMu17_v5:prescale"),
        hltIsoMu17v6bit       = cms.string("HLT_IsoMu17_v6:bit"),
        hltIsoMu17v6prescale  = cms.string("HLT_IsoMu17_v6:prescale"),
        hltIsoMu17v8bit       = cms.string("HLT_IsoMu17_v8:bit"),
        hltIsoMu17v8prescale  = cms.string("HLT_IsoMu17_v8:prescale"),
        hltIsoMu17v9bit       = cms.string("HLT_IsoMu17_v9:bit"),
        hltIsoMu17v9prescale  = cms.string("HLT_IsoMu17_v9:prescale"),
        hltIsoMu17v11bit      = cms.string("HLT_IsoMu17_v11:bit"),
        hltIsoMu17v11prescale = cms.string("HLT_IsoMu17_v11:prescale"),
        hltMu24v1bit          = cms.string("HLT_Mu24_v1:bit"),
        hltMu24v1prescale     = cms.string("HLT_Mu24_v1:prescale"),
        hltMu24v2bit          = cms.string("HLT_Mu24_v2:bit"),
        hltMu24v2prescale     = cms.string("HLT_Mu24_v2:prescale")
    ),

    # define L1 seeds of all HLT trigger paths used in (edm)Ntuple filling
    # NOTE: definition will become obsolete once associated between HLT paths and L1 seeds (algorithms)
    #       is implemented in pat::TriggerEvent
    l1Seeds = cms.PSet(
        HLT_Jet30_v1    = cms.string("L1_SingleJet16"),
        HLT_Jet60_v1    = cms.string("L1_SingleJet36"),
        HLT_Jet80_v1    = cms.string("L1_SingleJet52"),
        HLT_Jet110_v1   = cms.string("L1_SingleJet68"),
        HLT_Jet150_v1   = cms.string("L1_SingleJet92"),
        HLT_Jet190_v1   = cms.string("L1_SingleJet92"),
        HLT_Jet240_v1   = cms.string("L1_SingleJet92"),
        HLT_Mu8_v1      = cms.string("L1_SingleMu3"),
        HLT_Mu12_v1     = cms.string("L1_SingleMu7"),
        HLT_IsoMu12_v1  = cms.string("L1_SingleMu7"),
        HLT_Mu15_v1     = cms.string("L1_SingleMu10"),
        HLT_Mu15_v2     = cms.string("L1_SingleMu10"),
        HLT_IsoMu15_v5  = cms.string("L1_SingleMu10"),
        HLT_IsoMu17_v5  = cms.string("L1_SingleMu10"),
        HLT_IsoMu17_v6  = cms.string("L1_SingleMu10"),
        HLT_IsoMu17_v8  = cms.string("L1_SingleMu10"),
        HLT_IsoMu17_v9  = cms.string("L1_SingleMu10"),
        HLT_IsoMu17_v11 = cms.string("L1_SingleMu10"),
        HLT_Mu24_v1     = cms.string("L1_SingleMu12"),
        HLT_Mu24_v2     = cms.string("L1_SingleMu12")
    ),

    maxWarnings = cms.int32(1)
)    
