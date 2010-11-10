import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to Generator level jets, including
# Tau GenJets - GenJets containing the visible products of true tau decays
#--------------------------------------------------------------------------------

tauGenJets_genInfo = cms.PSet(
    pluginType = cms.string("GenJetVectorValExtractor"),
    vector=cms.bool(True),

    src = cms.InputTag("tauGenJets"),

    columns = cms.PSet(
        # generator level information
        genPt = cms.string("genPt"),
        genEta = cms.string("genEta"),
        genPhi = cms.string("genPhi"),
        
        genMass = cms.string("genMass"),
        
        genDecayMode = cms.string("genDecayMode")
    )
)    

genJets_genInfo = cms.PSet(
    pluginType = cms.string("GenJetVectorValExtractor"),
    vector=cms.bool(True),

    src = cms.InputTag("tauGenJets"),

    columns = cms.PSet(
        # generator level information
        genPt = cms.string("genPt"),
        genEta = cms.string("genEta"),
        genPhi = cms.string("genPhi")
    )
)    

