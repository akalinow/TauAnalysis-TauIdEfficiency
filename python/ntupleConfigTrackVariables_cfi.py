import FWCore.ParameterSet.Config as cms

tauTrackVariables_template = cms.PSet(
    vector = cms.bool(True),

    pluginType = cms.string("PATTauVectorTrackValExtractor"),

    src = cms.InputTag("patPFTausDijetTagAndProbeShrinkingCone"),
    vertexSrc = cms.InputTag("offlinePrimaryVertices"),

    collection = cms.string("leadTrack"),
    
    columns = cms.PSet(
        numValidHitsLeadTrack = cms.string("numValidHits"),
        numMissingHitsLeadTrack = cms.string("numMissingHits"),
        chi2LeadTrack = cms.string("chi2"),
        nDoFLeadTrack = cms.string("nDoF"),
        dzLeadTrack = cms.string("dz"),
        dxyLeadTrack = cms.string("dxy"),
        ptLeadTrack = cms.string("pt"),
        ptErrLeadTrack = cms.string("ptErr"),
        qualityLeadTrack = cms.string("quality")
    )
)
