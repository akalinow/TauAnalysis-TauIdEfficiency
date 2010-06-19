import FWCore.ParameterSet.Config as cms

extraTauCandVariables_template = cms.PSet(
    vector = cms.bool(True),

    pluginType = cms.string("PATTauVectorExtraValExtractor"),

    src = cms.InputTag("patPFTausDijetTagAndProbeShrinkingCone"),
    pfCandSrc = cms.InputTag("particleFlow"),
    jetSrc = cms.InputTag("ak5PFJets"),
    jetMinPt = cms.double(10.),
    jetMaxAbsEta = cms.double(2.5),
    
    columns = cms.PSet(
        numTracksOther = cms.string("numTracksOut"),
        numChargedHadrOther = cms.string("numChargedHadrOut"),
        numGammaOther = cms.string("numPhotonsOut"),
    
        nearestJetDR = cms.string("nearestJetDR"),
        nearestJetPt = cms.string("nearestJetPt"),
        nearestJetEta = cms.string("nearestJetEta"),
        nearestJetPhi = cms.string("nearestJetPhi"),
        nearestJetWidth = cms.string("nearestJetWidth"),
    )
)

jets_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("NumCandidateExtractor"),
    src = cms.InputTag("patJets"),

    columns = cms.PSet(
        numJets = cms.string("JetsMultiplicity")
    )  
)

met_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("PATMetValExtractor"),

    src = cms.InputTag("patPFMETs"),

    columns = cms.PSet(
        MEt = cms.string("p4().Et()"),
        sumEt = cms.string("sumEt()")
    )
)
