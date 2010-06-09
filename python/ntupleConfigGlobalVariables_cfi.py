import FWCore.ParameterSet.Config as cms

extraVarforTauCands_template =  cms.PSet(
    vector = cms.bool(True),

    pluginType = cms.string("PATTauVectorExtraValExtractor"),

    src = cms.InputTag("patPFTausDijetTagAndProbeShrinkingCone"),
    pfCandSrc = cms.InputTag("pfCandidates"),
    jetSrc = cms.InputTag("ak5PFJets"),
    jetMinPt = cms.double(10.),
    jetMaxAbsEta = cms.double(2.5),
    
    columns = cms.PSet(
    nTracksOther = cms.string("nTracksOut"),
    nChargedHadrOther = cms.string("nChargedHadrOut"),
    nGammaOther = cms.string("nGammaOut"),
    
    closestJetDR = cms.string("ClosestJetDR"),
    closestJetPt = cms.string("ClosestJetPt"),
    closestJetEta = cms.string("ClosestJetEta"),
    closestJetPhi = cms.string("ClosestJetPhi"),
    closestJetJetWidth = cms.string("ClosestJetJetWidth"),
    )
)

jets_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("NumCandidateExtractor"),
    src = cms.InputTag("patJets"),

    columns = cms.PSet(
    nJets = cms.string("JetsMultiplicity")
    )  
)


met_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("PATMetValExtractor"),

    src = cms.InputTag("patPFMETs"),

    columns = cms.PSet(
    Met = cms.string("p4().Et()"),
    sumEt = cms.string("sumEt()")
    )

)
