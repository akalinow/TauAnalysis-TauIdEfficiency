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

numGlobalMuons_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("NumCandidateExtractor"),
    src = cms.InputTag("patMuonsGlobal"),

    columns = cms.PSet(
        multiplicity = cms.string("numGlobalMuons")
    )  
)

numStandAloneMuons_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("NumCandidateExtractor"),
    src = cms.InputTag("patMuonsStandAlone"),

    columns = cms.PSet(
        multiplicity = cms.string("numStandAloneMuons")
    )  
)

jets_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("NumSelPATJetExtractor"),
    src = cms.InputTag("patJets"),
    srcNotToBeFiltered = cms.VInputTag(
        "cleanPatElectrons",
        "cleanPatMuons",
    ),
    dRmin = cms.double(0.5),

    columns = cms.PSet(
        numJetsEta25Pt10 = cms.string("abs(eta) < 2.5 & pt > 10."),
        numJetsEta25Pt15 = cms.string("abs(eta) < 2.5 & pt > 15."),
        numJetsEta25Pt20 = cms.string("abs(eta) < 2.5 & pt > 20."),
        numJetsEta25Pt25 = cms.string("abs(eta) < 2.5 & pt > 25."),
        numJetsEta25Pt30 = cms.string("abs(eta) < 2.5 & pt > 30.")
    )  
)

caloMet_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("PATMetValExtractor"),

    src = cms.InputTag("patMETs"),

    columns = cms.PSet(
        caloMEt = cms.string("p4().Et()"),
        caloSumEt = cms.string("sumEt()")
    )
)

pfMet_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("PATMetValExtractor"),

    src = cms.InputTag("patPFMETs"),

    columns = cms.PSet(
        pfMEt = cms.string("p4().Et()"),
        pfSumEt = cms.string("sumEt()")
    )
)

diTau_template = cms.PSet(
    vector = cms.bool(False),

    pluginType = cms.string("DiCandidatePairValExtractor"),

    src = cms.InputTag(""),

    columns = cms.PSet()
)
