
import FWCore.ParameterSet.Config as cms

ntupleProducer = cms.EDAnalyzer(
    "ObjValEDNtupleProducer",
    ntupleName = cms.string("exampleNtuple"),
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc
        hadronicTaus = cms.PSet(
            # Select multiplicy of object(s) to store
            vector = cms.bool(True), # Store a value for all objects in this collection
            #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects

            # Extractor plugin
            pluginType = cms.string("PATTauVectorValExtractor"),

            # Collection to extract from
            src = cms.InputTag("selectedPatTaus"),

            # Variables to compute for this source
            columns = cms.PSet(
                absEta = cms.string("abs(eta())"),
                pt = cms.string("pt()"),
                byLeadPionPt = cms.string('tauID("leadingPionPtCut")'), #NB quote format!
                byIsolation = cms.string('tauID("byIsolation")'),
                byTaNCfrOne = cms.string('tauID("byTaNCfrOnePercent")'),
            )
        ),
    )
)

