import FWCore.ParameterSet.Config as cms

#process = cms.Process("testNtuple")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
#process.load('Configuration/StandardSequences/Services_cff')
#process.load('FWCore/MessageService/MessageLogger_cfi')
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration/StandardSequences/GeometryIdeal_cff')
#process.load('Configuration/StandardSequences/MagneticField_cff')
#process.load('Configuration/StandardSequences/Reconstruction_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'MC_31X_V3::All'

#--------------------------------------------------------------------------------
# import sequence for PAT-tuple production
#process.load("TauAnalysis.Configuration.producePatTuple_cff")

# import event-content definition of products to be stored in patTuple
#from TauAnalysis.Configuration.patTupleEventContent_cff import *
#from TauAnalysis.Skimming.EventContent_cff import *
#--------------------------------------------------------------------------------
from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.load("PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi")

# print event content 
#process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Comment out if not at Davis T3
        'file:/mnt/hadoop/store/mc/Summer09/Ztautau/AODSIM/MC_31X_V3_AODSIM-v1/0008/08206CD7-7F7F-DE11-9C51-000AE488BBBA.root',
        # Uncomment to run on lxplus/lcp
        #'/store/relval/CMSSW_3_3_6/RelValZTT/GEN-SIM-RECO/STARTUP3X_V8H-v1/0009/C884F463-40E4-DE11-B272-003048679220.root',
        #'/store/relval/CMSSW_3_3_6/RelValZTT/GEN-SIM-RECO/STARTUP3X_V8H-v1/0009/B2466310-42E4-DE11-88CD-0030486792B6.root',
        #'/store/relval/CMSSW_3_3_6/RelValZTT/GEN-SIM-RECO/STARTUP3X_V8H-v1/0009/40201F20-41E4-DE11-9F5A-003048678B0E.root',
        #'/store/relval/CMSSW_3_3_6/RelValZTT/GEN-SIM-RECO/STARTUP3X_V8H-v1/0009/3AD94763-40E4-DE11-B871-002618FDA287.root',
        #'/store/relval/CMSSW_3_3_6/RelValZTT/GEN-SIM-RECO/STARTUP3X_V8H-v1/0009/24052BB1-9EE4-DE11-87C2-002618943957.root'
    ),
    skipEvents = cms.untracked.uint32(0)            
)

#--------------------------------------------------------------------------------

# produce ntuple
process.ntupleProducer = cms.EDAnalyzer(
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
            src = cms.InputTag("allLayer1Taus"),

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

# Save ntuple
process.out = cms.OutputModule(
    "PoolOutputModule",                                                                                                                                                        
    outputCommands = cms.untracked.vstring("drop *", "keep *_*ntupleProducer*_*_*" ),
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string("example_ntuple.root")      
)

process.p = cms.Path(
    process.makeAllLayer1Taus
    * process.selectedLayer1Taus
    * process.ntupleProducer )

process.end = cms.EndPath(process.out)
