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
        #'file:/mnt/hadoop/store/mc/Summer09/Ztautau/AODSIM/MC_31X_V3_AODSIM-v1/0008/08206CD7-7F7F-DE11-9C51-000AE488BBBA.root',
        # Uncomment to run on lxplus/lcp
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/FAB3991D-8039-DF11-8E2E-002618FDA277.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/DC514536-9639-DF11-B076-001BFCDBD160.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/C4855F52-9439-DF11-A454-0018F3D0962E.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/B6ABAECD-8E39-DF11-A23B-00261894389E.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/B287B900-9739-DF11-AB6D-0018F3D095EC.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/9238AEC8-8F39-DF11-94CF-003048678FE4.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/8CC8BCCB-9439-DF11-8572-0018F3D0968A.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/7A13CC45-9139-DF11-883B-00304867915A.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/56EEA948-9039-DF11-9782-0030486790FE.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/4C4DE0CB-9539-DF11-B29F-001BFCDBD154.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/441C2548-9539-DF11-9E04-001BFCDBD160.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/427DC648-8F39-DF11-8B87-002618943821.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/1CDE1F22-533A-DF11-A12A-00261894388A.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/08B28CC8-8F39-DF11-B0F8-0026189438D7.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0007/F6AC0AB4-1138-DF11-A8F8-00261894383A.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0006/C884454F-D337-DF11-8BCA-002618943935.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0006/70362C4B-D337-DF11-AF72-00304867BED8.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0006/50717D2B-D237-DF11-84BD-00304867BED8.root',
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0006/2AE4C62F-D837-DF11-8718-002354EF3BDB.root'
    ),
    skipEvents = cms.untracked.uint32(0)            
)

#--------------------------------------------------------------------------------

# dijet tag and probe production test
process.load("TauAnalysis.TauIdEfficiency.TauIdDijetTagAndProbeProducer_cfi")
process.dijetTagAndProbes.source = cms.InputTag("selectedPatTaus")

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
            src = cms.InputTag("dijetTagAndProbes", "highestPtProbe"),

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
    process.makePatTaus
    * process.selectedPatTaus
    * process.dijetTagAndProbes
    * process.ntupleProducer )

process.end = cms.EndPath(process.out)
