import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('produceTauPtResPATTuple')

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.MessageLogger.suppressInfo = cms.untracked.vstring()
process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V11C::All')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ##'file:/data1/veelken/CMSSW_5_2_x/skims/selEvents_debugPATTaus_forSimon_AOD.root'
        'file:/data1/veelken/CMSSW_5_2_x/skims/selEvents_debugPATTaus_forSimon_discrAgainstMuonsFailed_AOD.root'                               
    ),
    #eventsToProcess = cms.untracked.VEventRange(
    #    '1:5:4262773',
    #    '1:5:4141629'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load('RecoTauTag.Configuration.RecoPFTauTag_cff')
## import RecoTauTag.RecoTau.RecoTauPiZeroBuilderPlugins_cfi as builders
## modStrips = copy.deepcopy(builders.strips)
## modStrips.plugin = cms.string('RecoTauPiZeroStripPlugin2')
## #modStrips.qualityCuts.signalQualityCuts.minGammaEt = cms.double(0.)
## modStrips.minGammaEtStripSeed = cms.double(0.5)
## modStrips.minGammaEtStripAdd = cms.double(0.)
## modStrips.minStripEt = cms.double(1.0)
## modStrips.updateStripAfterEachDaughter = cms.bool(False)
## modStrips.maxStripBuildIterations = cms.int32(-1)
## process.ak5PFJetsLegacyHPSPiZeros.builders = cms.VPSet(
##     modStrips
## )

process.load('PhysicsTools.PatAlgos.patSequences_cff')

# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    doJTA = True,
    doBTagging = True,
    jetCorrLabel = ( 'AK5PF', cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute') ),
    doType1MET = False,
    doJetID = True,
    jetIdLabel = "ak5",
    outputModules = []
)

# configure pat::Tau production
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
process.patTaus.userIsolation = cms.PSet()
process.patTaus.isoDeposits = cms.PSet()

## import PhysicsTools.PatAlgos.tools.helpers as configtools
## process.allTauPtResSequences = cms.Sequence()

process.tauPtResSequence = cms.Sequence(
    process.recoTauClassicHPSSequence
   + process.tauMatch + process.tauGenJetMatch + process.patTaus
)

## tauPtResOptions = {
##     'Default' : {
##         'applyElecTrackQcuts'          : True,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.5,
##         'minStripEt'                   : 0.5,
##         'updateStripAfterEachDaughter' : True,
##         'maxStripBuildIterations'      :  1
##     },
##     'Seed05Add00Strip05' : {
##         'applyElecTrackQcuts'          : True,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.0,
##         'minStripEt'                   : 0.5,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsDefault' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.5,
##         'minStripEt'                   : 0.5,
##         'updateStripAfterEachDaughter' : True,
##         'maxStripBuildIterations'      :  1
##     },
##     'noEleTrackQcutsSeed00Add00Strip05' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.0,
##         'minGammaEtStripAdd'           : 0.0,
##         'minStripEt'                   : 0.5,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsSeed00Add00Strip10' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.0,
##         'minGammaEtStripAdd'           : 0.0,
##         'minStripEt'                   : 1.0,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsSeed05Add00Strip05' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.0,
##         'minStripEt'                   : 0.5,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsSeed05Add00Strip10' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.0,
##         'minStripEt'                   : 1.0,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsSeed05Add00Strip15' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.0,
##         'minStripEt'                   : 1.5,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsSeed05Add00Strip20' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.0,
##         'minStripEt'                   : 2.0,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsSeed05Add00Strip25' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.0,
##         'minStripEt'                   : 2.5,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsSeed05Add05Strip05' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.5,
##         'minStripEt'                   : 0.5,
##         'updateStripAfterEachDaughter' : False,
##         'maxStripBuildIterations'      : -1
##     },
##     'noEleTrackQcutsSeed05Add05Strip10' : {
##         'applyElecTrackQcuts'          : False,
##         'minGammaEtStripSeed'          : 0.5,
##         'minGammaEtStripAdd'           : 0.5,
##         'minStripEt'                   : 1.0,
##         'updateStripAfterEachDaughter' : True,
##         'maxStripBuildIterations'      : -1
##     }
## }

## for tauPtResName, tauPtResConfig in tauPtResOptions.items():

##     modStrips.applyElecTrackQcuts = cms.bool(tauPtResConfig['applyElecTrackQcuts'])
##     modStrips.minGammaEtStripSeed = cms.double(tauPtResConfig['minGammaEtStripSeed'])
##     modStrips.minGammaEtStripAdd = cms.double(tauPtResConfig['minGammaEtStripAdd'])
##     modStrips.minStripEt = cms.double(tauPtResConfig['minStripEt'])
##     modStrips.updateStripAfterEachDaughter = cms.bool(tauPtResConfig['updateStripAfterEachDaughter'])
##     modStrips.maxStripBuildIterations = cms.int32(tauPtResConfig['maxStripBuildIterations'])
    
##     configtools.cloneProcessingSnippet(process, process.tauPtResSequence, tauPtResName)

##     process.allTauPtResSequences *= getattr(process, "%s%s" % ("tauPtResSequence", tauPtResName))

## process.ak5PFJetsLegacyHPSPiZeros.builders = cms.VPSet(
##     builders.strips
## )

# before starting to process 1st event, print event content
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")
process.filterFirstEvent = cms.EDFilter("EventCountFilter",
    numEvents = cms.int32(1)
)

##process.o = cms.Path(process.filterFirstEvent + process.printEventContent)

process.p = cms.Path(
    process.PFTau
   + process.patDefaultSequence
   ##+ process.allTauPtResSequences
)

process.patTupleOutputModule = cms.OutputModule("PoolOutputModule",
    cms.PSet(
        outputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_patTaus*_*_*',
            'keep *_cleanPatMuons_*_*',
            'keep *_genParticles_*_*',
            'keep *_tauGenJets_*_*',
            'keep *_tauGenJetsSelectorAllHadrons_*_*',                          
            'keep *_offlinePrimaryVertices_*_*',
            'keep *_offlinePrimaryVerticesWithBS_*_*',
            'keep *_ak5PFJets_*_*',
            'keep *_particleFlow_*_*',
            'keep *_generalTracks_*_*',
            'keep *_electronGsfTracks_*_*',
            'keep *_towerMaker_*_*'
        )               
    ),
    fileName = cms.untracked.string("/data1/veelken/CMSSW_5_2_x/skims/selEvents_debugPATTaus_forSimon_tauPtResPATtuple.root")
)

process.q = cms.EndPath(process.patTupleOutputModule)

processDumpFile = open('produceTauPtResPATTuple.dump' , 'w')
print >> processDumpFile, process.dumpPython()
