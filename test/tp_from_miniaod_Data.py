import FWCore.ParameterSet.Config as cms

import subprocess

process = cms.Process("TagProbe")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

import os
if "CMSSW_7_6_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15')
    process.source.fileNames = [
        'file:///home/akalinow/scratch/CMS/TauID/Data/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD/00006301-CAA8-E511-AD39-549F35AD8BC9.root'
    ]

else: raise RuntimeError, "Unknown CMSSW version %s" % os.environ['CMSSW_VERSION']

## SELECT WHAT DATASET YOU'RE RUNNING ON
TRIGGER="SingleMu"

## ==== Fast Filters ====
process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
    filter = cms.bool(True),
)

process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")

if TRIGGER == "SingleMu":
    process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_IsoMu17_eta2p1_v*', 'HLT_IsoMu18_v*',
                                                                  'HLT_IsoMu20_v*')
    
else:
    raise RuntimeError, "TRIGGER must be 'SingleMu' or 'DoubleMu'"

process.triggerResultsFilter.l1tResults = "gtDigis"
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")
process.fastFilter     = cms.Sequence(process.goodVertexFilter + process.triggerResultsFilter)

##    __  __                       
##   |  \/  |_   _  ___  _ __  ___ 
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##                                 
from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 15 && "+MuonIDFlags.Tight2012.value()+
                     " && pfIsolationR04().sumChargedHadronPt/pt < 0.2"),
)

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

##
## Taus
##
process.mergedTaus = cms.EDProducer("PFTauMerger",
    mergeTracks = cms.bool(True),
    taus     = cms.InputTag("slimmedTaus"), 
    tracks    = cms.InputTag("packedPFCandidates"),
    ## Apply some minimal pt cut
    tausCut     = cms.string("pt > 10"),
    tracksCut    = cms.string("pt > 10"),
)

process.probeTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("mergedTaus"),
    cut = cms.string(''),
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    #cut = cms.string('60 < mass < 140 && abs(daughter(0).vz - daughter(1).vz) < 4'),
    cut = cms.string('60 < mass && abs(daughter(0).vz - daughter(1).vz) < 4'),
    decay = cms.string('tagMuons@+ probeTaus@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    # probe variables: all useful ones
    variables = cms.PSet(
        KinematicVariables,
        decayModeFindingNewDMs =  cms.string('tauID("decayModeFindingNewDMs")'),
        againstMuonLoose3 =  cms.string('tauID("againstMuonLoose3")'),
        againstMuonTight3 =  cms.string('tauID("againstMuonTight3")')
    ),
    flags = cms.PSet(
    ),
    tagVariables = cms.PSet(    
        KinematicVariables,
    ),
    tagFlags = cms.PSet(HighPtTriggerFlags,HighPtTriggerFlagsDebug),
    pairVariables = cms.PSet(
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
        probeMultiplicity = cms.InputTag("probeMultiplicity"),       
    ),
    pairFlags = cms.PSet(
        BestZ = cms.InputTag("bestPairByZMass"),
    ),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(False),
)

process.nverticesModule.objects = cms.InputTag("offlineSlimmedPrimaryVertices")

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeTaus +
    process.tpPairs    +
    process.onePair    +
    process.nverticesModule +
    process.probeMultiplicities + 
    process.bestPairByZMass + 
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.fastFilter +
    process.mergedTaus +                 
    process.tnpSimpleSequence
)


process.schedule = cms.Schedule(
   process.tagAndProbe
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ_Data.root"))
