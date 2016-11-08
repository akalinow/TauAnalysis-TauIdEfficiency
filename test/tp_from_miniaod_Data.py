import FWCore.ParameterSet.Config as cms

import commands
import subprocess

process = cms.Process("TagProbe")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

import os
if "CMSSW_7_6_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15')
    process.source.fileNames = [
        'file:///home/akalinow/scratch/CMS/TauID/Data/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD/1/00006301-CAA8-E511-AD39-549F35AD8BC9.root',
    ]
elif "CMSSW_8_0_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('80X_dataRun2_Prompt_ICHEP16JEC_v0')
    process.source.fileNames = [
        'file:///home/akalinow/scratch/CMS/TauID/Data/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD/1/00006301-CAA8-E511-AD39-549F35AD8BC9.root',
    ]
else: raise RuntimeError, "Unknown CMSSW version %s" % os.environ['CMSSW_VERSION']

'''
dataPath = "/scratch_local/akalinow/CMS/TauID/Data/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD"
command = "ls "+dataPath+"/*.root"
fileList = commands.getoutput(command).split("\n")   
process.source.fileNames =  cms.untracked.vstring()
for aFile in fileList:
    process.source.fileNames.append('file:'+aFile)
'''

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
    process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_IsoMu18_v*', 'HLT_IsoMu20_v*', 'HLT_IsoMu22_v*)
    
else:
    raise RuntimeError, "TRIGGER must be 'SingleMu' or 'DoubleMu'"

process.triggerResultsFilter.l1tResults = "gtStage2Digis"
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")

process.muonFilter = cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("slimmedMuons"),
                                  minNumber = cms.uint32(1),
                                  maxNumber = cms.uint32(2))

process.fastFilter     = cms.Sequence(process.goodVertexFilter + process.triggerResultsFilter +  process.muonFilter)

##    __  __                       
##   |  \/  |_   _  ___  _ __  ___ 
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##                                 
from MuonAnalysis.TagAndProbe.common_variables_cff import *
from TauAnalysis.TauIdEfficiency.common_variables_tau_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

'''
process.tagMuons = cms.EDProducer("PFMuonMerger",
    mergeTracks = cms.bool(False),
    muons     = cms.InputTag("slimmedMuons"),
    photons     = cms.InputTag("slimmedPhotons"), 
    tracks    = cms.InputTag("packedPFCandidates"),
    ## Apply some minimal pt cut
    muonsCut     = cms.string("pt > 25 && abs(eta)<2.1 && "+ MuonIDFlags2016.Tight2016.value()+
                              " && " + MuonIDFlags2016.Isolation2016.value()),
    tracksCut    = cms.string("pt > 25 && abs(eta)<2.1"),
)
'''

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 25 && abs(eta)<2.1 && "+ MuonIDFlags2016.Tight2016.value()+
                     " && " + MuonIDFlags2016.Isolation2016.value())
)

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.tagTriggerMatchModule = cms.EDProducer("TriggerObjectStandAloneMatch", 
    tags   = cms.InputTag("tagMuons"),
    objects = cms.InputTag("selectedPatTrigger"),
    objectSelection = cms.string('hasFilterLabel("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09")'),
    maxTagObjDR   = cms.double(0.1),
)

##
## Taus
##
process.mergedTaus = cms.EDProducer("PFTauMerger",
    mergeTracks = cms.bool(True),
    taus     = cms.InputTag("slimmedTaus"),
    photons     = cms.InputTag("slimmedPhotons"), 
    tracks    = cms.InputTag("packedPFCandidates"),
    ## Apply some minimal pt cut
    tausCut     = cms.string("pt > 20 && abs(eta)<2.3"),
    tracksCut    = cms.string("abs(eta)<2.3"),
)

process.probeTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("mergedTaus"),
    cut = cms.string('tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")==1'),
    minNumber = cms.uint32(1)
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('40 < mass < 160 && abs(daughter(0).vz - daughter(1).vz) < 4'),
    decay = cms.string('tagMuons@+ probeTaus@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

process.njets30Module.objects = cms.InputTag("slimmedJets")

process.pairMETPtBalanceModule = cms.EDProducer("CandMETPairPtBalance", 
    pairs   = cms.InputTag("tpPairs"),
    objects = cms.InputTag("slimmedMETs"),
    objectSelection = cms.string(""), 
    minTagObjDR   = cms.double(0.0),
    minProbeObjDR = cms.double(0.0),
)

process.tagMTModule = cms.EDProducer("CandMTPair", 
    pairs   = cms.InputTag("tpPairs"),
    objects = cms.InputTag("slimmedMETs"),
    objectSelection = cms.string(""), 
    minTagObjDR   = cms.double(0.0),
    minProbeObjDR = cms.double(0.0),
    useProbe = cms.bool(False)                              
)

process.probeMTModule = cms.EDProducer("CandMTPair", 
    pairs   = cms.InputTag("tpPairs"),
    objects = cms.InputTag("slimmedMETs"),
    objectSelection = cms.string(""), 
    minTagObjDR   = cms.double(0.0),
    minProbeObjDR = cms.double(0.0),
    useProbe = cms.bool(True)                              
)

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    # probe variables: all useful ones
    variables = cms.PSet(
        KinematicVariables,
        decayMode =  cms.string('tauID("decayMode")'),
        alternatLorentzVectPt =  cms.string("alternatLorentzVect.pt"),
        alternatLorentzVectEta =  cms.string("alternatLorentzVect.eta"),
    ),
    flags = cms.PSet(
        decayModeFinding =  cms.string('tauID("decayModeFinding")'),
        decayModeFindingNewDMs =  cms.string('tauID("decayModeFindingNewDMs")'),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.string('tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")'),
        byTightIsolationMVArun2v1DBoldDMwLT = cms.string('tauID("byTightIsolationMVArun2v1DBoldDMwLT")'),
        againstMuonLoose3 =  cms.string('tauID("againstMuonLoose3")'),
        againstMuonTight3 =  cms.string('tauID("againstMuonTight3")')
    ),
    tagVariables = cms.PSet(    
        KinematicVariables,
        dB = cms.string("dB"),
        triggerMatch = cms.InputTag("tagTriggerMatchModule")
    ),
    tagFlags = cms.PSet(),
    pairVariables = cms.PSet(
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
        probeMultiplicity = cms.InputTag("probeMultiplicity"),
        nJets30 = cms.InputTag("njets30Module"),
        MET = cms.InputTag("pairMETPtBalanceModule"),
        MTprobe = cms.InputTag("probeMTModule"),
        MTtag = cms.InputTag("tagMTModule"),
        alternativeMass = cms.InputTag("pairAlternativeMass"),
    ),
    pairFlags = cms.PSet(
        BestZ = cms.InputTag("bestPairByZMass"),
    ),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

process.nverticesModule.objects = cms.InputTag("offlineSlimmedPrimaryVertices")

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.tagTriggerMatchModule + 
    process.probeTaus +
    process.tpPairs    +
    process.onePair    +
    process.nverticesModule +
    process.njets30Module +     
    process.pairMETPtBalanceModule + process.tagMTModule +  process.probeMTModule +
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
