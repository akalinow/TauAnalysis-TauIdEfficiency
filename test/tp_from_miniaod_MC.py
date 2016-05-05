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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

import os
if "CMSSW_7_6_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('76X_mcRun2_asymptotic_v12')
    process.source.fileNames = [
        'file:///home/akalinow/scratch/CMS/TauID/Data/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM/0C4FAF25-DAC7-E511-A535-BCAEC54E98B3.root'
    ]
else: raise RuntimeError, "Unknown CMSSW version %s" % os.environ['CMSSW_VERSION']

dataPath = "/home/akalinow/scratch/CMS/TauID/Data/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"
command = "ls "+dataPath+"/*.root"
fileList = commands.getoutput(command).split("\n")   
process.source.fileNames =  cms.untracked.vstring()
for aFile in fileList:
    process.source.fileNames.append('file:'+aFile)


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
    process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_IsoMu18_v*')
    
else:
    raise RuntimeError, "TRIGGER must be 'SingleMu' or 'DoubleMu'"

process.triggerResultsFilter.l1tResults = "gtDigis"
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")
process.metSelector = cms.EDFilter("PATMETSelector",
    src = cms.InputTag("slimmedMETs"),
    cut = cms.string("pt < 25"),
    filter =  cms.bool(True),                                
)
process.fastFilter     = cms.Sequence(process.goodVertexFilter + process.triggerResultsFilter + process.metSelector)

##    __  __                       
##   |  \/  |_   _  ___  _ __  ___ 
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##                                 
from MuonAnalysis.TagAndProbe.common_variables_cff import *
from TauAnalysis.TauIdEfficiency.common_variables_tau_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 25 && abs(eta)<2.1 && "+ MuonIDFlags2016.Tight2016.value()+
                     " && " + MuonIDFlags2016.Isolation2016.value())
)



process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.tagTriggerMatchModule = cms.EDProducer("TriggerObjectStandAloneMatch", 
    tags   = cms.InputTag("tagMuons"),
    objects = cms.InputTag("selectedPatTrigger"),
    objectSelection = cms.string('hasFilterLabel("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09")'),
    maxTagObjDR   = cms.double(0.1),
)

process.mtSelector = cms.EDProducer("CandViewShallowCloneCombiner",
     decay = cms.string("slimmedMETs tagMuons"),
     checkCharge = cms.bool(False),                           
     cut   = cms.string("sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi)))<40"),
)
process.mtFilter = cms.EDFilter("CandViewCountFilter",
                           src = cms.InputTag("mtSelector"),
                           minNumber = cms.uint32(1)
                           )

##
## Taus
##
process.mergedTaus = cms.EDProducer("PFTauMerger",
    mergeTracks = cms.bool(True),
    taus     = cms.InputTag("slimmedTaus"), 
    tracks    = cms.InputTag("packedPFCandidates"),
    ## Apply some minimal pt cut
    tausCut     = cms.string("pt > 20 && abs(eta)<2.3"),
    tracksCut    = cms.string("pt > 20 && abs(eta)<2.3"),
)

process.probeTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("mergedTaus"),
    cut = cms.string(''),
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 120 && abs(daughter(0).vz - daughter(1).vz) < 4'),
    decay = cms.string('tagMuons@+ probeTaus@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

process.goodGenMuons.src = cms.InputTag("prunedGenParticles")

process.njets30Module.objects = cms.InputTag("slimmedJets")

process.pairMETPtBalanceModule = cms.EDProducer("CandMETPairPtBalance", 
    pairs   = cms.InputTag("tpPairs"),
    objects = cms.InputTag("slimmedMETs"),
    objectSelection = cms.string(""), 
    minTagObjDR   = cms.double(0.0),
    minProbeObjDR = cms.double(0.0),
)

process.tagMuonsMCMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("tagMuons"), # RECO objects to match
    matched = cms.InputTag("goodGenMuons"),   # mc-truth particle collection
    mcPdgId     = cms.vint32(13),  # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(False), # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(1),      # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.1),   # Minimum deltaR for the match
    maxDPtRel = cms.double(0.5),   # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True), # False = just match input in order; True = pick lowest deltaR pair first
)
process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(src = "probeTaus", maxDeltaR = 0.1, maxDPtRel = 1.0, resolveAmbiguities = False,  resolveByMatchQuality = False)

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    # probe variables: all useful ones
    variables = cms.PSet(
        KinematicVariables,        
    ),
    flags = cms.PSet(
        decayModeFinding =  cms.string('tauID("decayModeFinding")'),
        decayModeFindingNewDMs =  cms.string('tauID("decayModeFindingNewDMs")'),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.string('tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")'),
        againstMuonLoose3 =  cms.string('tauID("againstMuonLoose3")'),
        againstMuonTight3 =  cms.string('tauID("againstMuonTight3")')
    ),
    tagVariables = cms.PSet(    
        KinematicVariables,
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
        METPtBalance = cms.InputTag("pairMETPtBalanceModule"),
    ),
    pairFlags = cms.PSet(
        BestZ = cms.InputTag("bestPairByZMass"),
    ),
    isMC           = cms.bool(True),
    addRunLumiInfo = cms.bool(True),
    tagMatches       = cms.InputTag("tagMuonsMCMatch"),
    probeMatches     = cms.InputTag("probeMuonsMCMatch"),
    motherPdgId      = cms.vint32(23),
    makeMCUnbiasTree       = cms.bool(False), 
    checkMotherInUnbiasEff = cms.bool(True),
)

process.nverticesModule.objects = cms.InputTag("offlineSlimmedPrimaryVertices")

process.tnpSimpleSequence = cms.Sequence(
    process.goodGenMuons +
    process.tagMuons*process.tagMuonsMCMatch +  
    process.oneTag     +
    process.tagTriggerMatchModule + 
    process.mtSelector*process.mtFilter +
    process.probeTaus*process.probeMuonsMCMatch +
    process.tpPairs    +
    process.onePair    +
    process.nverticesModule +
    process.njets30Module +
    process.pairMETPtBalanceModule +
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

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ_MC.root"))
