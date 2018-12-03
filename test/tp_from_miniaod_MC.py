import FWCore.ParameterSet.Config as cms

import commands
import subprocess

#Checkout DNN training files if those are not present in the working area
#(dut to size the DNN training files should not be send with input sandbox)
import os

DNN_data_path = os.environ["CMSSW_BASE"]+"/external/"+os.environ["SCRAM_ARCH"]+"/data/RecoTauTag/TrainingFiles/data"
isDNN_present = os.path.exists(DNN_data_path)

print "isDNN_present:",isDNN_present

if not isDNN_present:
    command = "cd $CMSSW_BASE/src; "
    command += "git clone https://github.com/cms-tau-pog/RecoTauTag-TrainingFiles -b nonQuantizedDNN RecoTauTag/TrainingFiles/data; "
    command += "cd -;"
    os.system(command)

process = cms.Process("TagProbe")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

import os
if "CMSSW_8_0_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v7')
elif "CMSSW_9_2_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('92X_upgrade2017_realistic_v10')
elif "CMSSW_9_4_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v15')
else: raise RuntimeError, "Unknown CMSSW version %s" % os.environ['CMSSW_VERSION']


dataPath = "/home/akalinow/scratch/CMS/TauID/Data/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/"
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
    process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_IsoMu27_v*')    
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

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 28 && abs(eta)<2.1 && "+ MuonIDFlags2016.Tight2016.value()+
                     " && " + MuonIDFlags2016.Isolation2016.value())
)

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.unpackedPatTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
     patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
     triggerResults              = cms.InputTag('TriggerResults::HLT'),
     unpackFilterLabels = cms.bool(True)
)


process.tagTriggerMatchModule = cms.EDProducer("TriggerObjectStandAloneMatch", 
    tags   = cms.InputTag("tagMuons"),
    objects = cms.InputTag("unpackedPatTrigger"),
    objectSelection = cms.string('hasFilterLabel("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")'),
    maxTagObjDR   = cms.double(0.1),
)

##
## Taus
##
process.mergedTaus = cms.EDProducer("PFTauMerger",
    mergeTracks = cms.bool(True),
    taus     = cms.InputTag("slimmedTausNewMVAIDs"),
    photons     = cms.InputTag("slimmedPhotons"), 
    tracks    = cms.InputTag("packedPFCandidates"),
    ## Apply some minimal pt cut
    tausCut     = cms.string("pt > 20 && abs(eta)<2.3"),
    tracksCut    = cms.string("abs(eta)<2.3"),
)

process.probeTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("mergedTaus"),
    cut = cms.string('tauID("byTightIsolationMVArun2v1DBnewDMwLT2017v2")==1'),
    minNumber = cms.uint32(1)
)


process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('40 < mass < 200 && abs(daughter(0).vz - daughter(1).vz) < 4'),
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

process.pairAlternativeMass = cms.EDProducer("CandPairVariables", 
    pairs   = cms.InputTag("tpPairs"),
    genParticles = cms.InputTag("prunedGenParticles"),
    variableName =  cms.string("alternativeMass")
)

process.ZDecayMode = cms.EDProducer("CandPairVariables", 
    pairs   = cms.InputTag("tpPairs"),
    genParticles = cms.InputTag("prunedGenParticles"),
    variableName =  cms.string("ZDecayMode")
)

process.tagMuonsMCMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("tagMuons"), # RECO objects to match
    matched = cms.InputTag("goodGenMuons"),   # mc-truth particle collection
    mcPdgId     = cms.vint32(13),  # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(1),      # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.1),   # Minimum deltaR for the match
    maxDPtRel = cms.double(1),   # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True), # False = just match input in order; True = pick lowest deltaR pair first
)
process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(src = "probeTaus", maxDeltaR = 0.1, maxDPtRel = 1, resolveAmbiguities = False,  resolveByMatchQuality = False)

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
        byTightIsolationMVArun2v1DBnewDMwLT2017v2 = cms.string('tauID("byTightIsolationMVArun2v1DBnewDMwLT2017v2")'),
        againstMuonLoose3 =  cms.string('tauID("againstMuonLoose3")'),
        againstMuonTight3 =  cms.string('tauID("againstMuonTight3")'),
    ),
    tagVariables = cms.PSet(    
        KinematicVariables,
        dB = cms.string("dB"),
        triggerMatch = cms.InputTag("tagTriggerMatchModule"),
        nVertices   = cms.InputTag("nverticesModule"),
        isolation = MuonIDFlags2016.IsolationValue2016
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
        nVertices   = cms.InputTag("nverticesModule"),
        ## Gen related variables
        ZDecayMode = cms.InputTag("ZDecayMode"),        
        genWeight    = cms.InputTag("genAdditionalInfo", "genWeight"),
        truePileUp   = cms.InputTag("genAdditionalInfo", "truePileUp"),
        actualPileUp = cms.InputTag("genAdditionalInfo", "actualPileUp"),
        ###
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
process.genAdditionalInfo.pileUpInfoTag = cms.InputTag("slimmedAddPileupInfo")

###Rerun tauID
from tauIdRerun import *
addMVA_WPs_run2_2017(process)
##############

process.tnpSimpleSequence = cms.Sequence(
    process.goodGenMuons +
    process.tagMuons*process.tagMuonsMCMatch +  
    process.oneTag     +
    process.unpackedPatTrigger +
    process.tagTriggerMatchModule + 
    process.probeTaus*process.probeMuonsMCMatch +
    process.tpPairs    +
    process.onePair    +
    process.nverticesModule +
    process.njets30Module +
    process.pairMETPtBalanceModule + process.tagMTModule +  process.probeMTModule +
    process.pairAlternativeMass + process.ZDecayMode +
    process.probeMultiplicities + 
    process.bestPairByZMass + 
    process.genAdditionalInfo +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.fastFilter +
    process.newTauMVAIDsSeq + 
    process.mergedTaus +                 
    process.tnpSimpleSequence
)


process.schedule = cms.Schedule(
    process.tagAndProbe
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ_MC.root"))
