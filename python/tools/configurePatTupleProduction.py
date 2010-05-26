import FWCore.ParameterSet.Config as cms
import copy

from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching

from TauAnalysis.Configuration.tools.metTools import *

# Get the files to support embedding of TaNC inputs
from RecoTauTag.TauTagTools.PFTauMVAInputDiscriminatorTranslator_cfi import \
        loadMVAInputsIntoPatTauDiscriminants

from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildDijetTauSequence

def configurePatTupleProduction(process, addGenInfo = False):

    #--------------------------------------------------------------------------------
    # produce PAT objects
    #--------------------------------------------------------------------------------

    process.load("PhysicsTools.PatAlgos.patSequences_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cff")
    process.load("TauAnalysis.CandidateTools.muTauPairProduction_cff")

    process.load("TauAnalysis.TauIdEfficiency.patConfiguration.matchingPrototypes")
    patCaloTauMatchProtoType = copy.deepcopy(process.dijetCleanerPrototype)
    patCaloTauMatchProtoType.checkOverlaps.TagJet.src = cms.InputTag("caloJetsTagAndProbes", "tagObject")
    patCaloTauMatchProtoType.checkOverlaps.HighestPtProbeJet.src = cms.InputTag("caloJetsTagAndProbes", "highestPtProbe")
    patCaloTauMatchProtoType.checkOverlaps.SecondHighestPtProbe.src = cms.InputTag("caloJetsTagAndProbes", "secondHighestPtProbe")

    patPFTauMatchProtoType = copy.deepcopy(process.dijetCleanerPrototype)
    patPFTauMatchProtoType.checkOverlaps.TagJet.src = cms.InputTag("pfJetsTagAndProbes", "tagObject")
    patPFTauMatchProtoType.checkOverlaps.HighestPtProbeJet.src = cms.InputTag("pfJetsTagAndProbes", "highestPtProbe")
    patPFTauMatchProtoType.checkOverlaps.SecondHighestPtProbe.src = cms.InputTag("pfJetsTagAndProbes", "secondHighestPtProbe")

    if not addGenInfo:
        removeMCMatching(process)

    #--------------------------------------------------------------------------------
    # configure PAT trigger matching
    process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
    
    patTauTriggerMatchHLTsingleJet15UprotoType = cms.EDFilter("PATTriggerMatcherDRDPtLessByR",
        src                   = cms.InputTag("cleanLayer1Taus"),
        matched               = cms.InputTag("patTrigger"),
        andOr                 = cms.bool(False),
        filterIdsEnum         = cms.vstring('*'),
        filterIds             = cms.vint32(0),
        filterLabels          = cms.vstring('*'),
        pathNames             = cms.vstring('HLT_Jet15U'),
        collectionTags        = cms.vstring('*'),
        maxDPtRel             = cms.double(0.5),
        maxDeltaR             = cms.double(0.5),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(False)
    )                                                             
    #--------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------- 
    #
    # produce combinations of muon + tau-jet pairs
    # for collection of pat::Tau objects representing CaloTaus 
    #
    switchToCaloTau(process)    
    patCaloTauProducer = copy.deepcopy(process.patTaus)

    process.caloTauSequence = buildDijetTauSequence(
        process,
        collectionName = [ "patCaloTaus", "" ],
        patTauProducerPrototype = patCaloTauProducer,
        triggerMatcherProtoType = patTauTriggerMatchHLTsingleJet15UprotoType,
        addGenInfo = addGenInfo
    )

    process.patMuonCaloTauPairs = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuons'),
        srcLeg2 = cms.InputTag('patCaloTausDijetTagAndProbe'),
        srcMET = cms.InputTag('patMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collection of pat::Tau objects representing PFTaus
    # reconstructed by fixed signal cone algorithm
    # (plus combinations of muon + tau-jet pairs)
    #
    switchToPFTauFixedCone(process)
    patPFTauProducerFixedCone = copy.deepcopy(process.patTaus)

    process.pfTauSequenceFixedCone = buildDijetTauSequence(
        process,
        collectionName = [ "patPFTaus", "FixedCone" ],
        patTauProducerPrototype = patPFTauProducerFixedCone,
        triggerMatcherProtoType = patTauTriggerMatchHLTsingleJet15UprotoType,
        addGenInfo = addGenInfo
    )

    process.patMuonPFTauPairsFixedCone = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuons'),
        srcLeg2 = cms.InputTag('patPFTausDijetTagAndProbeFixedCone'),
        srcMET = cms.InputTag('patPFMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collection of pat::Tau objects representing PFTaus
    # reconstructed by shrinking signal cone algorithm
    # (plus combinations of muon + tau-jet pairs) 
    #
    switchToPFTauShrinkingCone(process)
    patPFTauProducerShrinkingCone = copy.deepcopy(process.patTaus)

    # Load TaNC inputs into pat::Tau
    process.load("RecoTauTag.Configuration.ShrinkingConePFTaus_cfi")
    process.load("RecoTauTag.TauTagTools.PFTauMVAInputDiscriminatorTranslator_cfi")
    loadMVAInputsIntoPatTauDiscriminants(patPFTauProducerShrinkingCone)
    # Enable embedding of decay mode from PFT Decay Mode data format
    patPFTauProducerShrinkingCone.addDecayMode = cms.bool(True)

    process.pfTauSequenceShrinkingCone = buildDijetTauSequence(
        process,
        collectionName = [ "patPFTaus", "ShrinkingCone" ],
        patTauProducerPrototype = patPFTauProducerShrinkingCone,
        triggerMatcherProtoType = patTauTriggerMatchHLTsingleJet15UprotoType,
        addGenInfo = addGenInfo
    )

    process.patMuonPFTauPairsShrinkingCone = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuons'),
        srcLeg2 = cms.InputTag('patPFTausDijetTagAndProbeShrinkingCone'),
        srcMET = cms.InputTag('patPFMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collection of pat::Tau objects representing PFTaus
    # reconstructed by hadron + strips (HPS) algorithm
    # (plus combinations of muon + tau-jet pairs) 
    #
    switchToPFTauHPS(process)
    patPFTauProducerHPS = copy.deepcopy(process.patTaus)

    process.pfTauSequenceHPS = buildDijetTauSequence(
        process,
        collectionName = [ "patPFTaus", "HPS" ],
        patTauProducerPrototype = patPFTauProducerHPS,
        triggerMatcherProtoType = patTauTriggerMatchHLTsingleJet15UprotoType,
        addGenInfo = addGenInfo
    )

    process.patMuonPFTauPairsHPS = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuons'),
        srcLeg2 = cms.InputTag('patPFTausDijetTagAndProbeHPS'),
        srcMET = cms.InputTag('patPFMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # replace caloJets by pfJets
    switchJetCollection(process, jetCollection = cms.InputTag("iterativeCone5PFJets"))
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # add pfMET
    # set Boolean swich to true in order to apply type-1 corrections
    addPFMet(process, correct = False)
    #--------------------------------------------------------------------------------

    process.patTupleProductionSequence = cms.Sequence(
        process.patDefaultSequence
       + process.patTrigger
       + process.caloTauSequence
       # Recomputed decay modes and embed TaNC inputs
       + process.shrinkingConePFTauDecayModeProducer               
       + process.produceTancMVAInputDiscriminators
       + process.pfTauSequenceFixedCone + process.pfTauSequenceShrinkingCone + process.pfTauSequenceHPS
       # EKF temporarily turning off, breaks runtime
       #+ process.patMuonCaloTauPairs
       #+ process.patMuonPFTauPairsFixedCone + process.patMuonPFTauPairsShrinkingCone + process.patMuonPFTauPairsHPS
    )
