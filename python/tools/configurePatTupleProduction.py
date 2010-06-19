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

def configurePatTupleProduction(process, patSequenceBuilder = None, 
                                patPFTauCleanerPrototype = None, 
                                patCaloTauCleanerPrototype = None,
                                addGenInfo = False):

    # check that patSequenceBuilder and patTauCleanerPrototype are defined and non-null
    if patSequenceBuilder is None:
        raise ValueError("Undefined patSequenceBuilder Parameter !!")
    if patPFTauCleanerPrototype is None or patCaloTauCleanerPrototype is None:
        raise ValueError("Undefined patTauCleanerPrototype Parameter !!")

    #--------------------------------------------------------------------------------
    # produce PAT objects
    #--------------------------------------------------------------------------------

    process.load("PhysicsTools.PatAlgos.patSequences_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cff")
    process.load("TauAnalysis.CandidateTools.muTauPairProduction_cff")

    if not addGenInfo:
        removeMCMatching(process)

    #--------------------------------------------------------------------------------
    # configure PAT trigger matching
    process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
    
    process.patTauTriggerMatchHLTsingleJet15UprotoType = cms.EDFilter("PATTriggerMatcherDRDPtLessByR",
        src                   = cms.InputTag("cleanLayer1Taus"),
        matched               = cms.InputTag("patTrigger"),
        andOr                 = cms.bool(False),
        filterIdsEnum         = cms.vstring('*'),
        filterIds             = cms.vint32(0),
        filterLabels          = cms.vstring('hlt1jet15U'),
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
    process.patCaloTauProducer = copy.deepcopy(process.patTaus)

    retVal_caloTau = patSequenceBuilder(
        process,
        collectionName = [ "patCaloTaus", "" ],
        patTauProducerPrototype = process.patCaloTauProducer,
        patTauCleanerPrototype = patCaloTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTsingleJet15UprotoType,
        addGenInfo = addGenInfo
    )
    process.caloTauSequence = retVal_caloTau["sequence"]

    process.patMuonCaloTauPairs = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuons'),
        srcLeg2 = cms.InputTag(retVal_caloTau["collection"]),
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
    process.patPFTauProducerFixedCone = copy.deepcopy(process.patTaus)

    retVal_pfTauFixedCone = patSequenceBuilder(
        process,
        collectionName = [ "patPFTaus", "FixedCone" ],
        patTauProducerPrototype = process.patPFTauProducerFixedCone,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTsingleJet15UprotoType,
        addGenInfo = addGenInfo
    )
    process.pfTauSequenceFixedCone = retVal_pfTauFixedCone["sequence"]

    process.patMuonPFTauPairsFixedCone = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuons'),
        srcLeg2 = cms.InputTag(retVal_pfTauFixedCone["collection"]),
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
    process.patPFTauProducerShrinkingCone = copy.deepcopy(process.patTaus)

    # add "transformed" TaNC output to pat::Tau
    ##setattr(process.patPFTauProducerShrinkingCone.tauIDSources,
    ##        "transformedTaNCoutput", cms.InputTag("shrinkingConePFTauTancCVTransform"))
        
    # load TaNC inputs into pat::Tau
    process.load("RecoTauTag.Configuration.ShrinkingConePFTaus_cfi")
    process.load("RecoTauTag.TauTagTools.PFTauMVAInputDiscriminatorTranslator_cfi")
    loadMVAInputsIntoPatTauDiscriminants(process.patPFTauProducerShrinkingCone)
    # enable embedding of decay mode from PFTauDecayMode data format
    process.patPFTauProducerShrinkingCone.addDecayMode = cms.bool(True)

    retVal_pfTauShrinkingCone = patSequenceBuilder(
        process,
        collectionName = [ "patPFTaus", "ShrinkingCone" ],
        patTauProducerPrototype = process.patPFTauProducerShrinkingCone,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTsingleJet15UprotoType,
        addGenInfo = addGenInfo
    )
    process.pfTauSequenceShrinkingCone = retVal_pfTauShrinkingCone["sequence"]

    process.patMuonPFTauPairsShrinkingCone = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuons'),
        srcLeg2 = cms.InputTag(retVal_pfTauShrinkingCone["collection"]),
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
    process.patPFTauProducerHPS = copy.deepcopy(process.patTaus)

    retVal_pfTauHPS = patSequenceBuilder(
        process,
        collectionName = [ "patPFTaus", "HPS" ],
        patTauProducerPrototype = process.patPFTauProducerHPS,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTsingleJet15UprotoType,
        addGenInfo = addGenInfo
    )
    process.pfTauSequenceHPS = retVal_pfTauHPS["sequence"]

    process.patMuonPFTauPairsHPS = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuons'),
        srcLeg2 = cms.InputTag(retVal_pfTauHPS["collection"]),
        srcMET = cms.InputTag('patPFMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # replace caloJets by pfJets
    switchJetCollection(process, jetCollection = cms.InputTag("iterativeCone5PFJets"))
    #
    # NOTE: need to delete empty sequence produced by call to "switchJetCollection"
    #       in order to avoid error when calling "process.dumpPython"
    #      ( cf. https://hypernews.cern.ch/HyperNews/CMS/get/physTools/1688/1/1/1/1/1.html )
    #
    del process.patJetMETCorrections
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
       # recompute decay modes and embed TaNC inputs
       + process.shrinkingConePFTauDecayModeProducer               
       + process.produceTancMVAInputDiscriminators
       + process.pfTauSequenceFixedCone + process.pfTauSequenceShrinkingCone + process.pfTauSequenceHPS
       + process.patMuonCaloTauPairs
       + process.patMuonPFTauPairsFixedCone + process.patMuonPFTauPairsShrinkingCone + process.patMuonPFTauPairsHPS
    )

    # return names of "final" collections of CaloTaus/different types of PFTaus
    # to be used as InputTag for further processing
    retVal = {}
    retVal["caloTauCollection"] = retVal_caloTau["collection"]
    retVal["pfTauCollectionFixedCone"] = retVal_pfTauFixedCone["collection"]
    retVal["pfTauCollectionShrinkingCone"] = retVal_pfTauShrinkingCone["collection"]
    retVal["pfTauCollectionHPS"] = retVal_pfTauHPS["collection"]
    return retVal
