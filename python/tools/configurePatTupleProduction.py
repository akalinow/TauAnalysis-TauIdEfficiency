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

from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildGenericTauSequence

def configurePatTupleProduction(process, patSequenceBuilder = buildGenericTauSequence, 
                                patPFTauCleanerPrototype = None, 
                                patCaloTauCleanerPrototype = None,
                                addSVfitInfo = False,
                                hltProcess = "HLT",
                                isMC = False,
                                applyTauVertexMatch = True):

    # check that patSequenceBuilder and patTauCleanerPrototype are defined and non-null
    if patSequenceBuilder is None:
        raise ValueError("Undefined 'patSequenceBuilder' Parameter !!")
    if patPFTauCleanerPrototype is None or patCaloTauCleanerPrototype is None:
        raise ValueError("Undefined 'patTauCleanerPrototype' Parameter !!")

    #--------------------------------------------------------------------------------
    # produce PAT objects
    #--------------------------------------------------------------------------------

    process.load("PhysicsTools.PatAlgos.patSequences_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff")
    process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cff")
    process.load("TauAnalysis.CandidateTools.muTauPairProduction_cff")

    # per default, do **not** run SVfit algorithm
    if not addSVfitInfo:
        process.allMuTauPairs.doSVreco = cms.bool(False)
        process.allMuTauPairs.doPFMEtSign = cms.bool(False)
        process.allMuTauPairsLooseMuonIsolation.doSVreco = cms.bool(False)
        process.allMuTauPairsLooseMuonIsolation.doPFMEtSign = cms.bool(False)

    if not isMC:
        removeMCMatching(process, outputInProcess = False)
    else:
        # match pat::Taus to all genJets
        # (including to genJets build from electrons/muons produced in tau --> e/mu decays)
        process.tauGenJetMatch.matched = cms.InputTag("tauGenJets")

    #--------------------------------------------------------------------------------
    # configure PAT trigger matching    
    switchOnTrigger(process, hltProcess = hltProcess, outputModule = '')
    #process.patTrigger.addL1Algos = cms.bool(True)
    # CV: disable L1Algos for now, to prevent error messages
    #
    #     %MSG-e L1GlobalTriggerObjectMapRecord:  PATTriggerProducer:patTrigger
    #
    #       ERROR: The requested algorithm name = L1_DoubleEG1
    #       does not exists in the trigger menu.
    #       Returning zero pointer for getObjectMap
    #
    #     to be printed for every event (06/05/2011)
    process.patTrigger.addL1Algos = cms.bool(False)

    process.patTauTriggerMatchHLTprotoType = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src                   = cms.InputTag("cleanLayer1Taus"),
        matched               = cms.InputTag("patTrigger"),
        matchedCuts           = cms.string('path("HLT_Jet30_v*")'),
        maxDPtRel             = cms.double(1.e+3),                                                    
        maxDeltaR             = cms.double(0.5),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True)
    )
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # compute Pt sum of charged + neutral hadrons and photons within isolation cones of size dR = 0.4/0.6

    process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")

    process.patMuonsLoosePFIsoEmbedded03 = cms.EDProducer("PATMuonPFIsolationEmbedder",
        src = cms.InputTag('patMuons'),                                       
        userFloatName = cms.string('pfLooseIsoPt03'),
        pfCandidateSource = cms.InputTag('pfNoPileUp'),
        chargedHadronIso = cms.PSet(
            ptMin = cms.double(0.5),        
            dRvetoCone = cms.double(-1.),
            dRisoCone = cms.double(0.3)
        ),
        neutralHadronIso = cms.PSet(
            ptMin = cms.double(1.0),        
            dRvetoCone = cms.double(0.08),        
            dRisoCone = cms.double(0.3)
        ),
        photonIso = cms.PSet(
            ptMin = cms.double(1.0),        
            dPhiVeto = cms.double(-1.),
            dEtaVeto = cms.double(-1.),
            dRvetoCone = cms.double(0.05),
            dRisoCone = cms.double(0.3)
        )
    )
    process.patMuonsLoosePFIsoEmbedded04 = process.patMuonsLoosePFIsoEmbedded03.clone(
        src = cms.InputTag('patMuonsLoosePFIsoEmbedded03'),
        userFloatName = cms.string('pfLooseIsoPt04'),
        pfCandidateSource = cms.InputTag('pfNoPileUp'),
        chargedHadronIso = process.patMuonsLoosePFIsoEmbedded03.chargedHadronIso.clone(
            dRisoCone = cms.double(0.4)
        ),
        neutralHadronIso = process.patMuonsLoosePFIsoEmbedded03.neutralHadronIso.clone(
            dRisoCone = cms.double(0.4)
        ),
        photonIso = process.patMuonsLoosePFIsoEmbedded03.photonIso.clone(
            dRisoCone = cms.double(0.4)
        )
    )
    process.patMuonsLoosePFIsoEmbedded06 = process.patMuonsLoosePFIsoEmbedded03.clone(
        src = cms.InputTag('patMuonsLoosePFIsoEmbedded04'),
        userFloatName = cms.string('pfLooseIsoPt06'),
        pfCandidateSource = cms.InputTag('pfNoPileUp'),
        chargedHadronIso = process.patMuonsLoosePFIsoEmbedded03.chargedHadronIso.clone(
            dRisoCone = cms.double(0.6)
        ),
        neutralHadronIso = process.patMuonsLoosePFIsoEmbedded03.neutralHadronIso.clone(
            dRisoCone = cms.double(0.6)
        ),
        photonIso = process.patMuonsLoosePFIsoEmbedded03.photonIso.clone(
            dRisoCone = cms.double(0.6)
        )
    )

    process.patMuonsLoosePFIsoEmbedded = cms.Sequence(
        process.pfNoPileUpSequence
       * process.patMuonsLoosePFIsoEmbedded03 * process.patMuonsLoosePFIsoEmbedded04 * process.patMuonsLoosePFIsoEmbedded06
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
        jetCollectionName = "patJetsAK5Calo",
        patTauProducerPrototype = process.patCaloTauProducer,
        patTauCleanerPrototype = patCaloTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
        addGenInfo = isMC,
        applyTauVertexMatch = applyTauVertexMatch
    )
    process.caloTauSequence = retVal_caloTau["sequence"]

    process.patMuonCaloTauPairs = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuonsLoosePFIsoEmbedded06'),
        srcLeg2 = cms.InputTag(retVal_caloTau["collection"]),
        srcMET = cms.InputTag('patMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False),
        doPFMEtSign = cms.bool(False)
    )
    if hasattr(process.patMuonCaloTauPairs, "nSVfit"):
        delattr(process.patMuonCaloTauPairs, "nSVfit")
    if hasattr(process.patMuonCaloTauPairs, "pfMEtSign"):
        delattr(process.patMuonCaloTauPairs, "pfMEtSign")
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
        jetCollectionName = "patJetsAK5PF",
        patTauProducerPrototype = process.patPFTauProducerFixedCone,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
        addGenInfo = isMC,
        applyTauVertexMatch = applyTauVertexMatch
    )
    process.pfTauSequenceFixedCone = retVal_pfTauFixedCone["sequence"]

    process.patMuonPFTauPairsFixedCone = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuonsLoosePFIsoEmbedded06'),
        srcLeg2 = cms.InputTag(retVal_pfTauFixedCone["collection"]),
        srcMET = cms.InputTag('patPFMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    if hasattr(process.patMuonPFTauPairsFixedCone, "nSVfit"):
        delattr(process.patMuonPFTauPairsFixedCone, "nSVfit")
    if hasattr(process.patMuonPFTauPairsFixedCone, "pfMEtSign"):
        delattr(process.patMuonPFTauPairsFixedCone, "pfMEtSign")
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collection of pat::Tau objects representing PFTaus
    # reconstructed by shrinking signal cone algorithm
    # (plus combinations of muon + tau-jet pairs) 
    #
    switchToPFTauShrinkingCone(process)
    process.patPFTauProducerShrinkingCone = copy.deepcopy(process.patTaus)

    # load TaNC inputs into pat::Tau
    process.load("RecoTauTag.TauTagTools.PFTauMVAInputDiscriminatorTranslator_cfi")
    loadMVAInputsIntoPatTauDiscriminants(process.patPFTauProducerShrinkingCone)

    retVal_pfTauShrinkingCone = patSequenceBuilder(
        process,
        collectionName = [ "patPFTaus", "ShrinkingCone" ],
        jetCollectionName = "patJetsAK5PF",
        patTauProducerPrototype = process.patPFTauProducerShrinkingCone,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
        addGenInfo = isMC,
        applyTauVertexMatch = applyTauVertexMatch
    )
    process.pfTauSequenceShrinkingCone = retVal_pfTauShrinkingCone["sequence"]

    process.patMuonPFTauPairsShrinkingCone = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuonsLoosePFIsoEmbedded06'),
        srcLeg2 = cms.InputTag(retVal_pfTauShrinkingCone["collection"]),
        srcMET = cms.InputTag('patPFMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    if hasattr(process.patMuonPFTauPairsShrinkingCone, "nSVfit"):
        delattr(process.patMuonPFTauPairsShrinkingCone, "nSVfit")
    if hasattr(process.patMuonPFTauPairsShrinkingCone, "pfMEtSign"):
        delattr(process.patMuonPFTauPairsShrinkingCone, "pfMEtSign")
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collection of pat::Tau objects representing PFTaus
    # reconstructed by hadron + strips (HPS) algorithm
    # (plus combinations of muon + tau-jet pairs) 
    #
    # NOTE: switchToPFTauHPS function overwrites process.cleanPatTaus.preselection using HPS specific discriminators;
    #       undo overwriting, in order to prevent run-time errors in case of subsequence _switchToPFTau call,
    #       arising from the fact that HPS specific discriminators are not available for all tau types
    #
    switchToPFTauHPS(process)
    process.cleanPatTaus.preselection = cms.string('')
    process.patPFTauProducerHPS = copy.deepcopy(process.patTaus)

    retVal_pfTauHPS = patSequenceBuilder(
        process,
        collectionName = [ "patPFTaus", "HPS" ],
        jetCollectionName = "patJetsAK5PF",
        patTauProducerPrototype = process.patPFTauProducerHPS,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
        addGenInfo = isMC,
        applyTauVertexMatch = applyTauVertexMatch
    )
    process.pfTauSequenceHPS = retVal_pfTauHPS["sequence"]

    process.patMuonPFTauPairsHPS = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuonsLoosePFIsoEmbedded06'),
        srcLeg2 = cms.InputTag(retVal_pfTauHPS["collection"]),
        srcMET = cms.InputTag('patPFMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    if hasattr(process.patMuonPFTauPairsHPS, "nSVfit"):
        delattr(process.patMuonPFTauPairsHPS, "nSVfit")
    if hasattr(process.patMuonPFTauPairsHPS, "pfMEtSign"):
        delattr(process.patMuonPFTauPairsHPS, "pfMEtSign")
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collection of pat::Tau objects representing PFTaus
    # reconstructed by HPS + TaNC combined tau id. algorithm
    # (plus combinations of muon + tau-jet pairs) 
    #
    switchToPFTauHPSpTaNC(process)
    process.cleanPatTaus.preselection = cms.string('')
    process.patPFTauProducerHPSpTaNC = copy.deepcopy(process.patTaus)

    retVal_pfTauHPSpTaNC = patSequenceBuilder(
        process,
        collectionName = [ "patPFTaus", "HPSpTaNC" ],
        jetCollectionName = "patJetsAK5PF",
        patTauProducerPrototype = process.patPFTauProducerHPSpTaNC,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
        addGenInfo = isMC,
        applyTauVertexMatch = applyTauVertexMatch
    )
    process.pfTauSequenceHPSpTaNC = retVal_pfTauHPSpTaNC["sequence"]

    process.patMuonPFTauPairsHPSpTaNC = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('patMuonsLoosePFIsoEmbedded06'),
        srcLeg2 = cms.InputTag(retVal_pfTauHPSpTaNC["collection"]),
        srcMET = cms.InputTag('patPFMETs'),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    if hasattr(process.patMuonPFTauPairsHPSpTaNC, "nSVfit"):
        delattr(process.patMuonPFTauPairsHPSpTaNC, "nSVfit")
    if hasattr(process.patMuonPFTauPairsHPSpTaNC, "pfMEtSign"):
        delattr(process.patMuonPFTauPairsHPSpTaNC, "pfMEtSign")
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collections of pat::Jets for CaloJets and PFJets
    #
    # NOTE: needed for evaluating jetId for Calo/TCTaus and PFTaus
    #
    jec = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
    if not isMC:
        jec.extend([ 'L2L3Residual' ])
    addJetCollection(process, cms.InputTag('ak5PFJets'),
                     'AK5', 'PF',
                     doJTA            = False,
                     doBTagging       = False,
                     jetCorrLabel     = ('AK5PF', cms.vstring(jec)),
                     doType1MET       = False,
                     genJetCollection = cms.InputTag("ak5GenJets"),
                     doJetID          = True,
                     jetIdLabel       = "ak5",
                     outputModule     = ''
    )

    addJetCollection(process, cms.InputTag('ak5CaloJets'),
                     'AK5', 'Calo',
                     doJTA            = False,
                     doBTagging       = False,
                     jetCorrLabel     = ('AK5Calo', cms.vstring(jec)),
                     doType1MET       = False,
                     genJetCollection = cms.InputTag("ak5GenJets"),
                     doJetID          = True,
                     jetIdLabel       = "ak5",
                     outputModule     = ''
    )
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # replace CaloJets by PFJets for "standard" pat::Jet collection
    switchJetCollection(process, jetCollection = cms.InputTag("ak5PFJets"), outputModule = '')
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
       ##+ process.patTrigger + process.patTriggerEvent
       + process.patMuonsLoosePFIsoEmbedded       
       + process.caloTauSequence
       # store TaNC inputs as discriminators
       + process.produceTancMVAInputDiscriminators
       + process.pfTauSequenceFixedCone + process.pfTauSequenceShrinkingCone + process.pfTauSequenceHPS
       + process.pfTauSequenceHPSpTaNC
       + process.patMuonCaloTauPairs
       + process.patMuonPFTauPairsFixedCone + process.patMuonPFTauPairsShrinkingCone + process.patMuonPFTauPairsHPS
       + process.patMuonPFTauPairsHPSpTaNC
    )

    # return names of "final" collections of CaloTaus/different types of PFTaus
    # to be used as InputTag for further processing
    retVal = {}
    retVal["caloTauCollection"] = retVal_caloTau["collection"]
    retVal["muonCaloTauCollection"] = process.patMuonCaloTauPairs.label()
    retVal["pfTauCollectionFixedCone"] = retVal_pfTauFixedCone["collection"]
    retVal["muonPFTauCollectionFixedCone"] = process.patMuonPFTauPairsFixedCone.label()
    retVal["pfTauCollectionShrinkingCone"] = retVal_pfTauShrinkingCone["collection"]
    retVal["muonPFTauCollectionShrinkingCone"] = process.patMuonPFTauPairsShrinkingCone.label()
    retVal["pfTauCollectionHPS"] = retVal_pfTauHPS["collection"]
    retVal["muonPFTauCollectionHPS"] = process.patMuonPFTauPairsHPS.label()
    retVal["pfTauCollectionHPSpTaNC"] = retVal_pfTauHPSpTaNC["collection"]
    retVal["muonPFTauCollectionHPSpTaNC"] = process.patMuonPFTauPairsHPSpTaNC.label()
    return retVal
