import FWCore.ParameterSet.Config as cms
import copy

from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
import PhysicsTools.PatAlgos.tools.helpers as patutils

from TauAnalysis.Configuration.tools.metTools import *

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
        removeMCMatching(process, ["All"], outputModules = [])
    else:
        # match pat::Taus to all genJets
        # (including to genJets build from electrons/muons produced in tau --> e/mu decays)
        process.tauGenJetMatch.matched = cms.InputTag("tauGenJets")

    #--------------------------------------------------------------------------------
    # configure PAT trigger matching    
    switchOnTrigger(process, hltProcess = hltProcess, outputModule = '')
    # CV: disable L1Algos in MC for now, to prevent error messages
    #
    #     %MSG-e L1GlobalTriggerObjectMapRecord:  PATTriggerProducer:patTrigger
    #
    #       ERROR: The requested algorithm name = L1_DoubleEG1
    #       does not exists in the trigger menu.
    #       Returning zero pointer for getObjectMap
    #
    #     to be printed for every event (06/05/2011)
    #
    #     for Data the L1Algos flag needs to be enabled,
    #     in order for the prescale computation/correction for Data
    #     implemented in TauAnalysis/TauIdEfficiency/bin/FWLiteTauIdEffAnalyzer.cc to work
    if isMC:
        process.patTrigger.addL1Algos = cms.bool(False)
    else:
        process.patTrigger.addL1Algos = cms.bool(True)
        
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

    process.load("RecoMuon/MuonIsolation/muonPFIsolation_cff")    
    patutils.massSearchReplaceAnyInputTag(process.muonPFIsolationDepositsSequence, cms.InputTag('muons1stStep'), cms.InputTag('muons'))
    process.patMuons.isoDeposits = cms.PSet(
        # CV: strings for IsoDeposits defined in PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc
        pfChargedHadrons = cms.InputTag("muPFIsoDepositCharged"),
        pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutral"),
        pfPhotons = cms.InputTag("muPFIsoDepositGamma"),
        user = cms.VInputTag(
            cms.InputTag("muPFIsoDepositChargedAll"),
            cms.InputTag("muPFIsoDepositPU")
        )
    )

    process.patMuons.userIsolation = cms.PSet(
        # CV: strings for Isolation values defined in PhysicsTools/PatAlgos/src/MultiIsolator.cc
        pfChargedHadron = cms.PSet(
            deltaR = cms.double(0.4),
            src = process.patMuons.isoDeposits.pfChargedHadrons,
            vetos = process.muPFIsoValueCharged04.deposits[0].vetos,
            skipDefaultVeto = process.muPFIsoValueCharged04.deposits[0].skipDefaultVeto
        ),
        pfNeutralHadron = cms.PSet(
            deltaR = cms.double(0.4),
            src = process.patMuons.isoDeposits.pfNeutralHadrons,
            vetos = process.muPFIsoValueNeutral04.deposits[0].vetos,
            skipDefaultVeto = process.muPFIsoValueNeutral04.deposits[0].skipDefaultVeto
        ),
        pfGamma = cms.PSet(
            deltaR = cms.double(0.4),
            src = process.patMuons.isoDeposits.pfPhotons,
            vetos = process.muPFIsoValueGamma04.deposits[0].vetos,
            skipDefaultVeto = process.muPFIsoValueGamma04.deposits[0].skipDefaultVeto
        ),
        user = cms.VPSet(
            cms.PSet(
                deltaR = cms.double(0.4),
                src = process.patMuons.isoDeposits.user[0],
                vetos = process.muPFIsoValueChargedAll04.deposits[0].vetos,
                skipDefaultVeto = process.muPFIsoValueChargedAll04.deposits[0].skipDefaultVeto
            ),
            cms.PSet(
                deltaR = cms.double(0.4),
                src = process.patMuons.isoDeposits.user[1],
                vetos = process.muPFIsoValuePU04.deposits[0].vetos,
                skipDefaultVeto = process.muPFIsoValuePU04.deposits[0].skipDefaultVeto
            )
        ) 
    )

    process.patMuonsWithinAcc = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag('patMuons'),
        cut = cms.string("pt > 15. & abs(eta) < 2.1"),
        filter = cms.bool(False)
    )
    process.selectedPatMuonsVBTFid = cms.EDFilter("PATMuonIdSelector",
        src = cms.InputTag('patMuonsWithinAcc'),
        vertexSource = cms.InputTag('selectedPrimaryVertexPosition'),
        beamSpotSource = cms.InputTag('offlineBeamSpot'),                                         
        filter = cms.bool(False)                                         
    )
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
                     outputModules    = []
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
                     outputModules    = []
    )
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # configure Jet Energy Corrections
    #
    process.load("TauAnalysis.Configuration.jetCorrectionParameters_cfi")
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # add pfMET
    process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
    if isMC:
        import PhysicsTools.PatAlgos.tools.helpers as configtools
        configtools.cloneProcessingSnippet(process, process.producePatPFMETCorrections, "NoSmearing")
        process.selectedPatJetsForMETtype1p2CorrNoSmearing.src = cms.InputTag('patJetsNotOverlappingWithLeptonsForMEtUncertainty')
        process.selectedPatJetsForMETtype2CorrNoSmearing.src = process.selectedPatJetsForMETtype1p2CorrNoSmearing.src 

    process.patMEtProductionSequence = cms.Sequence()
    process.patMEtProductionSequence += process.patDefaultSequence

    from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
    doSmearJets = None
    if isMC:
        doSmearJets = True
    else:
        doSmearJets = False
    runMEtUncertainties(
        process,
        electronCollection = '',
        photonCollection = '',
        muonCollection = cms.InputTag('selectedPatMuonsVBTFid'),
        tauCollection = '',
        jetCollection = cms.InputTag('patJetsAK5PF'),
        doSmearJets = doSmearJets,
        doApplyType0corr = True,
        # CV: shift Jet energy by 3 standard-deviations,
        #     so that template morphing remains an interpolation and no extrapolation is needed
        varyByNsigmas = 3.0, 
        addToPatDefaultSequence = False
    )

    if isMC:
        process.patPFMet.addGenMET = cms.bool(True)
        process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string("L3Absolute")
    
        process.patMEtProductionSequence += process.metUncertaintySequence
    else:
        process.patPFMet.addGenMET = cms.bool(False)
        process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string("L2L3Residual")
    
        process.patMEtProductionSequence += process.patJetsAK5PFnotOverlappingWithLeptonsForMEtUncertainty
        process.patMEtProductionSequence += process.producePatPFMETCorrections
    #--------------------------------------------------------------------------------

    pfJetCollection = 'patJetsAK5PFnotOverlappingWithLeptonsForMEtUncertainty'
    pfMEtCollection = 'patType1CorrectedPFMet'
    if isMC:
        pfJetCollection = 'smearedPatJetsAK5PF'

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
        applyTauJEC = False,
        applyTauVertexMatch = applyTauVertexMatch
    )
    process.caloTauSequence = retVal_caloTau["sequence"]
    
    process.patMuonCaloTauPairs = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('selectedPatMuonsVBTFid'),
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
    #switchToPFTauFixedCone(process)
    #process.patPFTauProducerFixedCone = copy.deepcopy(process.patTaus)
    #
    #retVal_pfTauFixedCone = patSequenceBuilder(
    #    process,
    #    collectionName = [ "patPFTaus", "FixedCone" ],
    #    jetCollectionName = pfJetCollection,
    #    patTauProducerPrototype = process.patPFTauProducerFixedCone,
    #    patTauCleanerPrototype = patPFTauCleanerPrototype,
    #    triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
    #    addGenInfo = isMC,
    #    applyTauJEC = False,
    #    applyTauVertexMatch = applyTauVertexMatch
    #)
    #process.pfTauSequenceFixedCone = retVal_pfTauFixedCone["sequence"]
    #
    #process.patMuonPFTauPairsFixedCone = process.allMuTauPairs.clone(
    #    srcLeg1 = cms.InputTag('selectedPatMuonsVBTFid'),
    #    srcLeg2 = cms.InputTag(retVal_pfTauFixedCone["collection"]),
    #    srcMET = cms.InputTag(pfMEtCollection),
    #    srcGenParticles = cms.InputTag(''),
    #    doSVreco = cms.bool(False)
    #)
    #if hasattr(process.patMuonPFTauPairsFixedCone, "nSVfit"):
    #    delattr(process.patMuonPFTauPairsFixedCone, "nSVfit")
    #if hasattr(process.patMuonPFTauPairsFixedCone, "pfMEtSign"):
    #    delattr(process.patMuonPFTauPairsFixedCone, "pfMEtSign")
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collection of pat::Tau objects representing PFTaus
    # reconstructed by shrinking signal cone algorithm
    # (plus combinations of muon + tau-jet pairs) 
    #
    #switchToPFTauShrinkingCone(process)
    #process.patPFTauProducerShrinkingCone = copy.deepcopy(process.patTaus)
    #
    #retVal_pfTauShrinkingCone = patSequenceBuilder(
    #    process,
    #    collectionName = [ "patPFTaus", "ShrinkingCone" ],
    #    jetCollectionName = pfJetCollection,
    #    patTauProducerPrototype = process.patPFTauProducerShrinkingCone,
    #    patTauCleanerPrototype = patPFTauCleanerPrototype,
    #    triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
    #    addGenInfo = isMC,
    #    applyTauJEC = False,
    #    applyTauVertexMatch = applyTauVertexMatch
    #)
    #process.pfTauSequenceShrinkingCone = retVal_pfTauShrinkingCone["sequence"]
    #
    #process.patMuonPFTauPairsShrinkingCone = process.allMuTauPairs.clone(
    #    srcLeg1 = cms.InputTag('selectedPatMuonsVBTFid'),
    #    srcLeg2 = cms.InputTag(retVal_pfTauShrinkingCone["collection"]),
    #    srcMET = cms.InputTag(pfMEtCollection),
    #    srcGenParticles = cms.InputTag(''),
    #    doSVreco = cms.bool(False)
    #)
    #if hasattr(process.patMuonPFTauPairsShrinkingCone, "nSVfit"):
    #    delattr(process.patMuonPFTauPairsShrinkingCone, "nSVfit")
    #if hasattr(process.patMuonPFTauPairsShrinkingCone, "pfMEtSign"):
    #    delattr(process.patMuonPFTauPairsShrinkingCone, "pfMEtSign")
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
        jetCollectionName = pfJetCollection,
        patTauProducerPrototype = process.patPFTauProducerHPS,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
        addGenInfo = isMC,
        applyTauJEC = True,
        applyTauVertexMatch = applyTauVertexMatch
    )
    process.pfTauSequenceHPS = retVal_pfTauHPS["sequence"]

    process.patMuonPFTauPairsHPS = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('selectedPatMuonsVBTFid'),
        srcLeg2 = cms.InputTag(retVal_pfTauHPS["collection"]),
        srcMET = cms.InputTag(pfMEtCollection),
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
        jetCollectionName = pfJetCollection,
        patTauProducerPrototype = process.patPFTauProducerHPSpTaNC,
        patTauCleanerPrototype = patPFTauCleanerPrototype,
        triggerMatcherProtoType = process.patTauTriggerMatchHLTprotoType,
        addGenInfo = isMC,
        applyTauJEC = True,
        applyTauVertexMatch = applyTauVertexMatch
    )
    process.pfTauSequenceHPSpTaNC = retVal_pfTauHPSpTaNC["sequence"]

    process.patMuonPFTauPairsHPSpTaNC = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag('selectedPatMuonsVBTFid'),
        srcLeg2 = cms.InputTag(retVal_pfTauHPSpTaNC["collection"]),
        srcMET = cms.InputTag(pfMEtCollection),
        srcGenParticles = cms.InputTag(''),
        doSVreco = cms.bool(False)
    )
    if hasattr(process.patMuonPFTauPairsHPSpTaNC, "nSVfit"):
        delattr(process.patMuonPFTauPairsHPSpTaNC, "nSVfit")
    if hasattr(process.patMuonPFTauPairsHPSpTaNC, "pfMEtSign"):
        delattr(process.patMuonPFTauPairsHPSpTaNC, "pfMEtSign")
    #--------------------------------------------------------------------------------

    process.patTupleProductionSequence = cms.Sequence(
        process.muonPFIsolationSequence
       + process.patDefaultSequence        
       ##+ process.patTrigger + process.patTriggerEvent
       + process.patMuonsWithinAcc + process.selectedPatMuonsVBTFid
       + process.patMEtProductionSequence 
       + process.caloTauSequence
       # store TaNC inputs as discriminators
       #+ process.produceTancMVAInputDiscriminators
       #+ process.pfTauSequenceFixedCone
       #+ process.pfTauSequenceShrinkingCone
       + process.pfTauSequenceHPS
       + process.pfTauSequenceHPSpTaNC     
       + process.patMuonCaloTauPairs
       #+ process.patMuonPFTauPairsFixedCone
       #+ process.patMuonPFTauPairsShrinkingCone
       + process.patMuonPFTauPairsHPS
       + process.patMuonPFTauPairsHPSpTaNC
    )

    # return names of "final" collections of CaloTaus/different types of PFTaus
    # to be used as InputTag for further processing
    retVal = {}
    retVal["caloTauCollection"] = retVal_caloTau["collection"]
    retVal["muonCaloTauCollection"] = process.patMuonCaloTauPairs.label()
    #retVal["pfTauCollectionFixedCone"] = retVal_pfTauFixedCone["collection"]
    #retVal["muonPFTauCollectionFixedCone"] = process.patMuonPFTauPairsFixedCone.label()
    #retVal["pfTauCollectionShrinkingCone"] = retVal_pfTauShrinkingCone["collection"]
    #retVal["muonPFTauCollectionShrinkingCone"] = process.patMuonPFTauPairsShrinkingCone.label()
    retVal["pfTauCollectionHPS"] = retVal_pfTauHPS["collection"]
    retVal["muonPFTauCollectionHPS"] = process.patMuonPFTauPairsHPS.label()
    retVal["pfTauCollectionHPSpTaNC"] = retVal_pfTauHPSpTaNC["collection"]
    retVal["muonPFTauCollectionHPSpTaNC"] = process.patMuonPFTauPairsHPSpTaNC.label()
    return retVal
