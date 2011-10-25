import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *
from TauAnalysis.CandidateTools.tools.objSelConfigurator import objSelConfigurator
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProduction import configurePatTupleProduction
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildGenericTauSequence
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilderTauIdEffMeasSpecific import buildSequenceTauIdEffMeasSpecific
import PhysicsTools.PatAlgos.tools.helpers as patutils

def configurePatTupleProductionTauIdEffMeasSpecific(process, patSequenceBuilder = buildGenericTauSequence,
                                                    hltProcess = "HLT",
                                                    isMC = False, isEmbedded = False,
                                                    applyZrecoilCorrection = False, runSVfit = True):

    # check that patSequenceBuilder and patTauCleanerPrototype are defined and non-null
    if patSequenceBuilder is None:
        raise ValueError("Undefined 'patSequenceBuilder' Parameter !!")

    #--------------------------------------------------------------------------------
    # produce "basic" PAT objects
    #--------------------------------------------------------------------------------

    process.load("PhysicsTools.PatAlgos.cleaningLayer1.tauCleaner_cfi")
    patCaloTauCleanerPrototype = process.cleanPatTaus.clone(
        preselection = cms.string(''),
        checkOverlaps = cms.PSet(),
        finalCut = cms.string(
            'caloTauTagInfoRef().jetRef().pt() > 5 & abs(caloTauTagInfoRef().jetRef().eta()) < 2.5'
        )
    )
    patPFTauCleanerPrototype = process.cleanPatTaus.clone(
        preselection = cms.string(''),
        checkOverlaps = cms.PSet(),
        finalCut = cms.string(
            'pfJetRef().pt() > 5 & abs(pfJetRef().eta()) < 2.5'
        )
    )

    patTupleConfig = configurePatTupleProduction(process, patSequenceBuilder,
                                                 patPFTauCleanerPrototype, patCaloTauCleanerPrototype, True, hltProcess, isMC, False)

    process.selectedPatMuonsForTauIdEffPFRelIso = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag('selectedPatMuonsVBTFid'),                                                      
        cut = cms.string(
            '(userIsolation("pat::User1Iso")' + \
            ' + max(0., userIsolation("pat::PfNeutralHadronIso") + userIsolation("pat::PfGammaIso")' + \
            '          - 0.5*userIsolation("pat::User2Iso"))) < 0.50*pt'
        ),
        filter = cms.bool(False)
    )
   
    process.producePatTupleTauIdEffMeasSpecific = cms.Sequence(
        process.patDefaultSequence
       + process.selectedPatMuonsForTauIdEffPFRelIso
       + process.caloTauSequence
       # store TaNC inputs as discriminators
       + process.produceTancMVAInputDiscriminators
       + process.pfTauSequenceFixedCone + process.pfTauSequenceShrinkingCone + process.pfTauSequenceHPS
       + process.pfTauSequenceHPSpTaNC
    )

    #--------------------------------------------------------------------------------
    # select "the" primary event Vertex
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.recoVertexSelection_cff")
    process.producePatTupleTauIdEffMeasSpecific += process.selectPrimaryVertex

    #--------------------------------------------------------------------------------
    # produce collections of "global" and "stand-alone" Muons for di-Muon veto
    #--------------------------------------------------------------------------------

    process.selectedPatMuonsLoose = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag('patMuons'),                                  
        cut = cms.string('(isGlobalMuon() | isStandAloneMuon() | isTrackerMuon()) & pt > 5.0'),
        filter = cms.bool(False)                                
    )
    process.producePatTupleTauIdEffMeasSpecific += process.selectedPatMuonsLoose
    
    #--------------------------------------------------------------------------------
    # define Muon momentum scale corrections
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patMuonMomentumCorrection_cfi")
    process.patMuonsMuScleFitCorrectedMomentum.MuonLabel = cms.InputTag('selectedPatMuonsForTauIdEffPFRelIso')
    process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentum

    process.load("TauAnalysis.RecoTools.patLeptonSystematics_cff")
    process.patMuonsMuScleFitCorrectedMomentumShiftUp.MuonLabel = process.patMuonsMuScleFitCorrectedMomentum.MuonLabel
    process.patMuonsMuScleFitCorrectedMomentumShiftDown.MuonLabel = process.patMuonsMuScleFitCorrectedMomentum.MuonLabel
    process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentumShiftUp
    process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentumShiftDown

    #--------------------------------------------------------------------------------
    # define Muon selection for di-Muon veto
    #--------------------------------------------------------------------------------

    process.selectedPatMuonsForTauIdEffZmumuHypotheses = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("patMuons"),
        cut = cms.string('isGlobalMuon() | isTrackerMuon() | isStandAloneMuon()'),
        filter = cms.bool(False)
    )
    
    process.allDiMuPairForTauIdEffZmumuHypotheses = cms.EDProducer("PATDiMuPairProducer",
        useLeadingTausOnly = cms.bool(False),
        srcLeg1 = cms.InputTag('selectedPatMuonsForTauIdEffPFRelIso'),
        srcLeg2 = cms.InputTag('selectedPatMuonsForTauIdEffZmumuHypotheses'),
        dRmin12 = cms.double(0.01),
        srcMET = cms.InputTag(''),
        recoMode = cms.string(""),
        verbosity = cms.untracked.int32(0)
    )

    process.selectedDiMuPairForTauIdEffZmumuHypotheses = cms.EDFilter("PATDiMuPairSelector",
        src = cms.InputTag("allDiMuPairForTauIdEffZmumuHypotheses"),                                   
        cut = cms.string('charge = 0'),
        filter = cms.bool(False)
    )

    process.produceDiMuPairsForTauIdEff = cms.Sequence(
        process.selectedPatMuonsForTauIdEffZmumuHypotheses
       * process.allDiMuPairForTauIdEffZmumuHypotheses * process.selectedDiMuPairForTauIdEffZmumuHypotheses
    )

    process.producePatTupleTauIdEffMeasSpecific += process.produceDiMuPairsForTauIdEff

    #--------------------------------------------------------------------------------
    # define Electron selection (needed for pat::Jet overlap removal)
    #--------------------------------------------------------------------------------
    
    process.load("TauAnalysis/RecoTools/patLeptonSelection_cff")
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatElectrons

    #--------------------------------------------------------------------------------
    # define Jet energy scale corrections
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patJetSystematics_cff")
    # CV: shift jet energy by 3 standard-deviations,
    #     so that template morphing remains an interpolation and no extrapolation is needed
    setattr(process.patJetsJECshiftUp,   "varyByNsigmas", cms.double(3.0))
    setattr(process.patJetsJECshiftDown, "varyByNsigmas", cms.double(3.0))
    process.producePatTupleTauIdEffMeasSpecific += process.prodSmearedJets 

    #--------------------------------------------------------------------------------
    # produce collections of Tau-Jets and of Muon + Tau-jet candidate pairs
    #--------------------------------------------------------------------------------

    # CV: do not run SVFit and PFMetSignificance algorithms,
    #     as they use a significant amount of time (2011/05/27)
    retVal_pfTauFixedCone = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "PFTau", "FixedCone" ], patTupleConfig["pfTauCollectionFixedCone"], False,
                                          None,
                                          'patPFMETs',
                                          isMC = isMC, isEmbedded = isEmbedded,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = False)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauFixedCone["sequence"]
    retVal_pfTauShrinkingCone = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "PFTau", "ShrinkingCone" ], patTupleConfig["pfTauCollectionShrinkingCone"], False,
                                          None,
                                          'patPFMETs',
                                          isMC = isMC, isEmbedded = isEmbedded,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = False)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauShrinkingCone["sequence"]
    # CV: save HPS taus passing the following tau id. discriminators
    #     for measurement of tau charge misidentification rate
    savePFTauHPS = \
        "pt > 15.0 & abs(eta) < 2.5 & " \
       + "tauID('decayModeFinding') > 0.5 & " \
       + "(tauID('byLooseIsolation') > 0.5 |" \
       + " tauID('byLooseIsolationDeltaBetaCorr') > 0.5 |" \
       + " tauID('byLooseCombinedIsolationDeltaBetaCorr') > 0.5) & " \
       + "tauID('againstElectronLoose') > 0.5 & " \
       + "tauID('againstMuonTight') > 0.5"                 
    retVal_pfTauHPS = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "PFTau", "HPS" ], patTupleConfig["pfTauCollectionHPS"], True,
                                          savePFTauHPS,
                                          'patPFMETs',
                                          isMC = isMC, isEmbedded = isEmbedded,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = runSVfit)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauHPS["sequence"]
    # CV: save HPS+TaNC taus passing the following tau id. discriminators
    #     for measurement of tau charge misidentification rate
    savePFTauHPSpTaNC = \
        "pt > 15.0 & abs(eta) < 2.5 & " \
       + "tauID('leadingTrackFinding') > 0.5 & " \
       + "tauID('byTaNCloose') > 0.5 & " \
       + "tauID('againstElectronLoose') > 0.5 & " \
       + "tauID('againstMuonTight') > 0.5"
    retVal_pfTauHPSpTaNC = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "PFTau", "HPSpTaNC" ], patTupleConfig["pfTauCollectionHPSpTaNC"], True,
                                          savePFTauHPSpTaNC,
                                          'patPFMETs',
                                          isMC = isMC, isEmbedded = isEmbedded,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = runSVfit)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauHPSpTaNC["sequence"]

    retVal = {}
    retVal["pfTauFixedCone"] = retVal_pfTauFixedCone
    retVal["pfTauShrinkingCone"] = retVal_pfTauShrinkingCone
    retVal["pfTauHPS"] = retVal_pfTauHPS
    retVal["pfTauHPSpTaNC"] = retVal_pfTauHPSpTaNC
    retVal["algorithms"] = [
        "pfTauFixedCone",
        "pfTauShrinkingCone",
        "pfTauHPS",
        "pfTauHPSpTaNC"
    ]    
    return retVal



