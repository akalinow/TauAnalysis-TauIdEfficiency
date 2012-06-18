import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName
from TauAnalysis.CandidateTools.tools.objSelConfigurator import objSelConfigurator
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProduction import configurePatTupleProduction
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildGenericTauSequence
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilderTauIdEffMeasSpecific import buildSequenceTauIdEffMeasSpecific
import PhysicsTools.PatAlgos.tools.helpers as patutils
from PhysicsTools.PatAlgos.tools.jetTools import *

def configurePatTupleProductionTauIdEffMeasSpecific(process, patSequenceBuilder = buildGenericTauSequence,
                                                    hltProcess = "HLT",
                                                    isMC = False, isEmbedded = False,
                                                    runSVfit = False):

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

    if isMC:
       addJetCollection(process, cms.InputTag('smearedAK5PFJets'),
                        'SmearedAK5', 'PF',
                        doJTA            = False,
                        doBTagging       = True,
                        jetCorrLabel     = None, # CV: jet corrections already applied on reco::PFJet input
                        doType1MET       = False,
                        genJetCollection = cms.InputTag("ak5GenJets"),
                        doJetID          = True,
                        jetIdLabel       = "ak5",
                        outputModules    = []
        )
    else:
        addJetCollection(process, cms.InputTag('calibratedAK5PFJets'),
                        'CalibratedAK5', 'PF',
                        doJTA            = False,
                        doBTagging       = True,
                        jetCorrLabel     = None, # CV: jet corrections already applied on reco::PFJet input
                        doType1MET       = False,
                        doJetID          = True,
                        jetIdLabel       = "ak5",
                        outputModules    = []
        )

    # add "raw" (uncorrected) CaloMET,
    # needed to parametrize efficiency turn-on of HLT_IsoMu15_L1ETM20 cross-trigger
    process.patCaloMet = process.patMETs.clone(
        metSource = cms.InputTag('met'),
        addMuonCorrections = cms.bool(False)
    )

    process.patCaloMetNoHF = process.patMETs.clone(
        metSource = cms.InputTag('metNoHF'),
        addMuonCorrections = cms.bool(False)
    )
        
    process.producePatTupleTauIdEffMeasSpecific = cms.Sequence(
        process.patTupleProductionSequence
       + process.caloTauSequence
       # store TaNC inputs as discriminators
       #+ process.produceTancMVAInputDiscriminators
       #+ process.pfTauSequenceFixedCone
       #+ process.pfTauSequenceShrinkingCone
       + process.pfTauSequenceHPS
       + process.pfTauSequenceHPSpTaNC
       + process.patCaloMet
       + process.patCaloMetNoHF
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
    if isMC:
        process.poolDBESSourceMuScleFitCentralValue.toGet[0].tag = cms.string('MuScleFit_Scale_Z_MC_Startup_innerTrack')
    else:
        process.poolDBESSourceMuScleFitCentralValue.toGet[0].tag = cms.string('MuScleFit_Scale_Z_36_invPb_innerTrack_Dec22_v1')
    process.patMuonsMuScleFitCorrectedMomentum.MuonLabel = cms.InputTag('selectedPatMuonsVBTFid')
    # CV: MuScleFit muon momentum corrections do not work in CMSSW_5_2_x (May 4th 2012)
    ##process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentum

    ##process.load("TauAnalysis.RecoTools.patLeptonSystematics_cff")
    ##process.patMuonsMuScleFitCorrectedMomentumShiftUp.MuonLabel = process.patMuonsMuScleFitCorrectedMomentum.MuonLabel
    ##process.patMuonsMuScleFitCorrectedMomentumShiftDown.MuonLabel = process.patMuonsMuScleFitCorrectedMomentum.MuonLabel
    ##process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentumShiftUp
    ##process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentumShiftDown

    #--------------------------------------------------------------------------------
    # define Muon selection
    #--------------------------------------------------------------------------------

    process.selectedPatMuonsForTauIdEffPFRelIso = cms.EDFilter("PATMuonSelector",
        cut = cms.string(
            '(userIsolation("pat::User1Iso")' + \
            ' + max(0., userIsolation("pat::PfNeutralHadronIso") + userIsolation("pat::PfGammaIso")' + \
            '          - 0.5*userIsolation("pat::User2Iso"))) < 0.50*pt'
        ),
        filter = cms.bool(False)
    )

    patMuonSelConfiguratorForTauIdEff = objSelConfigurator(
        [ process.selectedPatMuonsForTauIdEffPFRelIso ],
        # CV: MuScleFit muon momentum corrections do not work in CMSSW_5_2_x (May 4th 2012)
        #src = "patMuonsMuScleFitCorrectedMomentum",
        src = "selectedPatMuonsVBTFid",
        pyModuleName = __name__,
        doSelIndividual = False
    )
    setattr(patMuonSelConfiguratorForTauIdEff, "systematics", muonSystematics)
    process.selectPatMuonsForTauIdEff = patMuonSelConfiguratorForTauIdEff.configure(process = process)
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatMuonsForTauIdEff
 
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
        srcLeg1 = cms.InputTag(composeModuleName(['selectedPatMuonsForTauIdEffPFRelIso', "cumulative"])),
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
    #retVal_pfTauFixedCone = \
    #    buildSequenceTauIdEffMeasSpecific(process,
    #                                      'selectedPatMuonsForTauIdEffPFRelIso',
    #                                      [ "PFTau", "FixedCone" ], patTupleConfig["pfTauCollectionFixedCone"], False,
    #                                      None,
    #                                      'patType1CorrectedPFMet',
    #                                      isMC = isMC, isEmbedded = isEmbedded,
    #                                      runSVfit = False)
    #process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauFixedCone["sequence"]
    #retVal_pfTauShrinkingCone = \
    #    buildSequenceTauIdEffMeasSpecific(process,
    #                                      'selectedPatMuonsForTauIdEffPFRelIso',
    #                                      [ "PFTau", "ShrinkingCone" ], patTupleConfig["pfTauCollectionShrinkingCone"], False,
    #                                      None,
    #                                      'patType1CorrectedPFMet',
    #                                      isMC = isMC, isEmbedded = isEmbedded,
    #                                      runSVfit = False)
    #process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauShrinkingCone["sequence"]
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
                                          'selectedPatMuonsForTauIdEffPFRelIso',
                                          [ "PFTau", "HPS" ], patTupleConfig["pfTauCollectionHPS"], True,
                                          savePFTauHPS,
                                          'patType1CorrectedPFMet',
                                          isMC = isMC, isEmbedded = isEmbedded,
                                          runSVfit = runSVfit)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauHPS["sequence"]
    # CV: save HPS+TaNC taus passing the following tau id. discriminators
    #     for measurement of tau charge misidentification rate
    #savePFTauHPSpTaNC = \
    #    "pt > 15.0 & abs(eta) < 2.5 & " \
    #   + "tauID('leadingTrackFinding') > 0.5 & " \
    #   + "tauID('byTaNCloose') > 0.5 & " \
    #   + "tauID('againstElectronLoose') > 0.5 & " \
    #   + "tauID('againstMuonTight') > 0.5"
    #retVal_pfTauHPSpTaNC = \
    #    buildSequenceTauIdEffMeasSpecific(process,
    #                                      'selectedPatMuonsForTauIdEffPFRelIso',
    #                                      [ "PFTau", "HPSpTaNC" ], patTupleConfig["pfTauCollectionHPSpTaNC"], True,
    #                                      savePFTauHPSpTaNC,
    #                                      'patType1CorrectedPFMet',
    #                                      isMC = isMC, isEmbedded = isEmbedded,
    #                                      runSVfit = runSVfit)
    #process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauHPSpTaNC["sequence"]

    retVal = {}
    #retVal["pfTauFixedCone"] = retVal_pfTauFixedCone
    #retVal["pfTauShrinkingCone"] = retVal_pfTauShrinkingCone
    retVal["pfTauHPS"] = retVal_pfTauHPS
    #retVal["pfTauHPSpTaNC"] = retVal_pfTauHPSpTaNC
    retVal["algorithms"] = [
        #"pfTauFixedCone",
        #"pfTauShrinkingCone",
        "pfTauHPS",
        #"pfTauHPSpTaNC"
    ]    
    return retVal



