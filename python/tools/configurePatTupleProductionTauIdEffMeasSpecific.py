import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *
from TauAnalysis.CandidateTools.tools.objSelConfigurator import objSelConfigurator
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProduction import configurePatTupleProduction
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildGenericTauSequence
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilderTauIdEffMeasSpecific import buildSequenceTauIdEffMeasSpecific

def configurePatTupleProductionTauIdEffMeasSpecific(process, patSequenceBuilder = buildGenericTauSequence,
                                                    hltProcess = "HLT",
                                                    addGenInfo = False, applyZrecoilCorrection = False):

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
                                                 patPFTauCleanerPrototype, patCaloTauCleanerPrototype, True, hltProcess, addGenInfo)
   
    process.producePatTupleTauIdEffMeasSpecific = cms.Sequence(
        process.patDefaultSequence
       + process.patMuonsLoosePFIsoEmbedded
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

    process.patMuonsGlobal = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag('patMuons'),                                            
        cut = cms.string('isGlobalMuon()'),
        filter = cms.bool(False)
    )
    process.producePatTupleTauIdEffMeasSpecific += process.patMuonsGlobal

    process.patMuonsStandAlone = process.patMuonsGlobal.clone(
        cut = cms.string('isStandAloneMuon()')
    )
    process.producePatTupleTauIdEffMeasSpecific += process.patMuonsStandAlone   

    #--------------------------------------------------------------------------------
    # define Muon momentum scale corrections
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patMuonMomentumCorrection_cfi")
    process.patMuonsMuScleFitCorrectedMomentum.MuonLabel = cms.InputTag('patMuonsLoosePFIsoEmbedded06')
    process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentum

    process.load("TauAnalysis.RecoTools.patLeptonSystematics_cff")
    process.patMuonsMuScleFitCorrectedMomentumShiftUp.MuonLabel = process.patMuonsMuScleFitCorrectedMomentum.MuonLabel
    process.patMuonsMuScleFitCorrectedMomentumShiftDown.MuonLabel = process.patMuonsMuScleFitCorrectedMomentum.MuonLabel
    process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentumShiftUp
    process.producePatTupleTauIdEffMeasSpecific += process.patMuonsMuScleFitCorrectedMomentumShiftDown
            
    #--------------------------------------------------------------------------------
    # define Muon selection
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patLeptonSelection_cff")
    process.selectedPatMuonsForTauIdEffGlobal   = process.selectedPatMuonsGlobal.clone(
        cut = cms.string('isGlobalMuon()')
    )
    process.selectedPatMuonsForTauIdEffEta21    = process.selectedPatMuonsEta21.clone(
        cut = cms.string('abs(eta) < 2.1')
    )
    process.selectedPatMuonsForTauIdEffPt15     = process.selectedPatMuonsPt15.clone(
        cut = cms.string('pt > 15.')
    )
    process.selectedPatMuonsForTauIdEffVbTfId   = process.selectedPatMuonsVbTfId.clone(
        beamSpotSource = cms.InputTag("offlineBeamSpot")
    )
    process.selectedPatMuonsForTauIdEffPFRelIso = process.selectedPatMuonsPFRelIso.clone(
        chargedHadronIso = process.patMuonsLoosePFIsoEmbedded04.chargedHadronIso,
        neutralHadronIso = process.patMuonsLoosePFIsoEmbedded04.neutralHadronIso,
        photonIso        = process.patMuonsLoosePFIsoEmbedded04.photonIso,
        sumPtMax         = cms.double(0.30),
        sumPtMethod      = cms.string("relative")
    )
    process.selectedPatMuonsForTauIdEffTrk      = process.selectedPatMuonsTrk.clone(
        cut = cms.string('innerTrack.isNonnull')
    )
    process.selectedPatMuonsForTauIdEffTrkIP    = process.selectedPatMuonsTrkIP.clone(
        vertexSource = cms.InputTag("selectedPrimaryVertexHighestPtTrackSum"),
        IpMax = cms.double(0.05)
    )
    
    patMuonSelConfiguratorForTauIdEff = objSelConfigurator(
        [ process.selectedPatMuonsForTauIdEffGlobal,
          process.selectedPatMuonsForTauIdEffEta21,
          process.selectedPatMuonsForTauIdEffPt15,
          process.selectedPatMuonsForTauIdEffVbTfId,
          process.selectedPatMuonsForTauIdEffPFRelIso,
          process.selectedPatMuonsForTauIdEffTrk,
          process.selectedPatMuonsForTauIdEffTrkIP ],
        src = "patMuonsMuScleFitCorrectedMomentum",
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
        cut = cms.string('pt > 10 & (isGlobalMuon() | isTrackerMuon() | isStandAloneMuon())'),
        filter = cms.bool(False)
    )
    
    process.allDiMuPairForTauIdEffZmumuHypotheses = cms.EDProducer("PATDiMuPairProducer",
        useLeadingTausOnly = cms.bool(False),
        srcLeg1 = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
        srcLeg2 = cms.InputTag('selectedPatMuonsForTauIdEffZmumuHypotheses'),
        dRmin12 = cms.double(0.3),
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

    process.producePatTupleTauIdEffMeasSpecific += process.selectPatElectrons

    #--------------------------------------------------------------------------------
    # define Jet energy scale corrections
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patJetSystematics_cff")
    process.producePatTupleTauIdEffMeasSpecific += process.prodSmearedJets 

    #--------------------------------------------------------------------------------
    # produce collections of Tau-Jets and of Muon + Tau-jet candidate pairs
    #--------------------------------------------------------------------------------

    caloTauPreselectionCriteria = [
        [ "LeadTrk",         'tauID("leadingTrackFinding") > 0.5'                       ],
        [ "LeadTrkPt",       'tauID("leadingTrackPtCut") > 0.5'                         ],
        [ "MuonVeto",        'tauID("againstMuon") > 0.5'                               ],
        [ "ElectronVeto",    'tauID("againstElectron") > 0.5'                           ],
        [ "EcalCrackVeto",   'abs(eta) < 1.460 | abs(eta) > 1.558'                      ]
    ]

    pfTauConeIsoPreselectionCriteria = [
        [ "LeadTrk",         'tauID("leadingTrackFinding") > 0.5'                       ],
        [ "LeadTrkPt",       'tauID("leadingTrackPtCut") > 0.5'                         ],
        [ "MuonVeto",        'tauID("againstMuon") > 0.5'                               ],
        [ "ElectronVeto",    'leadPFCand().isNonnull() & leadPFCand().mva_e_pi() < 0.6' ],
        [ "EcalCrackVeto",   'abs(eta) < 1.460 | abs(eta) > 1.558'                      ]
    ]

    pfTauHPSpreselectionCriteria = [
        [ "DecayMode",       'tauID("decayModeFinding") > 0.5'                          ],
        [ "MuonVeto",        'tauID("againstMuonTight") > 0.5'                          ],
        [ "ElectronVeto",    'tauID("againstElectronLoose") > 0.5'                      ],
        [ "EcalCrackVeto",   'abs(eta) < 1.460 | abs(eta) > 1.558'                      ]
    ]

    pfTauHPSpTaNCpreselectionCriteria = [
        [ "LeadTrk",         'tauID("leadingTrackFinding") > 0.5'                       ],
        #[ "LeadTrkPt",       'tauID("leadingTrackPtCut") > 0.5'                         ],
        [ "LeadPionPt",      'tauID("leadingPionPtCut") > 0.5'                          ],
        [ "MuonVeto",        'tauID("againstMuonTight") > 0.5'                          ],
        #[ "CaloMuonVeto",    'tauID("againstCaloMuon") > 0.5'                           ],
        [ "ElectronVeto",    'tauID("againstElectronLoose") > 0.5'                      ],
        [ "EcalCrackVeto",   'abs(eta) < 1.460 | abs(eta) > 1.558'                      ]
    ]
    
    retVal_caloTau = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "CaloTau", "" ], patTupleConfig["caloTauCollection"],
                                          'patMETs',
                                          caloTauPreselectionCriteria,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = False)
    process.producePatTupleTauIdEffMeasSpecific += retVal_caloTau["sequence"]
    retVal_pfTauFixedCone = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "PFTau", "FixedCone" ], patTupleConfig["pfTauCollectionFixedCone"],
                                          'patPFMETs',
                                          pfTauConeIsoPreselectionCriteria,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = True)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauFixedCone["sequence"]
    retVal_pfTauShrinkingCone = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "PFTau", "ShrinkingCone" ], patTupleConfig["pfTauCollectionShrinkingCone"],
                                          'patPFMETs',
                                          pfTauConeIsoPreselectionCriteria,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = True)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauShrinkingCone["sequence"]
    retVal_pfTauHPS = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "PFTau", "HPS" ], patTupleConfig["pfTauCollectionHPS"],
                                          'patPFMETs',
                                          pfTauHPSpreselectionCriteria,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = True)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauHPS["sequence"]
    retVal_pfTauHPSpTaNC = \
        buildSequenceTauIdEffMeasSpecific(process,
                                          'selectedPatMuonsForTauIdEffTrkIP',
                                          [ "PFTau", "HPSpTaNC" ], patTupleConfig["pfTauCollectionHPSpTaNC"],
                                          'patPFMETs',
                                          pfTauHPSpTaNCpreselectionCriteria,
                                          applyZrecoilCorrection = applyZrecoilCorrection, runSVfit = True)
    process.producePatTupleTauIdEffMeasSpecific += retVal_pfTauHPSpTaNC["sequence"]

    retVal = {}
    retVal["caloTau"]= retVal_caloTau
    retVal["pfTauFixedCone"] = retVal_pfTauFixedCone
    retVal["pfTauShrinkingCone"] = retVal_pfTauShrinkingCone
    retVal["pfTauHPS"] = retVal_pfTauHPS
    retVal["pfTauHPSpTaNC"] = retVal_pfTauHPSpTaNC
    retVal["algorithms"] = [
        "caloTau",
        "pfTauFixedCone",
        "pfTauShrinkingCone",
        "pfTauHPS",
        "pfTauHPSpTaNC"
    ]    
    return retVal



