import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from PhysicsTools.PFCandProducer.pfNoPileUp_cff import *
from TauAnalysis.RecoTools.recoVertexSelection_cff import *
from TauAnalysis.RecoTools.patLeptonSelection_cff import *
from TauAnalysis.CandidateTools.muTauPairProduction_cff import *
from TauAnalysis.CandidateTools.diTauPairSelectionAllKinds_cff import *

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *

#--------------------------------------------------------------------------------
# define HLT trigger path
#--------------------------------------------------------------------------------

hltMu = cms.EDFilter("EventSelPluginFilter",
    selector = cms.PSet(
        pluginName = cms.string('hltMu'),             
        pluginType = cms.string('TriggerResultEventSelector'),
        src = cms.InputTag('TriggerResults::HLT'),
        triggerPaths = cms.vstring('HLT_Mu9', 'HLT_IsoMu9', 'HLT_Mu11', 'HLT_IsoMu13_v3', 'HLT_IsoMu13_v4', 'HLT_Mu15_v1')
    )
)

#--------------------------------------------------------------------------------
# define Muon selection
#--------------------------------------------------------------------------------

patMuonsLoosePFIsoEmbedded04 = cms.EDProducer("PATMuonPFIsolationEmbedder",
    src = cms.InputTag('patMuons'),                                       
    userFloatName = cms.string('pfLooseIsoPt04'),
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = cms.PSet(
        ptMin = cms.double(1.0),        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.4)
    ),
    neutralHadronIso = cms.PSet(
        ptMin = cms.double(1000.),        
        dRvetoCone = cms.double(0.15),        
        dRisoCone = cms.double(0.4)
    ),
    photonIso = cms.PSet(
        ptMin = cms.double(1.0),        
        dPhiVeto = cms.double(-1.), # asymmetric Eta x Phi veto region 
        dEtaVeto = cms.double(-1.), # to account for photon conversions in electron isolation case        
        dRvetoCone = cms.double(0.05),
        dRisoCone = cms.double(0.4)
    )
)

patMuonsLoosePFIsoEmbedded06 = patMuonsLoosePFIsoEmbedded04.clone(
    src = cms.InputTag('patMuonsLoosePFIsoEmbedded04'),
    userFloatName = cms.string('pfLooseIsoPt06'),
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = patMuonsLoosePFIsoEmbedded04.chargedHadronIso.clone(
        dRisoCone = cms.double(0.6)
    ),
    neutralHadronIso = patMuonsLoosePFIsoEmbedded04.neutralHadronIso.clone(
        dRisoCone = cms.double(0.6)
    ),
    photonIso = patMuonsLoosePFIsoEmbedded04.photonIso.clone(
        dRisoCone = cms.double(0.6)
    )
)

patMuonsLoosePFIsoEmbedded = cms.Sequence(patMuonsLoosePFIsoEmbedded04 * patMuonsLoosePFIsoEmbedded06)

selectedPatMuonsForTauIdEffGlobal = copy.deepcopy(selectedPatMuonsGlobal)
selectedPatMuonsForTauIdEffGlobal.cut = cms.string('isGlobalMuon()')
selectedPatMuonsForTauIdEffEta21 = copy.deepcopy(selectedPatMuonsEta21)
selectedPatMuonsForTauIdEffEta21.cut = cms.string('abs(eta) < 2.1')
selectedPatMuonsForTauIdEffPt15 = copy.deepcopy(selectedPatMuonsPt15)
selectedPatMuonsForTauIdEffPt15.cut = cms.string('pt > 15.')
selectedPatMuonsForTauIdEffVbTfId = copy.deepcopy(selectedPatMuonsVbTfId)
selectedPatMuonsForTauIdEffVbTfId.beamSpotSource = cms.InputTag("offlineBeamSpot")
selectedPatMuonsForTauIdEffPFRelIso = copy.deepcopy(selectedPatMuonsPFRelIso)
selectedPatMuonsForTauIdEffPFRelIso.chargedHadronIso.ptMin = patMuonsLoosePFIsoEmbedded04.chargedHadronIso.ptMin
selectedPatMuonsForTauIdEffPFRelIso.chargedHadronIso.dRvetoCone = patMuonsLoosePFIsoEmbedded04.chargedHadronIso.dRvetoCone
selectedPatMuonsForTauIdEffPFRelIso.chargedHadronIso.dRisoCone = patMuonsLoosePFIsoEmbedded04.chargedHadronIso.dRisoCone
selectedPatMuonsForTauIdEffPFRelIso.neutralHadronIso.ptMin = patMuonsLoosePFIsoEmbedded04.neutralHadronIso.ptMin 
selectedPatMuonsForTauIdEffPFRelIso.neutralHadronIso.dRvetoCone = patMuonsLoosePFIsoEmbedded04.neutralHadronIso.dRvetoCone
selectedPatMuonsForTauIdEffPFRelIso.neutralHadronIso.dRisoCone = patMuonsLoosePFIsoEmbedded04.neutralHadronIso.dRisoCone
selectedPatMuonsForTauIdEffPFRelIso.photonIso.ptMin = patMuonsLoosePFIsoEmbedded04.photonIso.ptMin
selectedPatMuonsForTauIdEffPFRelIso.photonIso.dRvetoCone = patMuonsLoosePFIsoEmbedded04.photonIso.dRvetoCone
selectedPatMuonsForTauIdEffPFRelIso.photonIso.dRisoCone = patMuonsLoosePFIsoEmbedded04.photonIso.dRisoCone
selectedPatMuonsForTauIdEffPFRelIso.sumPtMax = cms.double(0.30)
selectedPatMuonsForTauIdEffPFRelIso.sumPtMethod = cms.string("relative")
selectedPatMuonsForTauIdEffTrk = copy.deepcopy(selectedPatMuonsTrk)
selectedPatMuonsForTauIdEffTrk.cut = cms.string('innerTrack.isNonnull')
selectedPatMuonsForTauIdEffTrkIP = copy.deepcopy(selectedPatMuonsTrkIP)
selectedPatMuonsForTauIdEffTrkIP.vertexSource = cms.InputTag("selectedPrimaryVertexHighestPtTrackSum")
selectedPatMuonsForTauIdEffTrkIP.IpMax = cms.double(0.05)

patMuonSelConfiguratorForTauIdEff = objSelConfigurator(
    [ selectedPatMuonsForTauIdEffGlobal,
      selectedPatMuonsForTauIdEffEta21,
      selectedPatMuonsForTauIdEffPt15,
      selectedPatMuonsForTauIdEffVbTfId,
      selectedPatMuonsForTauIdEffPFRelIso,
      selectedPatMuonsForTauIdEffTrk,
      selectedPatMuonsForTauIdEffTrkIP ],
    src = "patMuonsLoosePFIsoEmbedded06",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectPatMuonsForTauIdEff = patMuonSelConfiguratorForTauIdEff.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------
# define loose Tau-jet candidate selection
#--------------------------------------------------------------------------------

patTausLoosePFIsoEmbedded04 = cms.EDProducer("PATTauPFIsolationEmbedder",
    src = cms.InputTag('patTaus'),
    userFloatName = cms.string('pfLooseIsoPt04'),
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = cms.PSet(
        ptMin = cms.double(1.0),        
        dRvetoCone = cms.double(0.15),
        dRisoCone = cms.double(0.4)
    ),
    neutralHadronIso = cms.PSet(
        ptMin = cms.double(1000.),        
        dRvetoCone = cms.double(0.15),        
        dRisoCone = cms.double(0.4)
    ),
    photonIso = cms.PSet(
        ptMin = cms.double(1.5),        
        dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
        dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
        dRvetoCone = cms.double(0.15),
        dRisoCone = cms.double(0.4)
    )
)

patTausLoosePFIsoEmbedded06 = patTausLoosePFIsoEmbedded04.clone(
    src = cms.InputTag('patTausLoosePFIsoEmbedded04'),
    userFloatName = cms.string('pfLooseIsoPt06'),
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = patTausLoosePFIsoEmbedded04.chargedHadronIso.clone(
        dRisoCone = cms.double(0.6)
    ),
    neutralHadronIso = patTausLoosePFIsoEmbedded04.neutralHadronIso.clone(
        dRisoCone = cms.double(0.6)
    ),
    photonIso = patTausLoosePFIsoEmbedded04.photonIso.clone(
        dRisoCone = cms.double(0.6)
    )
)

patTausLoosePFIsoEmbedded = cms.Sequence(patTausLoosePFIsoEmbedded04 * patTausLoosePFIsoEmbedded06)

selectedPatTausForTauIdEffAntiOverlapWithMuonsVeto = copy.deepcopy(selectedPatTausForMuTauAntiOverlapWithMuonsVeto)
selectedPatTausForTauIdEffAntiOverlapWithMuonsVeto.dRmin = cms.double(0.7)
selectedPatTausForTauIdEffAntiOverlapWithMuonsVeto.srcNotToBeFiltered = cms.VInputTag("selectedPatMuonsForTauIdEffGlobalCumulative")
selectedPatTausForTauIdEffEta23 = copy.deepcopy(selectedPatTausForMuTauEta23)
selectedPatTausForTauIdEffEta23.cut = cms.string("abs(eta) < 2.3")
selectedPatTausForTauIdEffPt20 = copy.deepcopy(selectedPatTausForMuTauPt20)
selectedPatTausForTauIdEffPt20.cut = cut = cms.string("pt > 20.")
selectedPatTausForTauIdEffLeadTrk = copy.deepcopy(selectedPatTausForMuTauLeadTrk)
selectedPatTausForTauIdEffLeadTrk.cut = cms.string('tauID("leadingTrackFinding") > 0.5')
selectedPatTausForTauIdEffLeadTrkPt = copy.deepcopy(selectedPatTausForMuTauLeadTrkPt)
selectedPatTausForTauIdEffLeadTrkPt.cut = cms.string('tauID("leadingTrackPtCut") > 0.5')
selectedPatTausForTauIdEffMuonVeto = copy.deepcopy(selectedPatTausForMuTauMuonVeto)
selectedPatTausForTauIdEffMuonVeto.cut = cms.string('tauID("againstMuon") > 0.5')
selectedPatTausForTauIdEffCaloMuonVeto = copy.deepcopy(selectedPatTausForMuTauCaloMuonVeto)
selectedPatTausForTauIdEffCaloMuonVeto.cut = cms.string('tauID("againstCaloMuon") > 0.5')
#selectedPatTausForTauIdEffCaloMuonVeto.cut = cms.string('tauID("againstCaloMuon") > -1.')
selectedPatTausForTauIdEffElectronVeto = copy.deepcopy(selectedPatTausForMuTauElectronVeto)
#selectedPatTausForTauIdEffElectronVeto.cut = cms.string('tauID("againstElectron") > 0.5')
selectedPatTausForTauIdEffElectronVeto.cut = cms.string('leadPFCand().isNonnull() & leadPFCand().mva_e_pi() < 0.6')
selectedPatTausForTauIdEffEcalCrackVeto = copy.deepcopy(selectedPatTausEcalCrackVeto)
selectedPatTausForTauIdEffEcalCrackVeto.cut = cms.string("abs(eta) < 1.460 | abs(eta) > 1.558")

patTauSelConfiguratorForTauIdEff = objSelConfigurator(
    [ selectedPatTausForTauIdEffAntiOverlapWithMuonsVeto,
      selectedPatTausForTauIdEffEta23,
      selectedPatTausForTauIdEffPt20,
      selectedPatTausForTauIdEffLeadTrk,
      selectedPatTausForTauIdEffLeadTrkPt,
      selectedPatTausForTauIdEffMuonVeto,
      selectedPatTausForTauIdEffCaloMuonVeto,
      selectedPatTausForTauIdEffElectronVeto,
      selectedPatTausForTauIdEffEcalCrackVeto ],
    src = "patTausLoosePFIsoEmbedded06",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectPatTausForTauIdEff = patTauSelConfiguratorForTauIdEff.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------
# define selection of Muon + Tau-jet candidate pairs
#--------------------------------------------------------------------------------

muTauPairsForTauIdEff = allMuTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatTausForTauIdEffElectronVetoCumulative'),
    srcMET = cms.InputTag('patPFMETs'),
    srcGenParticles = cms.InputTag('genParticles')
)    

selectedMuTauPairsForTauIdEffAntiOverlapVeto = selectedMuTauPairsAntiOverlapVeto.clone(
    cut = cms.string('dR12 > 0.7')
)    
selectedMuTauPairsForTauIdEffMt1MET = selectedMuTauPairsMt1MET.clone(
    cut = cms.string('mt1MET < 40.')
)    
selectedMuTauPairsForTauIdEffPzetaDiff = selectedMuTauPairsPzetaDiff.clone(
    cut = cms.string('(pZeta - 1.5*pZetaVis) > -20.')
)       
selectedMuTauPairsForTauIdEffZeroCharge = selectedMuTauPairsZeroCharge.clone(            
    cut = cms.string('(leg1.charge + leg2.leadPFChargedHadrCand.charge) = 0')
)

patMuTauPairSelConfiguratorForTauIdEff = objSelConfigurator(
    [ selectedMuTauPairsForTauIdEffAntiOverlapVeto,
      selectedMuTauPairsForTauIdEffMt1MET,
      selectedMuTauPairsForTauIdEffPzetaDiff,
      selectedMuTauPairsForTauIdEffZeroCharge ],
    src = "muTauPairsForTauIdEff",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectMuTauPairsForTauIdEff = patMuTauPairSelConfiguratorForTauIdEff.configure(pyNameSpace = locals())

produceMuTauPairsForTauIdEff = cms.Sequence(muTauPairsForTauIdEff * selectMuTauPairsForTauIdEff)

selectedMuTauPairFilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuTauPairsForTauIdEffAntiOverlapVetoCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                   
)

#--------------------------------------------------------------------------------
# veto events containing di-Muon pairs
# (the hypothesis being that the pair results from a Z --> mu+ mu- decay)
#--------------------------------------------------------------------------------

selectedPatMuonsForZmumuHypothesesMuonTrack = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string('isGlobalMuon() | isTrackerMuon() | isStandAloneMuon()'),
    filter = cms.bool(False)
)

selectedPatMuonsForZmumuHypothesesPt10 = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsForZmumuHypothesesMuonTrack"),
    cut = cms.string('pt > 10'),
    filter = cms.bool(False)
)

allDiMuPairZmumuHypotheses = cms.EDProducer("PATDiMuPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatMuonsForZmumuHypothesesPt10'),
    dRmin12 = cms.double(0.3),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)
)

selectedDiMuPairZmumuHypotheses = cms.EDFilter("PATDiMuPairSelector",
    src = cms.InputTag("allDiMuPairZmumuHypotheses"),                                   
    cut = cms.string('charge = 0'),
    filter = cms.bool(False)
)

produceDiMuPairs = cms.Sequence(
    selectedPatMuonsForZmumuHypothesesMuonTrack * selectedPatMuonsForZmumuHypothesesPt10
   * allDiMuPairZmumuHypotheses * selectedDiMuPairZmumuHypotheses
)

diMuPairZmumuHypothesisVeto = cms.EDFilter("PATCandViewMaxFilter",
    src = cms.InputTag('selectedDiMuPairZmumuHypotheses'),
    maxNumber = cms.uint32(0)                                                            
)

#--------------------------------------------------------------------------------
# define event selection sequence
#--------------------------------------------------------------------------------

muonPFTauSkimPath = cms.Path(
    hltMu
   + selectPrimaryVertex
   + patMuonsLoosePFIsoEmbedded + selectPatMuonsForTauIdEff
   + patTausLoosePFIsoEmbedded + selectPatTausForTauIdEff
   + produceMuTauPairsForTauIdEff + selectedMuTauPairFilter
   + produceDiMuPairs + diMuPairZmumuHypothesisVeto
   + dataQualityFilters
)

tauIdEffSampleEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'muonPFTauSkimPath'
        )
    )
)
