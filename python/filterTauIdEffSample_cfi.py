import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from PhysicsTools.PFCandProducer.pfNoPileUp_cff import *
from TauAnalysis.RecoTools.recoVertexSelection_cff import *
from TauAnalysis.RecoTools.patMuonSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *
from TauAnalysis.CandidateTools.muTauPairProduction_cff import *
from TauAnalysis.CandidateTools.muTauPairSelection_cfi import *

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
    src = cms.InputTag('cleanPatMuons'),                                       
    userFloatName = cms.string('pfLooseIsoPt04'),
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = cms.PSet(
        ptMin = cms.double(0.5),        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.4)
    ),
    neutralHadronIso = cms.PSet(
        ptMin = cms.double(1000.),        
        dRvetoCone = cms.double(0.08),        
        dRisoCone = cms.double(0.4)
    ),
    photonIso = cms.PSet(
        ptMin = cms.double(0.5),        
        dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
        dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
        dRvetoCone = cms.double(-1.),
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

selectedPatMuonsGlobal.cut = cms.string('isGlobalMuon()')
selectedPatMuonsEta21.cut = cms.string('abs(eta) < 2.1')
selectedPatMuonsPt15.cut = cms.string('pt > 15.')
selectedPatMuonsVbTfId.beamSpotSource = cms.InputTag("offlineBeamSpot")
selectedPatMuonsPFRelIso.chargedHadronIso.ptMin = patMuonsLoosePFIsoEmbedded04.chargedHadronIso.ptMin
selectedPatMuonsPFRelIso.chargedHadronIso.dRvetoCone = patMuonsLoosePFIsoEmbedded04.chargedHadronIso.dRvetoCone
selectedPatMuonsPFRelIso.chargedHadronIso.dRisoCone = patMuonsLoosePFIsoEmbedded04.chargedHadronIso.dRisoCone
selectedPatMuonsPFRelIso.neutralHadronIso.ptMin = patMuonsLoosePFIsoEmbedded04.neutralHadronIso.ptMin 
selectedPatMuonsPFRelIso.neutralHadronIso.dRvetoCone = patMuonsLoosePFIsoEmbedded04.neutralHadronIso.dRvetoCone
selectedPatMuonsPFRelIso.neutralHadronIso.dRisoCone = patMuonsLoosePFIsoEmbedded04.neutralHadronIso.dRisoCone
selectedPatMuonsPFRelIso.photonIso.ptMin = patMuonsLoosePFIsoEmbedded04.photonIso.ptMin
selectedPatMuonsPFRelIso.photonIso.dRvetoCone = patMuonsLoosePFIsoEmbedded04.photonIso.dRvetoCone
selectedPatMuonsPFRelIso.photonIso.dRisoCone = patMuonsLoosePFIsoEmbedded04.photonIso.dRisoCone
selectedPatMuonsPFRelIso.sumPtMax = cms.double(0.30)
selectedPatMuonsPFRelIso.sumPtMethod = cms.string("relative")
selectedPatMuonsTrk.cut = cms.string('innerTrack.isNonnull')
selectedPatMuonsTrkIP.vertexSource = cms.InputTag("selectedPrimaryVertexHighestPtTrackSum")
selectedPatMuonsTrkIP.IpMax = cms.double(0.05)

patMuonSelConfigurator = objSelConfigurator(
    [ selectedPatMuonsGlobal,
      selectedPatMuonsEta21,
      selectedPatMuonsPt15,
      selectedPatMuonsVbTfId,
      selectedPatMuonsPFRelIso,
      selectedPatMuonsTrk,
      selectedPatMuonsTrkIP ],
    src = "patMuonsLoosePFIsoEmbedded06",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectPatMuons = patMuonSelConfigurator.configure(pyNameSpace = locals())

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

selectedPatTausForMuTauAntiOverlapWithMuonsVeto.dRmin = cms.double(0.7)
selectedPatTausForMuTauAntiOverlapWithMuonsVeto.srcNotToBeFiltered = cms.VInputTag("selectedPatMuonsGlobalCumulative")
selectedPatTausForMuTauEta23.cut = cms.string("abs(eta) < 2.3")
selectedPatTausForMuTauPt20.cut = cut = cms.string("pt > 20.")
selectedPatTausForMuTauLeadTrk.cut = cms.string('tauID("leadingTrackFinding") > 0.5')
selectedPatTausForMuTauLeadTrkPt.cut = cms.string('tauID("leadingTrackPtCut") > 0.5')
selectedPatTausForMuTauMuonVeto.cut = cms.string('tauID("againstMuon") > 0.5')
selectedPatTausForMuTauCaloMuonVeto.cut = cms.string('tauID("againstCaloMuon") > 0.5')
selectedPatTausForMuTauElectronVeto.cut = cms.string('tauID("againstElectron") > 0.5')
selectedPatTausForMuTauEcalCrackVeto = copy.deepcopy(selectedPatTausEcalCrackVeto)
selectedPatTausForMuTauEcalCrackVeto.cut = cms.string("abs(eta) < 1.460 | abs(eta) > 1.558")

patTauSelConfiguratorForMuTau = objSelConfigurator(
    [ selectedPatTausForMuTauAntiOverlapWithMuonsVeto,
      selectedPatTausForMuTauEta23,
      selectedPatTausForMuTauPt20,
      selectedPatTausForMuTauLeadTrk,
      selectedPatTausForMuTauLeadTrkPt,
      selectedPatTausForMuTauMuonVeto,
      selectedPatTausForMuTauCaloMuonVeto,
      selectedPatTausForMuTauElectronVeto,
      selectedPatTausForMuTauEcalCrackVeto ],
    src = "patTausLoosePFIsoEmbedded06",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectPatTausForMuTau = patTauSelConfiguratorForMuTau.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------
# define selection of Muon + Tau-jet candidate pairs
#--------------------------------------------------------------------------------

muTauPairsForTauIdEff = allMuTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatTausForMuTauEcalCrackVetoCumulative'),
    srcMET = cms.InputTag('patPFMETs'),
    srcGenParticles = cms.InputTag('')
)    

selectedMuTauPairsAntiOverlapVeto.cut = cms.string('dR12 > 0.7')
selectedMuTauPairsMt1MET.cut = cms.string('mt1MET < 40.')
selectedMuTauPairsPzetaDiff.cut = cms.string('(pZeta - 1.5*pZetaVis) > -20.')
selectedMuTauPairsZeroCharge.cut = cms.string('(leg1.charge + leg2.leadPFChargedHadrCand.charge) = 0')

patMuTauPairSelConfigurator = objSelConfigurator(
    [ selectedMuTauPairsAntiOverlapVeto,
      selectedMuTauPairsMt1MET,
      selectedMuTauPairsPzetaDiff,
      selectedMuTauPairsZeroCharge ],
    src = "muTauPairsForTauIdEff",
    pyModuleName = __name__,
    doSelIndividual = True
)

selectMuTauPairs = patMuTauPairSelConfigurator.configure(pyNameSpace = locals())

selectedMuTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedMuTauPairsAntiOverlapVetoCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                   
)

#--------------------------------------------------------------------------------
# define event selection sequence
#--------------------------------------------------------------------------------

muonPFTauSkimPath = cms.Path(
    hltMu
   + selectPrimaryVertex
   + patMuonsLoosePFIsoEmbedded + selectPatMuons
   + patTausLoosePFIsoEmbedded + selectPatTausForMuTau
   + muTauPairsForTauIdEff + selectMuTauPairs
   + selectedMuTauPairFilter
   + dataQualityFilters
)

tauIdEffSampleEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'muonPFTauSkimPath'
        )
    )
)
