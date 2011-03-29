import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from CommonTools.ParticleFlow.pfNoPileUp_cff import *
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

patMuonsLoosePFIsoEmbedded = cms.EDProducer("PATMuonPFIsolationEmbedder",
    src = cms.InputTag('cleanPatMuons'),                                       
    userFloatName = cms.string('pfLooseIsoPt'),
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = cms.PSet(
        ptMin = cms.double(0.5),        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.6)
    ),
    neutralHadronIso = cms.PSet(
        ptMin = cms.double(1000.),        
        dRvetoCone = cms.double(0.08),        
        dRisoCone = cms.double(0.6)
    ),
    photonIso = cms.PSet(
        ptMin = cms.double(0.5),        
        dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
        dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.6)
    )
)

selectedPatMuonsGlobal.cut = cms.string('isGlobalMuon()')
selectedPatMuonsEta21.cut = cms.string('abs(eta) < 2.1')
selectedPatMuonsPt15.cut = cms.string('pt > 15.')
selectedPatMuonsVbTfId.beamSpotSource = cms.InputTag("offlineBeamSpot")
selectedPatMuonsPFRelIso.chargedHadronIso.ptMin = patMuonsLoosePFIsoEmbedded.chargedHadronIso.ptMin
selectedPatMuonsPFRelIso.chargedHadronIso.dRvetoCone = patMuonsLoosePFIsoEmbedded.chargedHadronIso.dRvetoCone
selectedPatMuonsPFRelIso.chargedHadronIso.dRisoCone = patMuonsLoosePFIsoEmbedded.chargedHadronIso.dRisoCone
selectedPatMuonsPFRelIso.neutralHadronIso.ptMin = patMuonsLoosePFIsoEmbedded.neutralHadronIso.ptMin 
selectedPatMuonsPFRelIso.neutralHadronIso.dRvetoCone = patMuonsLoosePFIsoEmbedded.neutralHadronIso.dRvetoCone
selectedPatMuonsPFRelIso.neutralHadronIso.dRisoCone = patMuonsLoosePFIsoEmbedded.neutralHadronIso.dRisoCone
selectedPatMuonsPFRelIso.photonIso.ptMin = patMuonsLoosePFIsoEmbedded.photonIso.ptMin
selectedPatMuonsPFRelIso.photonIso.dRvetoCone = patMuonsLoosePFIsoEmbedded.photonIso.dRvetoCone
selectedPatMuonsPFRelIso.photonIso.dRisoCone = patMuonsLoosePFIsoEmbedded.photonIso.dRisoCone
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
    src = "patMuonsLoosePFIsoEmbedded",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectPatMuons = patMuonSelConfigurator.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------
# define loose Tau-jet candidate selection
#--------------------------------------------------------------------------------

patTausLoosePFIsoEmbedded = cms.EDProducer("PATTauPFIsolationEmbedder",
    src = cms.InputTag('patTaus'),                                       
    userFloatName = cms.string('pfLooseIsoPt'),
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = cms.PSet(
        ptMin = cms.double(1.0),        
        dRvetoCone = cms.double(0.15),
        dRisoCone = cms.double(0.6)
    ),
    neutralHadronIso = cms.PSet(
        ptMin = cms.double(1000.),        
        dRvetoCone = cms.double(0.15),        
        dRisoCone = cms.double(0.6)
    ),
    photonIso = cms.PSet(
        ptMin = cms.double(1.5),        
        dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
        dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
        dRvetoCone = cms.double(0.15),
        dRisoCone = cms.double(0.6)
    )
)

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
    src = "patTausLoosePFIsoEmbedded",
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
