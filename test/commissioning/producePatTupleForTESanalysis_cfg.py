import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('producePatTupleForTESanalysis')

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START311_V2::All')

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/k/khurana/public/For-Christian/pickevents_1_1_nV0.root'
    )
    #skipBadFiles = cms.untracked.bool(True) 
)

#--------------------------------------------------------------------------------
# run Vertex reconstruction via Deterministic annealing algorithm
#
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
process.offlinePrimaryVerticesDAwithBS = process.offlinePrimaryVerticesDA.clone()
process.offlinePrimaryVerticesDAwithBS.useBeamConstraint = cms.bool(True)
process.offlinePrimaryVerticesDAwithBS.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
process.offlinePrimaryVerticesDAwithBS.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# rerun PFTau reconstruction
#
process.load('RecoTauTag/Configuration/RecoPFTauTag_cff')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce PAT objects
#
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
process.cleanPatTaus.preselection = cms.string(
    "pt > 15 & abs(eta) < 2.3 & tauID('decayModeFinding') > 0.5 & tauID('byLooseIsolation') > 0.5" \
  + " & tauID('againstMuonTight') > 0.5 & tauID('againstElectronLoose') > 0.5"
)
#--------------------------------------------------------------------------------

process.p = cms.Path(
    process.offlinePrimaryVerticesDAwithBS
   + process.PFTau
   + process.patDefaultSequence
)

# use offlinePrimaryVerticesDAwithBS Vertex collection as input for PFTau reconstruction
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag 
massSearchReplaceAnyInputTag(process.p, cms.InputTag("offlinePrimaryVerticesDA"), cms.InputTag('offlinePrimaryVerticesDAwithBS'))

process.patTupleEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'drop *',
        # keep pat::Electrons, Muons and Taus
        # ( cf. PhysicsTools/PatAlgos/python/patEventContent_cff.py )
        #'keep *_cleanPatElectrons*_*_*',
        'keep *_cleanPatMuons*_*_*',
        'keep *_cleanPatTaus*_*_*',
        # keep generator level information 
        'keep *_genParticles*_*_*',
        # keep all PFCandidates
        # ( cf. RecoParticleFlow/Configuration/python/RecoParticleFlow_EventContent_cff.py )
        'keep recoPFCandidates_particleFlow_*_*',
        # keep reco::PFJets reconstructed by anti-kT algorithm with dR = 0.5
        # ( cf. RecoJets/Configuration/python/RecoJets_EventContent_cff.py )
        'keep *_ak5PFJets_*_*',
        # keep reco::PFTau objects reconstructed by HPS algorithm,
        # pi0 candidates and tau id. discriminators
        # ( cf. RecoTauTag/Configuration/python/RecoTauTag_EventContent_cff.py )
        #'keep *_ak5PFJetsRecoTauPiZeros_*_*',
        #'keep *_hpsPFTauProducer*_*_*',
        #'keep *_hpsPFTauDiscrimination*_*_*'
    )
)

# drop unnecessary pat::Electron data-members
process.patElectrons.addEfficiencies = cms.bool(False)
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.addGenMatch = cms.bool(True)
process.patElectrons.addResolutions = cms.bool(False)
process.patElectrons.efficiencies = cms.PSet()
process.patElectrons.embedGenMatch = cms.bool(True)
process.patElectrons.embedGsfElectronCore = cms.bool(True)
process.patElectrons.embedGsfTrack = cms.bool(True)
process.patElectrons.embedHighLevelSelection = cms.bool(False)
process.patElectrons.embedPFCandidate = cms.bool(False)
process.patElectrons.embedSuperCluster = cms.bool(False)
process.patElectrons.embedTrack = cms.bool(False)
process.patElectrons.isoDeposits = cms.PSet()
process.patElectrons.resolutions = cms.PSet()
process.patElectrons.userIsolation = cms.PSet()

# drop unnecessary pat::Muon data-members
process.patMuons.addEfficiencies = cms.bool(False)
process.patMuons.addGenMatch = cms.bool(True)
process.patMuons.addResolutions = cms.bool(False)
process.patMuons.efficiencies = cms.PSet()
process.patMuons.embedCaloMETMuonCorrs = cms.bool(False)
process.patMuons.embedCombinedMuon = cms.bool(True)
process.patMuons.embedGenMatch = cms.bool(True)
process.patMuons.embedHighLevelSelection = cms.bool(False)
process.patMuons.embedPFCandidate = cms.bool(False)
process.patMuons.embedPickyMuon = cms.bool(False)
process.patMuons.embedStandAloneMuon = cms.bool(True)
process.patMuons.embedTcMETMuonCorrs = cms.bool(False)
process.patMuons.embedTpfmsMuon = cms.bool(False)
process.patMuons.embedTrack = cms.bool(True)
process.patMuons.isoDeposits = cms.PSet()
process.patMuons.resolutions = cms.PSet()
process.patMuons.userIsolation = cms.PSet()

# drop unnecessary pat::Tau data-members
process.patTaus.addEfficiencies = cms.bool(False)
process.patTaus.addGenJetMatch = cms.bool(True)
process.patTaus.addGenMatch = cms.bool(True)
process.patTaus.addResolutions = cms.bool(False)
process.patTaus.addTauID = cms.bool(True)
process.patTaus.efficiencies = cms.PSet()
process.patTaus.embedGenJetMatch = cms.bool(True)
process.patTaus.embedGenMatch = cms.bool(True)
process.patTaus.embedIsolationPFChargedHadrCands = cms.bool(False)
process.patTaus.embedIsolationPFGammaCands = cms.bool(False)
process.patTaus.embedIsolationPFNeutralHadrCands = cms.bool(False)
process.patTaus.embedIsolationTracks = cms.bool(False)
process.patTaus.embedLeadTrack = cms.bool(False)
process.patTaus.embedLeadPFChargedHadrCand = cms.bool(False)
process.patTaus.embedLeadPFCand = cms.bool(False)
process.patTaus.embedLeadPFNeutralCand = cms.bool(False)
process.patTaus.embedSignalPFCands = cms.bool(False)
process.patTaus.embedSignalPFChargedHadrCands = cms.bool(False)
process.patTaus.embedSignalPFGammaCands = cms.bool(False)
process.patTaus.embedSignalPFNeutralHadrCands = cms.bool(False)
process.patTaus.embedSignalTracks = cms.bool(False)
process.patTaus.isoDeposits = cms.PSet()
process.patTaus.resolutions = cms.PSet()
process.patTaus.userIsolation = cms.PSet()

process.patTupleOutputModule = cms.OutputModule("PoolOutputModule",
    process.patTupleEventContent,
    fileName = cms.untracked.string("patTupleForTESanalysis.root")
)

process.o = cms.EndPath(process.patTupleOutputModule)

# print-out all python configuration parameter information
#print process.dumpPython()

