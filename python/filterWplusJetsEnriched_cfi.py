import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from TauAnalysis.RecoTools.recoVertexSelection_cff import *

#--------------------------------------------------------------------------------
# select W + jets events
# from which to determine tau fake-rates
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# define HLT trigger path
#
hltMu = cms.EDFilter("EventSelPluginFilter",
    selector = cms.PSet(
        pluginName = cms.string('hltMu'),             
        pluginType = cms.string('TriggerResultEventSelector'),
        src = cms.InputTag('TriggerResults::HLT'),
        triggerPaths = cms.vstring(
            'HLT_IsoMu17_v5', # Summer'11 MC
            'HLT_IsoMu17_v6',
            'HLT_IsoMu17_v8',
            'HLT_IsoMu17_v9',
            'HLT_IsoMu17_v11',
            'HLT_IsoMu17_v12',
            'HLT_IsoMu17_v13',
            'HLT_IsoMu17_v14',
            'HLT_IsoMu24_v1',
            'HLT_IsoMu24_v2',
            'HLT_IsoMu24_v4',
            'HLT_IsoMu24_v5',
            'HLT_IsoMu24_v6',
            'HLT_IsoMu24_v7',
            'HLT_IsoMu24_v8',
            'HLT_IsoMu24_v9',
            'HLT_IsoMu24_v12',
            'HLT_IsoMu24_v13'
        )
    )
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# define Vertex selection
#
goodVertex = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("isValid & ndof >= 4 & abs(z) < 24 & abs(position.Rho) < 2"),
    filter = cms.bool(True)
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# define selection of "loose" Muons
#
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer

from CommonTools.ParticleFlow.pfNoPileUp_cff import *
pfPileUp.Enable = cms.bool(True)
pfPileUp.checkClosestZVertex = cms.bool(True)

from CommonTools.ParticleFlow.pfParticleSelection_cff import *

from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import patMuons

# compute muon IsoDeposits and add muon isolation sums to pat::Muon objects
from RecoMuon.MuonIsolation.muonPFIsolation_cff import *
import PhysicsTools.PatAlgos.tools.helpers as patutils
patutils.massSearchReplaceAnyInputTag(muonPFIsolationDepositsSequence, cms.InputTag('muons1stStep'), cms.InputTag('muons'))
patMuons.isoDeposits = cms.PSet(
    # CV: strings for IsoDeposits defined in PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc
    pfChargedHadrons = cms.InputTag("muPFIsoDepositCharged"),
    pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutral"),
    pfPhotons = cms.InputTag("muPFIsoDepositGamma"),
    user = cms.VInputTag(
        cms.InputTag("muPFIsoDepositChargedAll"),
        cms.InputTag("muPFIsoDepositPU")
    )
)

patMuons.userIsolation = cms.PSet(
    # CV: strings for Isolation values defined in PhysicsTools/PatAlgos/src/MultiIsolator.cc
    pfChargedHadron = cms.PSet(
        deltaR = cms.double(0.4),
        src = patMuons.isoDeposits.pfChargedHadrons,
        vetos = muPFIsoValueCharged04.deposits[0].vetos,
        skipDefaultVeto = muPFIsoValueCharged04.deposits[0].skipDefaultVeto
    ),
    pfNeutralHadron = cms.PSet(
        deltaR = cms.double(0.4),
        src = patMuons.isoDeposits.pfNeutralHadrons,
        vetos = muPFIsoValueNeutral04.deposits[0].vetos,
        skipDefaultVeto = muPFIsoValueNeutral04.deposits[0].skipDefaultVeto
    ),
    pfGamma = cms.PSet(
        deltaR = cms.double(0.4),
        src = patMuons.isoDeposits.pfPhotons,
        vetos = muPFIsoValueGamma04.deposits[0].vetos,
        skipDefaultVeto = muPFIsoValueGamma04.deposits[0].skipDefaultVeto
    ),
    user = cms.VPSet(
        cms.PSet(
            deltaR = cms.double(0.4),
            src = patMuons.isoDeposits.user[0],
            vetos = muPFIsoValueChargedAll04.deposits[0].vetos,
            skipDefaultVeto = muPFIsoValueChargedAll04.deposits[0].skipDefaultVeto
        ),
        cms.PSet(
            deltaR = cms.double(0.4),
            src = patMuons.isoDeposits.user[1],
            vetos = muPFIsoValuePU04.deposits[0].vetos,
            skipDefaultVeto = muPFIsoValuePU04.deposits[0].skipDefaultVeto
        )
    )
)

patMuons.addGenMatch = cms.bool(False)
patMuons.embedHighLevelSelection = cms.bool(True)
patMuons.usePV = cms.bool(False) # compute transverse impact parameter wrt. beamspot (not event vertex)

# Cuts for both muons, no isolation cuts applied
selectedMuons = cms.EDFilter("PATMuonSelector",
  src = cms.InputTag("patMuons"),
    cut = cms.string(
      'pt > 20 & abs(eta) < 2.5 & isGlobalMuon' \
     + ' & innerTrack.hitPattern.numberOfValidTrackerHits > 9 & innerTrack.hitPattern.numberOfValidPixelHits > 0' \
     + ' & abs(dB) < 0.2 & globalTrack.normalizedChi2 < 10' \
     + ' & globalTrack.hitPattern.numberOfValidMuonHits > 0 & numberOfMatches > 1' 
    ),
    filter = cms.bool(True)
)

# Cuts for muon leg with isolation cut applied
selectedIsoMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("selectedMuons"),
    cut = cms.string(
        '(userIsolation("pat::User1Iso")' + \
        ' + max(0., userIsolation("pat::PfNeutralHadronIso") + userIsolation("pat::PfGammaIso")' + \
        '          - 0.5*userIsolation("pat::User2Iso"))) < 0.06*pt'
    ),                                
    filter = cms.bool(False)
)

globalMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag('patMuons'),
    cut = cms.string("isGlobalMuon"),
    filter = cms.bool(False)
)

diMuonVeto = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('globalMuons'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# define loose CaloTau candidate/CaloJet selection
#
selectedCaloTaus = cms.EDFilter("CaloTauSelector",
    src = cms.InputTag('caloRecoTauProducer'),
    discriminators = cms.VPSet(),
    cut = cms.string("abs(caloTauTagInfoRef().jetRef().eta) < 2.5 & caloTauTagInfoRef().jetRef().pt > 10."),
    filter = cms.bool(False)
)

muonCaloTauPairs = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedIsoMuons'),
    srcLeg2 = cms.InputTag('selectedCaloTaus'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag('metJESCorAK5CaloJetMuons'),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                  
    verbosity = cms.untracked.int32(0)                                       
)

selectedMuonCaloTauPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('muonCaloTauPairs'),
    cut = cms.string("dR12 > 0.7 & mt1MET > 50."),
    filter = cms.bool(False)                                     
)

produceMuonCaloTauPairs = cms.Sequence(
    goodVertex
   * pfNoPileUpSequence
   * pfParticleSelectionSequence
   * muonPFIsolationDepositsSequence
   * patMuons * selectedMuons * selectedIsoMuons
   * selectedCaloTaus * muonCaloTauPairs * selectedMuonCaloTauPairs
)

selectedMuonCaloTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedMuonCaloTauPairs'),
    minNumber = cms.uint32(1)
)

muonCaloTauSkimPath = cms.Path(
    hltMu
   + produceMuonCaloTauPairs
   + globalMuons + diMuonVeto + selectedMuonCaloTauPairFilter
   + dataQualityFilters
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# define loose PFTau candidate/PFJet selection
#
selectedPFTaus = cms.EDFilter("PFTauSelector",
    src = cms.InputTag('shrinkingConePFTauProducer'),
    discriminators = cms.VPSet(),                          
    cut = cms.string("abs(jetRef().eta) < 2.5 & jetRef().pt > 10."),
    filter = cms.bool(False)
)

muonPFTauPairs = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedIsoMuons'),
    srcLeg2 = cms.InputTag('selectedPFTaus'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag('pfMet'),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)                                       
)

selectedMuonPFTauPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('muonPFTauPairs'),
    cut = cms.string("dR12 > 0.7 & mt1MET > 50."),
    filter = cms.bool(False)                                     
)

produceMuonPFTauPairs = cms.Sequence(
    goodVertex
   * pfNoPileUpSequence
   * pfParticleSelectionSequence
   * muonPFIsolationDepositsSequence
   * patMuons * selectedMuons * selectedIsoMuons
   * selectedPFTaus * muonPFTauPairs * selectedMuonPFTauPairs
)

selectedMuonPFTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedMuonPFTauPairs'),
    minNumber = cms.uint32(1)
)

muonPFTauSkimPath = cms.Path(    
    hltMu
   + produceMuonPFTauPairs
   + globalMuons + diMuonVeto + selectedMuonPFTauPairFilter
   + dataQualityFilters
)
#--------------------------------------------------------------------------------

wPlusJetsEnrichedEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            ##'muonCaloTauSkimPath',
            'muonPFTauSkimPath'
        )
    )
)
