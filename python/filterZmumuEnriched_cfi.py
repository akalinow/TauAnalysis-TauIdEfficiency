import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from TauAnalysis.RecoTools.recoVertexSelection_cff import *
from TauAnalysis.Skimming.goldenZmmSelectionVBTFnoMuonIsolation_cfi import *

#--------------------------------------------------------------------------------
# select Z --> mu+ mu- events
# from which to determine tau fake-rates
#--------------------------------------------------------------------------------

goodIsoMuons.chargedHadronIso.dRisoCone = cms.double(0.4)
goodIsoMuons.neutralHadronIso.dRisoCone = cms.double(0.4)
goodIsoMuons.photonIso.dRisoCone = cms.double(0.4)
goodIsoMuons.sumPtMax = cms.double(0.15)

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
    srcLeg1 = cms.InputTag('goodIsoMuons'),
    srcLeg2 = cms.InputTag('selectedCaloTaus'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag('metJESCorAK5CaloJetMuons'),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                  
    verbosity = cms.untracked.int32(0)                                       
)

selectedMuonCaloTauPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('muonCaloTauPairs'),
    cut = cms.string("dR12 > 0.7"),
    filter = cms.bool(False)                                     
)

produceMuonCaloTauPairs = cms.Sequence(
    patMuons * goodMuons * goodIsoMuons
   * selectedCaloTaus * muonCaloTauPairs * selectedMuonCaloTauPairs
)

selectedMuonCaloTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedMuonCaloTauPairs'),
    minNumber = cms.uint32(1)
)

muonCaloTauSkimPath = cms.Path(
    zmmHLTFilter
   + pfNoPileUpSequence
   + produceMuonCaloTauPairs
   + selectedMuonCaloTauPairFilter
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
    cut = cms.string("abs(pfTauTagInfoRef().pfjetRef().eta) < 2.5 & pfTauTagInfoRef().pfjetRef().pt > 10."),
    filter = cms.bool(False)
)

muonPFTauPairs = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('goodIsoMuons'),
    srcLeg2 = cms.InputTag('selectedPFTaus'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag('pfMet'),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)                                       
)

selectedMuonPFTauPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('muonPFTauPairs'),
    cut = cms.string("dR12 > 0.7"),
    filter = cms.bool(False)                                     
)

produceMuonPFTauPairs = cms.Sequence(
    patMuons * goodMuons * goodIsoMuons
   * selectedPFTaus * muonPFTauPairs * selectedMuonPFTauPairs
)

selectedMuonPFTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedMuonPFTauPairs'),
    minNumber = cms.uint32(1)
)

muonPFTauSkimPath = cms.Path(
    zmmHLTFilter
   + pfNoPileUpSequence
   + produceMuonPFTauPairs
   + selectedMuonPFTauPairFilter
   + dataQualityFilters
)
#--------------------------------------------------------------------------------

zMuMuEnrichedEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'muonCaloTauSkimPath',
            'muonPFTauSkimPath'
        )
    )
)
