import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from TauAnalysis.RecoTools.recoVertexSelection_cff import *
from TauAnalysis.Skimming.goldenZmmSelectionVBTFnoMuonIsolation_cfi import *

#--------------------------------------------------------------------------------
# select Z --> mu+ mu- events
# from which to determine tau fake-rates
#--------------------------------------------------------------------------------

zmmHLTFilter.HLTPaths = [ 'HLT_Mu15_v1', 'HLT_Mu15_v2', 'HLT_IsoMu15_v5', 'HLT_IsoMu17_v5', 'HLT_Mu24_v2',
                          'HLT_DoubleMu6_v1', 'HLT_DoubleMu6_v2', 'HLT_DoubleMu7_v1', 'HLT_DoubleMu7_v2' ]

goldenZmumuCandidateFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('goldenZmumuCandidatesGe2IsoMuons'),
    minNumber = cms.uint32(1)
)

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

selectedCaloTausAntiOverlapWithLeptonsVeto = cms.EDFilter("CaloTauAntiOverlapSelector",
    src = cms.InputTag('selectedCaloTaus'),                                                      
    srcNotToBeFiltered = cms.VInputTag('goodIsoMuons'),
    dRmin = cms.double(0.7),
    filter = cms.bool(False)
)

selectedCaloTauAntiOverlapWithLeptonsFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedCaloTausAntiOverlapWithLeptonsVeto'),
    minNumber = cms.uint32(1)
)

diMuonCaloTauSkimPath = cms.Path(
     pfNoPileUpSequence
   + goldenZmumuSelectionSequence + goldenZmumuCandidateFilter
   + selectedCaloTaus + selectedCaloTausAntiOverlapWithLeptonsVeto + selectedCaloTauAntiOverlapWithLeptonsFilter
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

selectedPFTausAntiOverlapWithLeptonsVeto = cms.EDFilter("PFTauAntiOverlapSelector",
    src = cms.InputTag('selectedPFTaus'),                                                      
    srcNotToBeFiltered = cms.VInputTag('goodIsoMuons'),
    dRmin = cms.double(0.7),
    filter = cms.bool(False)
)

selectedPFTauAntiOverlapWithLeptonsFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedPFTausAntiOverlapWithLeptonsVeto'),
    minNumber = cms.uint32(1)
)

diMuonPFTauSkimPath = cms.Path(
    pfNoPileUpSequence
   + goldenZmumuSelectionSequence + goldenZmumuCandidateFilter
   + selectedPFTaus + selectedPFTausAntiOverlapWithLeptonsVeto + selectedPFTauAntiOverlapWithLeptonsFilter
   + dataQualityFilters
)
#--------------------------------------------------------------------------------

zMuMuEnrichedEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'diMuonCaloTauSkimPath',
            'diMuonPFTauSkimPath'
        )
    )
)
