import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from TauAnalysis.RecoTools.recoVertexSelection_cff import *
from TauAnalysis.Skimming.goldenZmmSelectionVBTFnoMuonIsolation_cfi import *

#--------------------------------------------------------------------------------
# select Z --> mu+ mu- events
# from which to determine tau fake-rates
#--------------------------------------------------------------------------------

zmmHLTFilter.HLTPaths = [
    'HLT_IsoMu17_v5', # Summer'11 MC
    'HLT_IsoMu17_v6',
    'HLT_IsoMu17_v8',
    'HLT_IsoMu17_v9',
    'HLT_IsoMu17_v11',
    'HLT_IsoMu17_v12',
    'HLT_IsoMu17_v13',
    'HLT_IsoMu17_v14',
    'HLT_DoubleMu7_v1', # Summer'11 MC
    'HLT_Mu13_Mu8_v7',
    'HLT_Mu13_Mu8_v10',
    'HLT_Mu17_Mu8_v7',
    'HLT_Mu17_Mu8_v10'
]

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
    src = cms.InputTag('hpsPFTauProducer'),
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
            ##'diMuonCaloTauSkimPath',
            'diMuonPFTauSkimPath'
        )
    )
)
