import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from TauAnalysis.RecoTools.recoVertexSelection_cff import *
from TauAnalysis.Skimming.goldenZmmSelectionVBTFnoMuonIsolation_cfi import *

#--------------------------------------------------------------------------------
# select Z --> mu+ mu- events
# from which to determine tau fake-rates
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patJetSelection_cff import *
#
# define loose CaloTau candidate/CaloJet selection
#
selectedCaloTaus = cms.EDFilter("CaloTauSelector",
    src = cms.InputTag('caloRecoTauProducer'),
    discriminators = cms.VPSet(),
    cut = cms.string("abs(caloTauTagInfoRef().jetRef().eta) < 2.5 & caloTauTagInfoRef().jetRef().pt > 10."),
    filter = cms.bool(False)
)

selectedPatJetsAntiOverlapWithLeptonsVeto = cms.EDFilter("PATJetAntiOverlapSelector",
  srcNotToBeFiltered = cms.VInputTag( "goodIsoMuons") ,
 dRmin = cms.double(0.7),
 filter = cms.bool(False)
 )

produceDiMuonCaloTauPairs = cms.Sequence(
    goldenZmumuSelectionSequence
   * selectedCaloTaus * selectedPatJetsAntiOverlapWithLeptonsVeto
)

selectedDiMuonCaloTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedDiMuonCaloTauPairs'),
    minNumber = cms.uint32(1)
)

dimuonCaloTauSkimPath = cms.Path(
     pfNoPileUpSequence
   + produceDiMuonCaloTauPairs
   + selectedDiMuonCaloTauPairFilter
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

selectedPatJetsAntiOverlapWithLeptonsVeto = cms.EDFilter("PATJetAntiOverlapSelector",
  srcNotToBeFiltered = cms.VInputTag( "goodIsoMuons") ,
 dRmin = cms.double(0.7),
 filter = cms.bool(False)
 )

produceDiMuonPFTauPairs = cms.Sequence(
    goldenZmumuSelectionSequence
   * selectedPFTaus * selectedPatJetsAntiOverlapWithLeptonsVeto
)

selectedDiMuonPFTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedDiMuonPFTauPairs'),
    minNumber = cms.uint32(1)
)

dimuonPFTauSkimPath = cms.Path(
    pfNoPileUpSequence
   + produceDiMuonPFTauPairs
   + selectedDiMuonPFTauPairFilter
   + dataQualityFilters
)
#--------------------------------------------------------------------------------

zMuMuEnrichedEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'dimuonCaloTauSkimPath',
            'dimuonPFTauSkimPath'
        )
    )
)
