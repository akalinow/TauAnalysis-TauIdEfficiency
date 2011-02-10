import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProductionTauIdEffMeasSpecific import \
     configurePatTupleProductionTauIdEffMeasSpecific

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
# produce and select collections of Muon, Tau and Muon + Tau pair objects
#--------------------------------------------------------------------------------

# CV: done in TauAnalysis/TauIdEfficiency/test/commissioning/produceTauIdEffMeasNtuple_cfg.py

#--------------------------------------------------------------------------------
# select events containing at least one Muon + Tau pair passing selection
#--------------------------------------------------------------------------------

muonCaloTauFilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuCaloTauPairsForTauIdEffCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                   
)

muonPFTauFixedConeFilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuPFTauFixedConePairsForTauIdEffCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                   
)

muonPFTauShrinkingConeFilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuPFTauShrinkingConePairsForTauIdEffCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                   
)

muonPFTauHPSfilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuPFTauHPSpairsForTauIdEffCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                   
)

muonPFTauHPSpTaNCfilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                   
)

#--------------------------------------------------------------------------------
# veto events containing di-Muon pairs
# (the hypothesis being that the pair results from a Z --> mu+ mu- decay)
#--------------------------------------------------------------------------------

diMuPairZmumuHypothesisVeto = cms.EDFilter("PATCandViewMaxFilter",
    src = cms.InputTag('selectedDiMuPairForTauIdEffZmumuHypotheses'),
    maxNumber = cms.uint32(0)                                                            
)

#--------------------------------------------------------------------------------
# define event selection sequence
#--------------------------------------------------------------------------------

commonSkimSequence = cms.Sequence(
    hltMu
   * dataQualityFilters
   * diMuPairZmumuHypothesisVeto
)

DQMStore = cms.Service("DQMStore")

countEventsProcessed = cms.EDAnalyzer("DQMEventCounter",
    meName = cms.string('numEventsProcessed')                                         
)

countEventsPassed = cms.EDAnalyzer("DQMEventCounter",
    meName = cms.string('numEventsPassed')                                         
)

muonCaloTauSkimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonCaloTauFilter + countEventsPassed)
muonPFTauFixedConeSkimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonPFTauFixedConeFilter + countEventsPassed)
muonPFTauShrinkingConeSkimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonPFTauShrinkingConeFilter + countEventsPassed)
muonPFTauHPSskimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonPFTauHPSfilter + countEventsPassed)
muonPFTauHPSpTaNCskimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonPFTauHPSpTaNCfilter + countEventsPassed)

tauIdEffSampleEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'muonCaloTauSkimPath',
            'muonPFTauFixedConeSkimPath',
            'muonPFTauShrinkingConeSkimPath',
            'muonPFTauHPSskimPath',
            'muonPFTauHPSpTaNCskimPath'
        )
    )
)
