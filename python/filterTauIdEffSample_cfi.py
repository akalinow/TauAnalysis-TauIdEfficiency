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
        triggerPaths = cms.vstring(
            # use IsoMu17 for 2011 Run A Data/MC                             
            'HLT_IsoMu17_v5',
            'HLT_IsoMu17_v6',
            'HLT_IsoMu17_v8',
            'HLT_IsoMu17_v9',
            'HLT_IsoMu17_v10',
            'HLT_IsoMu17_v11',
            'HLT_IsoMu17_v13',
            'HLT_IsoMu17_v14',
            # use Mu15/IsoMu15 for 2011 Run B MC                        
            'HLT_Mu15_v5',                 
            'HLT_IsoMu15_v5',
            # use Mu15_L1ETM20/IsoMu15_L1ETM20 for 2011 Run B Data (NOTE: events are in MET dataset)
            'HLT_IsoMu15_L1ETM20_v3',
            'HLT_IsoMu15_L1ETM20_v4',                 
            'HLT_Mu15_L1ETM20_v3',
            'HLT_Mu15_L1ETM20_v4'                 
        )
    )
)

#--------------------------------------------------------------------------------
# produce and select collections of Muon, Tau and Muon + Tau pair objects
#--------------------------------------------------------------------------------

# CV: done in TauAnalysis/TauIdEfficiency/test/commissioning/produceTauIdEffMeasNtuple_cfg.py

#--------------------------------------------------------------------------------
# select events containing at least one Muon + Tau pair passing selection
#--------------------------------------------------------------------------------

muonPFTauFixedConeFilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuPFTauFixedConePairsDzForTauIdEffCumulative'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000)                                   
)

muonPFTauShrinkingConeFilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuPFTauShrinkingConePairsDzForTauIdEffCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000)                                   
)

muonPFTauHPSfilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuPFTauHPSpairsDzForTauIdEffCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000)                                   
)

muonPFTauHPSpTaNCfilter = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('selectedMuPFTauHPSpTaNCpairsDzForTauIdEffCumulative'),      
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000)                                   
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

muonPFTauFixedConeSkimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonPFTauFixedConeFilter + countEventsPassed)
muonPFTauShrinkingConeSkimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonPFTauShrinkingConeFilter + countEventsPassed)
muonPFTauHPSskimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonPFTauHPSfilter + countEventsPassed)
muonPFTauHPSpTaNCskimPath = cms.Path(countEventsProcessed + commonSkimSequence + muonPFTauHPSpTaNCfilter + countEventsPassed)

tauIdEffSampleEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            #'muonPFTauFixedConeSkimPath',
            #'muonPFTauShrinkingConeSkimPath',
            'muonPFTauHPSskimPath',
            #'muonPFTauHPSpTaNCskimPath'
        )
    )
)
