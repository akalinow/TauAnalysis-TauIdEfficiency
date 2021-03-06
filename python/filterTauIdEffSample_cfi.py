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
            # use Mu15/IsoMu15 for 2012 Run A Data/MC
            'HLT_IsoMu15_L1ETM20_v5',                
            'HLT_IsoMu15_eta2p1_L1ETM20_v1',
            'HLT_IsoMu15_eta2p1_L1ETM20_v3',
            'HLT_IsoMu15_eta2p1_L1ETM20_v4',
            'HLT_IsoMu15_eta2p1_L1ETM20_v5'
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

selectedEventsSkimmingPFTauFixedCone = cms.EDProducer("EventCountProducer")
selectedEventsSkimmingPFTauShrinkingCone = cms.EDProducer("EventCountProducer")
selectedEventsSkimmingPFTauHPS = cms.EDProducer("EventCountProducer")
selectedEventsSkimmingPFTauHPSpTaNC = cms.EDProducer("EventCountProducer")

muonPFTauFixedConeSkimPath = cms.Path(commonSkimSequence + muonPFTauFixedConeFilter + selectedEventsSkimmingPFTauFixedCone)
muonPFTauShrinkingConeSkimPath = cms.Path(commonSkimSequence + muonPFTauShrinkingConeFilter + selectedEventsSkimmingPFTauShrinkingCone)
muonPFTauHPSskimPath = cms.Path(commonSkimSequence + muonPFTauHPSfilter + selectedEventsSkimmingPFTauHPS)
muonPFTauHPSpTaNCskimPath = cms.Path(commonSkimSequence + muonPFTauHPSpTaNCfilter + selectedEventsSkimmingPFTauHPSpTaNC)

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
