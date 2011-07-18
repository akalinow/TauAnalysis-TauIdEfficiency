import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *
from TauAnalysis.RecoTools.recoVertexSelection_cff import *

#--------------------------------------------------------------------------------
# select QCD di-jet (/multi-jet) events
# from which to determine tau fake-rates
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# define HLT trigger path
#
hltSingleJet = cms.EDFilter("EventSelPluginFilter",
    selector = cms.PSet(
        pluginName = cms.string('hltSingleJet'),             
        pluginType = cms.string('TriggerResultEventSelector'),
        src = cms.InputTag('TriggerResults::HLT'),
        triggerPaths = cms.vstring(
            'HLT_Jet30_v1',  'HLT_Jet30_v2',  'HLT_Jet30_v3',  'HLT_Jet30_v4',  'HLT_Jet30_v5',  'HLT_Jet30_v6', 
            'HLT_Jet60_v1',  'HLT_Jet60_v2',  'HLT_Jet60_v3',  'HLT_Jet60_v4',  'HLT_Jet60_v5',  'HLT_Jet60_v6',
            'HLT_Jet80_v1',  'HLT_Jet80_v2',  'HLT_Jet80_v3',  'HLT_Jet80_v4',  'HLT_Jet80_v5',  'HLT_Jet80_v6', 
            'HLT_Jet110_v1', 'HLT_Jet110_v2', 'HLT_Jet110_v3', 'HLT_Jet110_v4', 'HLT_Jet110_v5', 'HLT_Jet110_v6',
            'HLT_Jet150_v1', 'HLT_Jet150_v2', 'HLT_Jet150_v3', 'HLT_Jet150_v4', 'HLT_Jet150_v5', 'HLT_Jet150_v6'
        )
    )
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
    filter = cms.bool(True)
)

caloTauPairs = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedCaloTaus'),
    srcLeg2 = cms.InputTag('selectedCaloTaus'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)                                       
)

selectedCaloTauPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('caloTauPairs'),
    cut = cms.string("dR12 > 0.3"),
    filter = cms.bool(True)                                     
)

caloTauSkimPath = cms.Path(
    selectedCaloTaus + caloTauPairs + selectedCaloTauPairs
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
    filter = cms.bool(True)
)

pfTauPairs = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPFTaus'),
    srcLeg2 = cms.InputTag('selectedPFTaus'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)                                       
)

selectedPFTauPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('pfTauPairs'),
    cut = cms.string("dR12 > 0.3"),
    filter = cms.bool(True)                                     
)

pfTauSkimPath = cms.Path(    
    selectedPFTaus + pfTauPairs + selectedPFTauPairs
   + dataQualityFilters
)
#--------------------------------------------------------------------------------

qcdDiJetEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            ##'caloTauSkimPath',
            'pfTauSkimPath'
        )
    )
)
