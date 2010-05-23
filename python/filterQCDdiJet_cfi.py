import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.filterDataQuality_cfi import *

#--------------------------------------------------------------------------------
# select QCD di-jet (/multi-jet) events
# from which to determine tau fake-rates
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# define HLT trigger path
#
hltJet15U = cms.EDFilter("EventSelPluginFilter",
    selector = cms.PSet(
        pluginName = cms.string('hltJet15U'),             
        pluginType = cms.string('TriggerResultEventSelector'),
        src = cms.InputTag('TriggerResults::HLT'),
        triggerPaths = cms.vstring('HLT_Jet15U')
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
    cut = cms.string("abs(caloTauTagInfoRef().calojetRef().eta) < 2.5 & caloTauTagInfoRef().calojetRef().pt > 10."),
    filter = cms.bool(True)
)

caloTauPairs = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedCaloTaus'),
    srcLeg2 = cms.InputTag('selectedCaloTaus'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                  
    verbosity = cms.untracked.int32(0)                                       
)

selectedCaloTauPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('caloTauPairs'),
    cut = cms.string("dR12 > 0.3"),
    filter = cms.bool(True)                                     
)

caloTauSkimPath = cms.Path(
   ## hltJet15U
   ##+
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
    cut = cms.string("abs(pfTauTagInfoRef().pfjetRef().eta) < 2.5 & pfTauTagInfoRef().pfjetRef().pt > 10."),
    filter = cms.bool(True)
)

pfTauPairs = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPFTaus'),
    srcLeg2 = cms.InputTag('selectedPFTaus'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                  
    verbosity = cms.untracked.int32(0)                                       
)

selectedPFTauPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('pfTauPairs'),
    cut = cms.string("dR12 > 0.3"),
    filter = cms.bool(True)                                     
)

pfTauSkimPath = cms.Path(    
   ## hltJet15U
   ##+
    selectedPFTaus + pfTauPairs + selectedPFTauPairs
   + dataQualityFilters
)
#--------------------------------------------------------------------------------

qcdDiJetEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'caloTauSkimPath',
            'pfTauSkimPath'
        )
    )
)
