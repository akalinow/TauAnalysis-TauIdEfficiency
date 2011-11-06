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
# define selection of "loose" Muons
#
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
#
# define selection of "tight" Muons
#
selectedMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag('patMuons'),
    cut = cms.string("isGlobalMuon & pt > 20. & abs(eta) < 2.1"),
    filter = cms.bool(False)
)

selectedIsoMuons = cms.EDFilter("PATMuonPFIsolationSelector",
    src = cms.InputTag('selectedMuons'),
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = cms.PSet(
        ptMin = cms.double(1.0),        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.4)
    ),
    neutralHadronIso = cms.PSet(
        ptMin = cms.double(1.0),        
        dRvetoCone = cms.double(0.08),        
        dRisoCone = cms.double(0.4)
    ),
    photonIso = cms.PSet(
        ptMin = cms.double(1.0),        
        dPhiVeto = cms.double(-1.),
        dEtaVeto = cms.double(-1.),
        dRvetoCone = cms.double(0.05),
        dRisoCone = cms.double(0.4)
    ),
    sumPtMax = cms.double(0.10),
    sumPtMethod = cms.string("relative"),
    filter = cms.bool(True)                    
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
    globalMuons * selectedMuons * selectedIsoMuons
   * selectedCaloTaus * muonCaloTauPairs * selectedMuonCaloTauPairs
)

selectedMuonCaloTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedMuonCaloTauPairs'),
    minNumber = cms.uint32(1)
)

muonCaloTauSkimPath = cms.Path(
    hltMu
   + produceMuonCaloTauPairs
   + diMuonVeto + selectedMuonCaloTauPairFilter
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
    globalMuons * selectedMuons * selectedIsoMuons
   * selectedPFTaus * muonPFTauPairs * selectedMuonPFTauPairs
)

selectedMuonPFTauPairFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedMuonPFTauPairs'),
    minNumber = cms.uint32(1)
)

muonPFTauSkimPath = cms.Path(    
    hltMu
   + produceMuonPFTauPairs
   + diMuonVeto + selectedMuonPFTauPairFilter
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
