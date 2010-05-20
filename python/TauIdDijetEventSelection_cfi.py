import FWCore.ParameterSet.Config as cms

'''

TauIdDijetEventSelection_cfi

Apply event level skimmming cuts.

Requirements:
    Jet trigger > 15 @ HLT
    At least two jets with |eta| < 2.5 & pt > 5
    At least one jet must satisfy the HLT trigger

TODO: add equivalent sequence for CaloTaus

'''

# The basic HLT requirement
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
jetTrigger = hltHighLevel.clone()
jetTrigger.HLTPaths = cms.vstring('HLT_Jet15U')

# Only look at PFJets relevant to taus
pfJetsInAcceptance = cms.EDProducer(
    'CandViewRefSelector',
    src = cms.InputTag('iterativeCone5PFJets'),
    cut = cms.string('abs(eta) < 2.5 && pt() > 5'),
)

# We must have at least two jets (one tag, >= 1 probe)
atLeastTwoPFJetsInAcceptance = cms.EDFilter(
    'CandViewCountFilter',
    src = cms.InputTag('pfJetsInAcceptance'),
    minNumber = cms.uint32(2)
)

from TauAnalysis.TauIdEfficiency.TauIdDijetTagAndProbeProducer_cfi \
        import dijetTagAndProbes

# Build tag and probe collections from jets
pfJetsTagAndProbes = dijetTagAndProbes.clone()
pfJetsTagAndProbes.source = cms.InputTag('pfJetsInAcceptance')
pfJetsTagAndProbes.pluginType = cms.string('CandidateVectorValExtractor')
pfJetsTagAndProbes.expression = cms.string('pt()')
pfJetsTagAndProbes.triggerThreshold = cms.double(15.0) # fixme ?

# Require at least one tag jet
atLeastOneTagPFJet = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("pfJetsTagAndProbes", "tagObject"),
    minNumber = cms.uint32(1),
)

# Sequence to produce PFJets
pfJetTagAndProbeSequence = cms.Sequence(
    jetTrigger
    * pfJetsInAcceptance
    * atLeastTwoPFJetsInAcceptance
    * pfJetsTagAndProbes
    * atLeastOneTagPFJet
)

