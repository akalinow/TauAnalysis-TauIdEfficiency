import FWCore.ParameterSet.Config as cms

'''
TauFakeRateDijetProbesProducer

Produce two unbiased collections corresponding the highest and second highest
PT objects from input collection <source> in the event.  The trigger selection
bias is removed by requiring an additional object satisfy the trigger requirement.

This producer produces three collections of RefToBaseVector<Candidate>, 
with the labels:

    highestPtProbe,
    secondHighestPtProbe,
    tagObject

Author: Evan Friis, UC Davis (friis@physics.ucdavis.edu)

'''

dijetTagAndProbes = cms.EDProducer(
    "TauIdDijetTagAndProbeProducer",
    source = cms.InputTag("probeTauJets"),
    pluginType = cms.string("PATTauVectorValExtractor"),
    expression = cms.string("pfTauTagInfoRef().pfjetRef().pt()"),
    triggerThreshold = cms.double(15.0),
    minimumThreshold = cms.double(1.0)
)
