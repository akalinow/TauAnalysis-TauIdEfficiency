import FWCore.ParameterSet.Config as cms
import copy

'''

tools.sequenceBuilder.py

Create pat::Taus from a specified collection, match
them to the tag and probe jets, and create ntuples
of the resulting collections.

Author: Evan K. Friis, Christian Veelken (UC Davis)

'''

##import TauAnalysis.TauIdEfficiency.patConfiguration.tauProductionPrototypes as tauProto
##import TauAnalysis.TauIdEfficiency.patConfiguration.matchingPrototypes as matchProto

def buildTauSequence(
    process, 
    collectionName = [ "patTaus", "" ],
    patTauProducerPrototype = None,
    patTauCleanerPrototype = None,
    triggerMatcherProtoType = None,
    addGenInfo = False):
    '''
    blah
    '''

    # build basic pat::Taus from input collection
    patTauProducer = copy.deepcopy(patTauProducerPrototype)
    patTauProducerName = "".join(collectionName)
    setattr(process, patTauProducerName, patTauProducer)

    # configure matching of basic pat::Tau collection
    # to generator level particles and jets
    outputSequence = None
    if addGenInfo:
        patTauGenParticleMatch = process.tauMatch.clone(
            src = patTauProducer.tauSource
        )
        patTauGenParticleMatchName = collectionName[0] + "GenParticleMatch" + collectionName[1]
        setattr(process, patTauGenParticleMatchName, patTauGenParticleMatch)

        patTauGenJetMatch = process.tauGenJetMatch.clone(
           src = patTauProducer.tauSource
        )
        patTauGenJetMatchName = collectionName[0] + "GenJetMatch" + collectionName[1]
        setattr(process, patTauGenJetMatchName, patTauGenJetMatch)

        patTauProducer.genParticleMatch = cms.InputTag(patTauGenParticleMatchName)
        patTauProducer.genJetMatch = cms.InputTag(patTauGenJetMatchName)
    
        outputSequence = cms.Sequence(
            patTauGenParticleMatch + patTauGenJetMatch
           + patTauProducer
        )
    else:
        patTauProducer.addGenMatch = cms.bool(False)
        patTauProducer.addGenJetMatch = cms.bool(False)
        
        outputSequence = cms.Sequence(
            patTauProducer
        )

    # configure matching of basic pat::Tau collection
    # to trigger primitives;
    # produce new collection of pat::Taus with trigger primitives embedded
    patTauTriggerMatch = triggerMatcherProtoType.clone(
        src = cms.InputTag(patTauProducerName)
    )
    patTauTriggerMatchName = collectionName[0]  + "TriggerMatched" + collectionName[1]
    setattr(process, patTauTriggerMatchName, patTauTriggerMatch)
    outputSequence += patTauTriggerMatch

    patTauTriggerEvent = process.patTriggerEvent.clone(
        patTriggerMatches = cms.VInputTag(patTauTriggerMatchName)
    )
    patTauTriggerEventName = collectionName[0] + "TriggerEvent" + collectionName[1]
    setattr(process, patTauTriggerEventName, patTauTriggerEvent)
    outputSequence += patTauTriggerEvent

    patTauTriggerEmbedder = cms.EDProducer("PATTriggerMatchTauEmbedder",
        src     = cms.InputTag(patTauProducerName),
        matches = cms.VInputTag(patTauTriggerMatchName)
    )
    patTauTriggerEmbedderName = collectionName[0] + "TriggerEmbedder" + collectionName[1]
    setattr(process, patTauTriggerEmbedderName, patTauTriggerEmbedder)
    outputSequence += patTauTriggerEmbedder

    # configure PATTauCleaner module
    # for removal of tau-jet candidates "overlapping" with electrons or muons
    patTauCleaner = copy.deepcopy(patTauCleanerPrototype)
    patTauCleaner.src = cms.InputTag(patTauTriggerEmbedderName)
    patTauCleanerName = collectionName[0] + "Cleaned" + collectionName[1]
    setattr(process, patTauCleanerName, patTauCleaner)
    outputSequence += patTauCleaner
    
    # return sequence for production of basic tau collection,
    # generator level particle and jet matches and trigger primitives embedded;
    # together with name of "cleaned" tau collection
    # to be used as InputTag for further processing
    retVal = {}
    retVal["sequence"] = outputSequence
    retVal["collection"] = patTauCleanerName
    return retVal

def buildQCDdiJetTauSequence(
    process, 
    collectionName = [ "patTaus", "" ],
    patTauProducerPrototype = None,
    patTauCleanerPrototype = None,
    triggerMatcherProtoType = None,
    addGenInfo = False):
    '''
    blah
    '''

    # produce collection of basic pat::Taus
    # matched to generator level particles and jets
    # and trigger primitives embedded
    retVal_tau = buildTauSequence(
        process, 
        collectionName = collectionName,
        patTauProducerPrototype = patTauProducerPrototype,
        patTauCleanerPrototype = patTauCleanerPrototype,
        triggerMatcherProtoType = triggerMatcherProtoType,
        addGenInfo = addGenInfo
    )

    outputSequence = retVal_tau["sequence"]

    # produce final collection of pat::Taus with
    # flags embedded to:
    #  o indicate index of tau-jet candidate in Pt-sorted collection of jets
    #  o indicate whether tau-jet candidate is tag/probe
    patTauDijetTagAndProbe = cms.EDProducer("TauIdTagAndProbeProducer",
        source = cms.InputTag(retVal_tau["collection"]),
        triggerPath = cms.string(triggerMatcherProtoType.pathNames.value()[0])
    )
    patTauDijetTagAndProbeName = collectionName[0] + "DijetTagAndProbe" + collectionName[1]
    setattr(process, patTauDijetTagAndProbeName, patTauDijetTagAndProbe)
    outputSequence += cms.Sequence(patTauDijetTagAndProbe)

    # return full sequence for production of basic tau collection,
    # generator level particle and jet matches, trigger matches,
    # Pt-sorting and tag/probe flags;
    # together with name of "cleaned" tau collection
    # to be used as InputTag for further processing
    retVal = {}
    retVal["sequence"] = outputSequence
    retVal["collection"] = patTauDijetTagAndProbeName
    return retVal

def buildQCDmuEnrichedTauSequence(
    process, 
    collectionName = [ "patTaus", "" ],
    patTauProducerPrototype = None,
    patTauCleanerPrototype = None,
    triggerMatcherProtoType = None,
    addGenInfo = False):
    '''
    blah
    '''

    # produce collection of basic pat::Taus
    # matched to generator level particles and jets
    # and trigger primitives embedded
    retVal_tau = buildTauSequence(
        process, 
        collectionName = collectionName,
        patTauProducerPrototype = patTauProducerPrototype,
        patTauCleanerPrototype = patTauCleanerPrototype,
        triggerMatcherProtoType = triggerMatcherProtoType,
        addGenInfo = addGenInfo
    )

    # return full sequence for production of basic tau collection,
    # generator level particle and jet matches and trigger matches
    return retVal_tau












