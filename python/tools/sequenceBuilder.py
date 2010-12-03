import FWCore.ParameterSet.Config as cms
import copy

'''

tools.sequenceBuilder.py

Create pat::Taus from a specified collection, match
them to the tag and probe jets, and create ntuples
of the resulting collections.

Author: Evan K. Friis, Christian Veelken (UC Davis)

'''

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
    #
    # NOTE: for matching generated to reconstructed tau-jets
    #       set matching cone size to dR = 0.5
    #       and disable requirement on transverse momentum
    #       of generated and reconstructed tau-jets 
    #
    outputSequence = None
    if addGenInfo:
        patTauGenParticleMatch = process.tauMatch.clone(
            src = patTauProducer.tauSource
        )
        patTauGenParticleMatchName = collectionName[0] + "GenParticleMatch" + collectionName[1]
        setattr(process, patTauGenParticleMatchName, patTauGenParticleMatch)

        patTauGenJetMatch = process.tauGenJetMatch.clone(
            src = patTauProducer.tauSource,
            maxDeltaR = cms.double(0.5),
            maxDPtRel = cms.double(999.9)
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
    outputSequence += getattr(process, patTauTriggerMatchName)

    patTauTriggerEvent = process.patTriggerEvent.clone(
        patTriggerMatches = cms.VInputTag(patTauTriggerMatchName)
    )
    patTauTriggerEventName = collectionName[0] + "TriggerEvent" + collectionName[1]
    setattr(process, patTauTriggerEventName, patTauTriggerEvent)
    outputSequence += getattr(process, patTauTriggerEventName)

    patTauTriggerEmbedder = cms.EDProducer("PATTriggerMatchTauEmbedder",
        src     = cms.InputTag(patTauProducerName),
        matches = cms.VInputTag(patTauTriggerMatchName)
    )
    patTauTriggerEmbedderName = collectionName[0] + "TriggerEmbedder" + collectionName[1]
    setattr(process, patTauTriggerEmbedderName, patTauTriggerEmbedder)
    outputSequence += getattr(process, patTauTriggerEmbedderName)

    # configure PATTauCleaner module
    # for removal of tau-jet candidates "overlapping" with electrons or muons
    patTauCleaner = copy.deepcopy(patTauCleanerPrototype)
    patTauCleaner.src = cms.InputTag(patTauTriggerEmbedderName)
    patTauCleanerName = collectionName[0] + "Cleaned" + collectionName[1]
    setattr(process, patTauCleanerName, patTauCleaner)
    outputSequence += getattr(process, patTauCleanerName)

    # embed loose PFIsolation sums into pat::Tau objects
    patTauLoosePFIsoEmbedder = cms.EDProducer("PATTauPFIsolationEmbedder",
        src = cms.InputTag(patTauCleanerName),                                       
        userFloatName = cms.string('pfLooseIsoPt'),
        pfCandidateSource = cms.InputTag('pfNoPileUp'),
        chargedHadronIso = cms.PSet(
            ptMin = cms.double(1.0),        
            dRvetoCone = cms.double(0.15),
            dRisoCone = cms.double(0.6)
        ),
        neutralHadronIso = cms.PSet(
            ptMin = cms.double(1000.),        
            dRvetoCone = cms.double(0.15),        
            dRisoCone = cms.double(0.6)
        ),
        photonIso = cms.PSet(
            ptMin = cms.double(1.5),        
            dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
            dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
            dRvetoCone = cms.double(0.15),
            dRisoCone = cms.double(0.6)
        )
    )
    patTauLoosePFIsoEmbedderName = collectionName[0] + "LoosePFIsoEmbedded" + collectionName[1]
    setattr(process, patTauLoosePFIsoEmbedderName, patTauLoosePFIsoEmbedder)
    outputSequence += getattr(process, patTauLoosePFIsoEmbedderName)

    # return sequence for production of basic tau collection,
    # generator level particle and jet matches and trigger primitives embedded;
    # together with name of "cleaned" tau collection
    # to be used as InputTag for further processing
    retVal = {}
    retVal["sequence"] = outputSequence
    retVal["collection"] = patTauLoosePFIsoEmbedderName
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
    outputSequence += getattr(process, patTauDijetTagAndProbeName)

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

def buildWplusJetsEnrichedTauSequence(
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

def buildZmumuEnrichedTauSequence(
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












