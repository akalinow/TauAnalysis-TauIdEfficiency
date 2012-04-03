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
    jetCollectionName = "patJets",
    patTauProducerPrototype = None,
    patTauCleanerPrototype = None,
    triggerMatcherProtoType = None,
    addGenInfo = False,
    applyTauJEC = False,
    applyTauVertexMatch = True):
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
    outputSequence = cms.Sequence()
    if addGenInfo:
        patTauGenParticleMatch = process.tauMatch.clone(
            src = patTauProducer.tauSource
        )
        patTauGenParticleMatchName = collectionName[0] + "GenParticleMatch" + collectionName[1]
        setattr(process, patTauGenParticleMatchName, patTauGenParticleMatch)

        patTauGenJetMatch = process.tauGenJetMatch.clone(
            src = patTauProducer.tauSource,
            checkCharge = cms.bool(False),
            maxDeltaR = cms.double(0.5),
            maxDPtRel = cms.double(999.9),
            resolveByMatchQuality = cms.bool(True)
        )
        patTauGenJetMatchName = collectionName[0] + "GenJetMatch" + collectionName[1]
        setattr(process, patTauGenJetMatchName, patTauGenJetMatch)

        patTauProducer.genParticleMatch = cms.InputTag(patTauGenParticleMatchName)
        patTauProducer.genJetMatch = cms.InputTag(patTauGenJetMatchName)
    
        outputSequence += getattr(process, patTauGenParticleMatchName)
        outputSequence += getattr(process, patTauGenJetMatchName)
    else:
        patTauProducer.addGenMatch = cms.bool(False)
        patTauProducer.addGenJetMatch = cms.bool(False)

    # configure tau-jet energy corrections
    ##if applyTauJEC:
    ##    patTauJetCorrFactors = process.patTauJetCorrFactors.clone(
    ##        src = patTauProducer.tauSource
    ##    )
    ##    patTauJetCorrFactorsName = collectionName[0] + "TauJetCorrFactors" + collectionName[1]
    ##    setattr(process, patTauJetCorrFactorsName, patTauJetCorrFactors)
    ##
    ##    print("enabling tau-JEC for %s" % patTauProducerName)
    ##    patTauProducer.tauJetCorrFactorsSource = cms.VInputTag(cms.InputTag(patTauJetCorrFactorsName))
    ##    patTauProducer.addTauJetCorrFactors = cms.bool(True)
    ##    
    ##    outputSequence += getattr(process, patTauJetCorrFactorsName)
    ##else:
    print("disabling tau-JEC for %s" % patTauProducerName)
    patTauProducer.addTauJetCorrFactors = cms.bool(False)

    # add pat::Tau producer module to sequence
    outputSequence += patTauProducer
         
    # configure matching of basic pat::Tau collection
    # to trigger primitives;
    # produce new collection of pat::Taus with trigger primitives embedded
    patTauTriggerMatch = triggerMatcherProtoType.clone(
        src = cms.InputTag(patTauProducerName),
        ##matchedCuts = cms.string(
        ##    'path("HLT_Jet30_v*")  |'
        ##    'path("HLT_Jet60_v*")  |'
        ##    'path("HLT_Jet80_v*")  |'
        ##    'path("HLT_Jet110_v*") |'
        ##    'path("HLT_Jet60_v*")'
        ##),
        matchedCuts = cms.string('type("TriggerJet")')
    )
    patTauTriggerMatchName = collectionName[0] + "TriggerMatched" + collectionName[1]
    setattr(process, patTauTriggerMatchName, patTauTriggerMatch)
    outputSequence += getattr(process, patTauTriggerMatchName)

    patTauTriggerEvent = process.patTriggerEvent.clone(
        patTriggerMatches = cms.VInputTag([patTauTriggerMatchName])
    )
    patTauTriggerEventName = collectionName[0] + "TriggerEvent" + collectionName[1]
    setattr(process, patTauTriggerEventName, patTauTriggerEvent)
    outputSequence += getattr(process, patTauTriggerEventName)

    patTauTriggerEmbedder = cms.EDProducer("PATTriggerMatchTauEmbedder",
        src     = cms.InputTag(patTauProducerName),
        matches = cms.VInputTag([patTauTriggerMatchName])
    )
    patTauTriggerEmbedderName = collectionName[0] + "TriggerEmbedder" + collectionName[1]
    setattr(process, patTauTriggerEmbedderName, patTauTriggerEmbedder)
    outputSequence += getattr(process, patTauTriggerEmbedderName)
    ##dumpPatTauTriggerEmbedded = cms.EDAnalyzer("DumpPATTaus",
    ##    src = cms.InputTag(patTauTriggerEmbedderName)
    ##)
    ##dumpPatTauTriggerEmbeddedName = "dump" + collectionName[0] + "TriggerEmbedded" + collectionName[1]
    ##setattr(process, dumpPatTauTriggerEmbeddedName, dumpPatTauTriggerEmbedded)
    ##outputSequence += getattr(process, dumpPatTauTriggerEmbeddedName)

    # configure PATTauCleaner module
    # for removal of tau-jet candidates "overlapping" with electrons or muons
    patTauCleaner = copy.deepcopy(patTauCleanerPrototype)
    patTauCleaner.src = cms.InputTag(patTauTriggerEmbedderName)
    patTauCleanerName = collectionName[0] + "Cleaned" + collectionName[1]
    setattr(process, patTauCleanerName, patTauCleaner)
    outputSequence += getattr(process, patTauCleanerName)
    ##dumpPatTauCleaned = cms.EDAnalyzer("DumpPATTaus",
    ##    src = cms.InputTag(patTauCleanerName)
    ##)
    ##dumpPatTauCleanedName = "dump" + collectionName[0] + "Cleaned" + collectionName[1]
    ##setattr(process, dumpPatTauCleanedName, dumpPatTauCleaned)
    ##outputSequence += getattr(process, dumpPatTauCleanedName)

    # require pat::Taus to originate from "the" orimary event vertex
    # (NOTE: per default, "the" orimary event vertex is the vertex of maximal sum(trackPt))
    patTauVertexMatcherName = None
    if applyTauVertexMatch:
        patTauVertexMatcher = cms.EDFilter("PATTauDzSelector",
            src = cms.InputTag(patTauCleanerName),                           
            vertexSource = cms.InputTag('selectedPrimaryVertexHighestPtTrackSum'),
            dzMax = cms.double(0.2), # [cm]                             
            filter = cms.bool(False)                                               
        )
        patTauVertexMatcherName = collectionName[0] + "VertexMatched" + collectionName[1]
        setattr(process, patTauVertexMatcherName, patTauVertexMatcher)
        outputSequence += getattr(process, patTauVertexMatcherName)
    else:
        patTauVertexMatcherName = patTauCleanerName

    # embed jetId quality flags
    patTauJetIdEmbedder = cms.EDProducer("PATTauJetIdEmbedder",
        src = cms.InputTag(patTauVertexMatcherName),
        srcJet = cms.InputTag(jetCollectionName)
    )                                         
    patTauJetIdEmbedderName = collectionName[0] + "JetIdEmbedded" + collectionName[1]
    setattr(process, patTauJetIdEmbedderName, patTauJetIdEmbedder)
    outputSequence += getattr(process, patTauJetIdEmbedderName)
    ##dumpPatTauJetIdEmbedded = cms.EDAnalyzer("DumpPATTaus",
    ##    src = cms.InputTag(patTauJetIdEmbedderName)
    ##)
    ##dumpPatTauJetIdEmbeddedName = "dump" + collectionName[0] + "JetIdEmbedded" + collectionName[1]
    ##setattr(process, dumpPatTauJetIdEmbeddedName, dumpPatTauJetIdEmbedded)
    ##outputSequence += getattr(process, dumpPatTauJetIdEmbeddedName)
    
    # return sequence for production of basic tau collection,
    # generator level particle and jet matches and trigger primitives embedded;
    # together with name of "cleaned" tau collection
    # to be used as InputTag for further processing
    retVal = {}
    retVal["sequence"] = outputSequence
    retVal["collection"] = patTauJetIdEmbedderName
    return retVal

def buildQCDdiJetTauSequence(
    process, 
    collectionName = [ "patTaus", "" ],
    jetCollectionName = "patJets",
    patTauProducerPrototype = None,
    patTauCleanerPrototype = None,
    triggerMatcherProtoType = None,
    addGenInfo = False,
    applyTauJEC = False,
    applyTauVertexMatch = True):
    '''
    blah
    '''

    # produce collection of basic pat::Taus
    # matched to generator level particles and jets
    # and trigger primitives embedded
    retVal_tau = buildTauSequence(
        process, 
        collectionName = collectionName,
        jetCollectionName = jetCollectionName,
        patTauProducerPrototype = patTauProducerPrototype,
        patTauCleanerPrototype = patTauCleanerPrototype,
        triggerMatcherProtoType = triggerMatcherProtoType,
        addGenInfo = addGenInfo,
        applyTauJEC = applyTauJEC,
        applyTauVertexMatch = applyTauVertexMatch
    )

    outputSequence = retVal_tau["sequence"]

    # produce final collection of pat::Taus with
    # flags embedded to:
    #  o indicate index of tau-jet candidate in Pt-sorted collection of jets
    #  o indicate whether tau-jet candidate is tag/probe
    patTauDijetTagAndProbe = cms.EDProducer("TauIdTagAndProbeProducer",
        source = cms.InputTag(retVal_tau["collection"]),
        triggerPaths = cms.PSet(
            HLT_Jet30  = cms.vstring('pt >  30.0'),
            HLT_Jet60  = cms.vstring('pt >  60.0'),
            HLT_Jet80  = cms.vstring('pt >  80.0'),
            HLT_Jet110 = cms.vstring('pt > 110.0'),
            HLT_Jet150 = cms.vstring('pt > 150.0')
        )
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

def buildGenericTauSequence(
    process, 
    collectionName = [ "patTaus", "" ],
    jetCollectionName = "patJets",
    patTauProducerPrototype = None,
    patTauCleanerPrototype = None,
    triggerMatcherProtoType = None,
    addGenInfo = False,
    applyTauJEC = False,
    applyTauVertexMatch = True):
    '''
    blah
    '''

    # produce collection of basic pat::Taus
    # matched to generator level particles and jets
    # and trigger primitives embedded
    retVal_tau = buildTauSequence(
        process, 
        collectionName = collectionName,
        jetCollectionName = jetCollectionName,
        patTauProducerPrototype = patTauProducerPrototype,
        patTauCleanerPrototype = patTauCleanerPrototype,
        triggerMatcherProtoType = triggerMatcherProtoType,
        addGenInfo = addGenInfo,
        applyTauJEC = applyTauJEC,
        applyTauVertexMatch = applyTauVertexMatch
    )

    # return full sequence for production of basic tau collection,
    # generator level particle and jet matches and trigger matches
    return retVal_tau











