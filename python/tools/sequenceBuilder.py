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
    addGenInfo = False,
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
    outputSequence = None
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
        src = cms.InputTag(patTauProducerName),
        filterLabels = cms.vstring(
            'hltSingleJet30',
            'hltSingleJet60Regional',
            'hltSingleJet80Regional',
            'hltSingleJet110Regional',
            'hltSingleJet150Regional'
        ),
        pathNames    = cms.vstring(
            'HLT_Jet30_v1',
            'HLT_Jet60_v1',
            'HLT_Jet80_v1',
            'HLT_Jet110_v1',
            'HLT_Jet150_v1'
        )
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

    # embed loose PFIsolation sums into pat::Tau objects
    patTauLoosePFIsoEmbedder04 = cms.EDProducer("PATTauPFIsolationEmbedder",
        src = cms.InputTag(patTauVertexMatcherName),                                       
        userFloatName = cms.string('pfLooseIsoPt04'),
        pfCandidateSource = cms.InputTag('pfNoPileUp'),
        chargedHadronIso = cms.PSet(
            ptMin = cms.double(1.0),        
            dRvetoCone = cms.double(0.15),
            dRisoCone = cms.double(0.4)
        ),
        neutralHadronIso = cms.PSet(
            ptMin = cms.double(1000.),        
            dRvetoCone = cms.double(0.15),        
            dRisoCone = cms.double(0.)
        ),
        photonIso = cms.PSet(
            ptMin = cms.double(1.5),        
            dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
            dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
            dRvetoCone = cms.double(0.15),
            dRisoCone = cms.double(0.4)
        ),
        direction = cms.string('track')
    )
    patTauLoosePFIsoEmbedder04Name = collectionName[0] + "LoosePFIsoEmbedded04" + collectionName[1]
    setattr(process, patTauLoosePFIsoEmbedder04Name, patTauLoosePFIsoEmbedder04)
    outputSequence += getattr(process, patTauLoosePFIsoEmbedder04Name)

    patTauLoosePFIsoEmbedder06 = patTauLoosePFIsoEmbedder04.clone(
        src = cms.InputTag(patTauLoosePFIsoEmbedder04Name),                                       
        userFloatName = cms.string('pfLooseIsoPt06'),
        chargedHadronIso = patTauLoosePFIsoEmbedder04.chargedHadronIso.clone(
            dRisoCone = cms.double(0.6)
        ),
        neutralHadronIso = patTauLoosePFIsoEmbedder04.neutralHadronIso.clone(
            dRisoCone = cms.double(0.6)
        ),
        photonIso = patTauLoosePFIsoEmbedder04.photonIso.clone(
            dRisoCone = cms.double(0.6)
        )
    )
    patTauLoosePFIsoEmbedder06Name = collectionName[0] + "LoosePFIsoEmbedded06" + collectionName[1]
    setattr(process, patTauLoosePFIsoEmbedder06Name, patTauLoosePFIsoEmbedder06)
    outputSequence += getattr(process, patTauLoosePFIsoEmbedder06Name)

    # return sequence for production of basic tau collection,
    # generator level particle and jet matches and trigger primitives embedded;
    # together with name of "cleaned" tau collection
    # to be used as InputTag for further processing
    retVal = {}
    retVal["sequence"] = outputSequence
    retVal["collection"] = patTauLoosePFIsoEmbedder06Name
    return retVal

def buildQCDdiJetTauSequence(
    process, 
    collectionName = [ "patTaus", "" ],
    patTauProducerPrototype = None,
    patTauCleanerPrototype = None,
    triggerMatcherProtoType = None,
    addGenInfo = False,
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
        triggerPaths = cms.vstring('HLT_Jet15U', 'HLT_Jet30U', 'HLT_Jet50U', 'HLT_Jet70U', 'HLT_Jet100U')
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
    patTauProducerPrototype = None,
    patTauCleanerPrototype = None,
    triggerMatcherProtoType = None,
    addGenInfo = False,
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
        patTauProducerPrototype = patTauProducerPrototype,
        patTauCleanerPrototype = patTauCleanerPrototype,
        triggerMatcherProtoType = triggerMatcherProtoType,
        addGenInfo = addGenInfo
    )

    # return full sequence for production of basic tau collection,
    # generator level particle and jet matches and trigger matches
    return retVal_tau











