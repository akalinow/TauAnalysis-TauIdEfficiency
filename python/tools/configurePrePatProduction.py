import FWCore.ParameterSet.Config as cms

def configurePrePatProduction(process, pfCandidateCollection = "particleFlow",
                              addGenInfo = False):

    #--------------------------------------------------------------------------------
    # PFCandidate pile-up removal
    process.load("PhysicsTools.PFCandProducer.pfNoPileUp_cff")
    process.prePatProductionSequence = process.pfNoPileUpSequence

    process.load("RecoJets.JetProducers.ak5PFJets_cfi")
    process.ak5PFJets.src = cms.InputTag(pfCandidateCollection)
    process.prePatProductionSequence += process.ak5PFJets
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # recreate collection of PFTaus
    # with latest tags of RecoTauTag package used by TauAnalysis software
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    process.pfRecoTauTagInfoProducer.PFCandidateProducer = cms.InputTag(
        pfCandidateCollection)
    tau_producers = ['shrinkingConePFTauProducer', 'combinatoricRecoTaus']
    # Update the pfCandidate collection for each builder
    for producer in tau_producers:
        if hasattr(process, producer):
            my_producer = getattr(process, producer)
            if hasattr(my_producer, 'builders'):
                # Loop over the tau builders
                builders = getattr(process, producer).builders
                for builder in builders:
                    builder.pfCandSrc = cms.InputTag(pfCandidateCollection)
    process.prePatProductionSequence += process.PFTau
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # recreate collection of CaloTaus
    # with four-momenta determined by TCTau instead of (regular) CaloTau algorithm
    #
    # NOTE: change leading track Pt threshold in order to see effect
    #       of leading track finding and and leading track Pt cut
    #       on tau identification efficiency and fake-rate
    #
    process.load("RecoTauTag.Configuration.RecoTCTauTag_cff")
    ##process.JPTCaloRecoTauProducer.LeadTrack_minPt = cms.double(0.5)

    process.prePatProductionSequence += process.tautagging
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # produce collection of PFTaus reconstructed by hadron + strips (HPS) algorithm
    # "on-the-fly", as it is not contained in data taken with CMSSW_3_5_x
    # EK: no longer necessary, runs in PFTau sequence
    #process.load("RecoTauTag.Configuration.HPSPFTaus_cfi")
    #process.prePatProductionSequence += process.produceAndDiscriminateHPSPFTaus
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # produce collection of PFTaus reconstructed using ellipse for photon isolation
    process.load("TauAnalysis.TauIdEfficiency.ShrinkingConePFTausEllipticPhotonIso_cfi")
    process.prePatProductionSequence += process.produceAndDiscriminateShrinkingConePFTausEllipticPhotonIso
    ##process.prePatProductionSequence += process.produceShrinkingConePFTauEllipticPhotonIsoDiscriminationByTauNeuralClassifier
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # if running on Monte Carlo, produce ak5GenJets collection "on-the-fly",
    # as it is needed for matching reconstructed particles to generator level information by PAT,
    # but not contained in Monte Carlo samples produced with CMSSW_3_5_x
    if addGenInfo:
        process.load("RecoJets.Configuration.GenJetParticles_cff")
        process.load("RecoJets.JetProducers.ak5GenJets_cfi")
        process.prePatProductionSequenceGen = cms.Sequence(process.genParticlesForJets * process.ak5GenJets)
        process.prePatProductionSequence += process.prePatProductionSequenceGen
    #--------------------------------------------------------------------------------
