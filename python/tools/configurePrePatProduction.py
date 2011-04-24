import FWCore.ParameterSet.Config as cms

def configurePrePatProduction(process, pfCandidateCollection = "particleFlow",
                              addGenInfo = False):

    #--------------------------------------------------------------------------------
    # PFCandidate pile-up removal
    # for CMSSW_4_2_0_pre8 and higher
    #process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")
    # for CMSSW_3_8_x and CMSSW_4_1_x release series
    process.load("PhysicsTools.PFCandProducer.pfNoPileUp_cff")
    process.pfPileUp.Enable = cms.bool(True)
    process.prePatProductionSequence = cms.Sequence(process.pfNoPileUpSequence)

    process.load("RecoJets.JetProducers.ak5PFJets_cfi")
    process.ak5PFJets.src = cms.InputTag(pfCandidateCollection)
    process.prePatProductionSequence += process.ak5PFJets
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # produce collections of kT PFJets for dR = 0.6 needed for rho (FastJet) pile-up corrections
    process.load("RecoJets.Configuration.RecoPFJets_cff")
    process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
    process.kt6PFJets.doRhoFastjet = True
    process.prePatProductionSequence += process.kt6PFJets
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # recreate collection of PFTaus
    # with latest tags of RecoTauTag package used by TauAnalysis software
    #
    # CV: switching from 'particleFlow' to 'pfNoPileUp' collection needs to be done in a different way
    #     for the new shrinkingCone/combinatoricTau code <-- FIXME
    #
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    process.pfRecoTauTagInfoProducer.PFCandidateProducer = cms.InputTag(pfCandidateCollection)
    tau_producers = [ 'shrinkingConePFTauProducer', 'combinatoricRecoTaus' ]
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
    
    # CV: discriminator against calo. muons currently disabled per default;
    #     add manually
    #process.prePatProductionSequence += process.hpsTancTausDiscriminationAgainstCaloMuon
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # add collection of fixed-cone PFTaus
    # (not produced per default anymore)

    process.prePatProductionSequence += process.recoTauClassicFixedConeSequence
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # add collection of shrinking-cone PFTaus
    # (not produced per default anymore)

    process.load("RecoTauTag.Configuration.ShrinkingConePFTaus_cff")

    process.prePatProductionSequence += process.produceAndDiscriminateShrinkingConePFTaus
    process.prePatProductionSequence += process.produceShrinkingConeDiscriminationByTauNeuralClassifier
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # recreate collection of CaloTaus
    # with four-momenta determined by TCTau instead of (regular) CaloTau algorithm
    #
    # NOTE:
    #      (1) change leading track Pt threshold in order to see effect
    #          of leading track finding and and leading track Pt cut
    #          on tau identification efficiency and fake-rate
    #      (2) Calo/TCTau sequence needs to be disabled when running on AOD
    #         (because production of CaloTauTagInfo objects requires TrackExtra objects)
    #
    #process.load("RecoTauTag.Configuration.RecoTCTauTag_cff")
    #
    #process.prePatProductionSequence += process.tautagging
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

    #--------------------------------------------------------------------------------
    # select collection of "good" primary event vertices with sum(trackPt) > 10 GeV,
    # used for vertex multiplicity reweighting
    process.load("TauAnalysis.RecoTools.recoVertexSelection_cff")
    process.prePatProductionSequence += process.selectedPrimaryVertexQuality
    process.prePatProductionSequence += process.selectedPrimaryVertexPosition
    process.prePatProductionSequence += process.selectedPrimaryVertexHighestPtTrackSum
    process.load("TauAnalysis.RecoTools.vertexMultiplicityReweight_cfi")
    process.prePatProductionSequence += process.selectedPrimaryVerticesTrackPtSumGt10
    #--------------------------------------------------------------------------------
