import FWCore.ParameterSet.Config as cms

import PhysicsTools.PatAlgos.tools.helpers as patutils

def configurePrePatProduction(process, pfCandidateCollection = "particleFlow",
                              addGenInfo = False):

    process.prePatProductionSequence = cms.Sequence()
    
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

    #--------------------------------------------------------------------------------
    # PFCandidate pile-up removal
    process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")
    process.pfPileUp.Enable = cms.bool(True)
    process.pfPileUp.checkClosestZVertex = cms.bool(True)
    process.pfPileUp.Vertices = cms.InputTag('selectedPrimaryVertexPosition')
    process.prePatProductionSequence += process.pfNoPileUpSequence

    process.load("CommonTools/ParticleFlow/pfParticleSelection_cff")
    process.prePatProductionSequence += process.pfParticleSelectionSequence

    process.load("JetMETCorrections/Type1MET/pfMETCorrections_cff")
    process.prePatProductionSequence += process.producePFMETCorrections

    if not hasattr(process, "ak5PFJets"):
        process.load("RecoJets.JetProducers.ak5PFJets_cfi")
        process.ak5PFJets.src = cms.InputTag(pfCandidateCollection)
        # CV: need to enable jet area computation for 'ak5PFJets' module,
        #     in order for L1FastJet corrections to work
        #    (L1FastjetCorrector::correction function always returns 1.0 if jet area is not set)
        process.ak5PFJets.doAreaFastjet = cms.bool(True)
        process.prePatProductionSequence += process.ak5PFJets

    # smear momenta of ak5PFJets,
    # to account for Data/MC difference in PFJet resolutions (cf. JME-10-014)
    if addGenInfo:
        process.ak5PFJetsL2L3 = cms.EDProducer('PFJetCorrectionProducer',
            src = cms.InputTag('ak5PFJets'),
            correctors = cms.vstring('ak5PFL2L3')
        )
        process.prePatProductionSequence += process.ak5PFJetsL2L3

        process.load("RecoJets/Configuration/GenJetParticles_cff")
        process.load("RecoJets/Configuration/RecoGenJets_cff")
        process.prePatProductionSequence += process.genParticlesForJetsNoNu
        process.prePatProductionSequence += process.ak5GenJetsNoNu
        
        process.smearedAK5PFJets = cms.EDProducer("SmearedPFJetProducer",
            src = cms.InputTag('ak5PFJetsL2L3'),
            inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
            lutName = cms.string('pfJetResolutionMCtoDataCorrLUT'),
            smearBy = cms.double(1.0),
            srcGenJets = cms.InputTag('ak5GenJetsNoNu'),
            dRmaxGenJetMatch = cms.double(0.5)
        )
        process.prePatProductionSequence += process.smearedAK5PFJets
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # produce collections of kT PFJets for dR = 0.6 needed for rho (FastJet) pile-up corrections
    if not hasattr(process, "kt6PFJets"):
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

    # switch input to PFTau reconstruction to collection of smearedAK5PFJets,
    # in order to account for Data/MC different in PFJet resolution
    if addGenInfo:
        patutils.massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('ak5PFJets'), cms.InputTag('smearedAK5PFJets'))
    
    # CV: discriminator against calo. muons currently disabled per default;
    #     add manually
    #process.prePatProductionSequence += process.hpsTancTausDiscriminationAgainstCaloMuon
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # add collection of fixed-cone and shrinking-cone PFTaus
    # (not produced per default anymore)
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    process.prePatProductionSequence += process.recoTauClassicFixedConeSequence
    process.prePatProductionSequence += process.recoTauClassicShrinkingConeSequence
    process.prePatProductionSequence += process.recoTauClassicShrinkingConeMVASequence
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

    
