import FWCore.ParameterSet.Config as cms

import PhysicsTools.PatAlgos.tools.helpers as patutils
import RecoMET.METProducers.METSigParams_cfi as jetResolutions

def configurePrePatProduction(process, pfCandidateCollection = "particleFlow",
                              isMC = False):
    
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

    # add reweighting factors to be applied to Monte Carlo simulated events
    # in order to match vertex multiplicity distribution in Data
    process.vertexMultiplicityReweight3dRunA = process.vertexMultiplicityReweight.clone(
        inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonMean_runs160404to173692.root"),
        type = cms.string("gen3d")
    )
    process.vertexMultiplicityReweight3dRunB = process.vertexMultiplicityReweight.clone(
        inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonMean_runs175832to180252.root"),
        type = cms.string("gen3d")
    )
    if isMC:
        process.prePatProductionSequence += process.vertexMultiplicityReweight3dRunA
        process.prePatProductionSequence += process.vertexMultiplicityReweight3dRunB
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
    ##process.dumpAK5PFJets = cms.EDAnalyzer("DumpPFJets",
    ##    src = cms.InputTag('ak5PFJets')
    ##)
    ##process.prePatProductionSequence += process.dumpAK5PFJets

    # apply jet energy corrections
    #
    # for MC   apply L1FastJet + L2 + L3 jet-energy corrections,
    # for Data apply L1FastJet + L2 + L3 + L2/L3 residual corrections
    #
    # CV: Ztautau samples produced via MCEmbedding technique are technically "Data',
    #     L2/L3 residual jet energy corrections **must not** be applied, however,
    #     since the tau-jet response is taken from the Monte Carlo simulation
    #
    process.load("JetMETCorrections/Configuration/JetCorrectionServices_cff")
    pfMEtCorrector = None
    if isMC:
        pfJetCorrector = "ak5PFL1FastL2L3"
    else:
        pfJetCorrector = "ak5PFL1FastL2L3Residual"
    process.calibratedAK5PFJets   = cms.EDProducer('PFJetCorrectionProducer',
        src = cms.InputTag('ak5PFJets'),
        correctors = cms.vstring(pfJetCorrector)
    )
    process.prePatProductionSequence += process.calibratedAK5PFJets
    ##process.dumpCalibratedAK5PFJets = cms.EDAnalyzer("DumpPFJets",
    ##    src = cms.InputTag('calibratedAK5PFJets')
    ##)
    ##process.prePatProductionSequence += process.dumpCalibratedAK5PFJets

    # smear momenta of ak5PFJets,
    # to account for Data/MC difference in PFJet resolutions (cf. JME-10-014)
    if isMC:
        process.load("RecoJets/Configuration/GenJetParticles_cff")
        process.load("RecoJets/Configuration/RecoGenJets_cff")
        process.prePatProductionSequence += process.genParticlesForJetsNoNu
        process.prePatProductionSequence += process.ak5GenJetsNoNu
        
        process.smearedAK5PFJets = cms.EDProducer("SmearedPFJetProducer",
            src = cms.InputTag('calibratedAK5PFJets'),
            dRmaxGenJetMatch = cms.string('TMath::Min(0.5, 0.1 + 0.3*TMath::Exp(-0.05*genJetPt - 10.))'),
            inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
            lutName = cms.string('pfJetResolutionMCtoDataCorrLUT'),
            jetResolutions = jetResolutions.METSignificance_params,
            skipRawJetPtThreshold = cms.double(10.), # GeV
            skipCorrJetPtThreshold = cms.double(1.e-2),
            srcGenJets = cms.InputTag('ak5GenJetsNoNu')
        )
        process.prePatProductionSequence += process.smearedAK5PFJets
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # produce collections of kT PFJets for dR = 0.6 needed for rho (FastJet) pile-up corrections
    if not hasattr(process, "kt6PFJets"):
        process.load("RecoJets.JetProducers.kt4PFJets_cfi")
        process.kt6PFJets = process.kt4PFJets.clone(
            rParam = cms.doube(0.6)
        )
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
    if isMC:
        patutils.massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('ak5PFJets'), cms.InputTag('smearedAK5PFJets'))
    else:
        patutils.massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('ak5PFJets'), cms.InputTag('calibratedAK5PFJets'))
    
    # CV: discriminator against calo. muons currently disabled per default;
    #     add manually
    #process.prePatProductionSequence += process.hpsTancTausDiscriminationAgainstCaloMuon
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # add collection of fixed-cone and shrinking-cone PFTaus
    # (not produced per default anymore)
    #process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    #process.prePatProductionSequence += process.recoTauClassicFixedConeSequence
    #process.prePatProductionSequence += process.recoTauClassicShrinkingConeSequence
    #process.prePatProductionSequence += process.recoTauClassicShrinkingConeMVASequence
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
    if isMC:
        process.load("RecoJets.Configuration.GenJetParticles_cff")
        process.load("RecoJets.JetProducers.ak5GenJets_cfi")
        process.prePatProductionSequenceGen = cms.Sequence(process.genParticlesForJets * process.ak5GenJets)
        process.prePatProductionSequence += process.prePatProductionSequenceGen
    #--------------------------------------------------------------------------------

    
