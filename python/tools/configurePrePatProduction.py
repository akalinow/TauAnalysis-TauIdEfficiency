import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.tools.configureParticleFlowInput import setParticleFlowTrackInput, _setInputTag

def configurePrePatProduction(process, applyTrackDowngrade = False, addGenInfo = False):

    #--------------------------------------------------------------------------------
    # rerun particle flow using as input collection of reco::Tracks
    # for which certain fraction of entries has been killed
    # (in order to estimate effect of track reconstruction inefficiency
    #  on PFChargedHadron and PFGammas based isolation of tau jet candidates)
    if applyTrackDowngrade:
        process.load("RecoLocalTracker/SiPixelRecHits/SiPixelRecHits_cfi")
        process.load("RecoLocalTracker/SiStripRecHitConverter/SiStripRecHitConverter_cfi")
        process.prePatProductionSequence = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)

        process.load("RecoParticleFlow/PFClusterProducer/particleFlowCluster_cff")
        process.prePatProductionSequence += process.particleFlowCluster
        
        process.load("RecoTracker/Configuration/RecoTracker_cff")
        process.prePatProductionSequence += process.ckftracks

        process.downgradedGeneralTracks = cms.EDAnalyzer("TrackAndTrajectoryDowngrade",
            src = cms.InputTag('generalTracks'),
            pDowngrade = cms.double(0.05)
        )                                                     

        process.prePatProductionSequence += process.downgradedGeneralTracks

        process.load("TrackingTools/GsfTracking/GsfElectronTracking_cff")
        process.prePatProductionSequence += process.electronGsfTracking
    
        process.load("RecoParticleFlow/Configuration/RecoParticleFlow_cff")
        setParticleFlowTrackInput(process, 'downgradedGeneralTracks')
        
        process.prePatProductionSequence += process.trackerDrivenElectronSeeds
        process.prePatProductionSequence += process.trackerOnlyConversions
        process.prePatProductionSequence += process.generalV0Candidates
        process.prePatProductionSequence += process.particleFlowReco
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # recreate collection of PFTaus
    # with latest tags of RecoTauTag package used by TauAnalysis software
    
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    if hasattr(process, "prePatProductionSequence"):
        process.prePatProductionSequence += process.PFTau
    else:
        process.prePatProductionSequence = cms.Sequence(process.PFTau)
    
    ##process.load("RecoTauTag.TauTagTools.TancCVTransform_cfi")
    ##process.prePatProductionSequence += process.shrinkingConePFTauTancCVTransform

    if applyTrackDowngrade:
        _setInputTag(process, "ak5PFJetTracksAssociatorAtVertex", "tracks", cms.InputTag('downgradedGeneralTracks'))
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

    if applyTrackDowngrade:
        _setInputTag(process, "trackExtrapolator", "trackSrc", cms.InputTag('downgradedGeneralTracks'))
        process.prePatProductionSequence += process.trackExtrapolator
        _setInputTag(process, "ak5JetTracksAssociatorAtVertex", "tracks", cms.InputTag('downgradedGeneralTracks'))
        _setInputTag(process, "ak5JetTracksAssociatorAtCaloFace", "tracks", cms.InputTag('downgradedGeneralTracks'))
        process.prePatProductionSequence += process.ak5JTA
        _setInputTag(process, "JPTAntiKt5JetTracksAssociatorAtVertex", "tracks", cms.InputTag('downgradedGeneralTracks'))
        _setInputTag(process, "caloRecoTauProducer", "TrackCollection", cms.InputTag('downgradedGeneralTracks'))
    
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
