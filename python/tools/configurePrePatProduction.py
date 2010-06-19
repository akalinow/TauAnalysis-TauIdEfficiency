import FWCore.ParameterSet.Config as cms
import copy

def configurePrePatProduction(process, addGenInfo = False):

    #--------------------------------------------------------------------------------
    # recreate collection of PFTaus
    # with latest tags of RecoTauTag package used by TauAnalysis software
    
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    process.prePatProductionSequence = cms.Sequence(process.PFTau)
    
    ##process.load("RecoTauTag.TauTagTools.TancCVTransform_cfi")
    ##process.prePatProductionSequence += process.shrinkingConePFTauTancCVTransform
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
    # if running on Monte Carlo, produce ak5GenJets collection "on-the-fly",
    # as it is needed for matching reconstructed particles to generator level information by PAT,
    # but not contained in Monte Carlo samples produced with CMSSW_3_5_x
    if addGenInfo:
        process.load("RecoJets.Configuration.GenJetParticles_cff")
        process.load("RecoJets.JetProducers.ak5GenJets_cfi")
        process.prePatProductionSequenceGen = cms.Sequence(process.genParticlesForJets * process.ak5GenJets)
        process.prePatProductionSequence += process.prePatProductionSequenceGen
    #--------------------------------------------------------------------------------    
