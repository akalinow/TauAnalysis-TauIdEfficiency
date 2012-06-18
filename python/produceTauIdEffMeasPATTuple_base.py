import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.tools.configurePrePatProduction import configurePrePatProduction
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProductionTauIdEffMeasSpecific import *
 
def produceTauIdEffMeasPATTuple_base(process, isMC, isEmbedded, HLTprocessName, pfCandidateCollection, runSVfit):

    # import of standard configurations for RECOnstruction
    # of electrons, muons and tau-jets with non-standard isolation cones
    process.load('Configuration/StandardSequences/Services_cff')
    process.load('FWCore/MessageService/MessageLogger_cfi')
    process.MessageLogger.cerr.FwkReport.reportEvery = 100
    #process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
    #process.MessageLogger.suppressInfo = cms.untracked.vstring()
    process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
    process.load('Configuration/StandardSequences/GeometryIdeal_cff')
    process.load('Configuration/StandardSequences/MagneticField_cff')
    process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

    #--------------------------------------------------------------------------------
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            ##'file:/data1/veelken/CMSSW_5_2_x/skims/simZplusJets_AOD_1_1_ZkM.root'
            ##'file:/data1/veelken/CMSSW_5_2_x/skims/tauIdEffSample_TTplusJets_madgraph2_2012May12_AOD_97_1_8KH.root'
            'file:/data1/veelken/CMSSW_5_2_x/skims/data2012runA_doubleMu_AOD_1_1_Fzg.root'
        )
    )

    # print event content 
    process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

    # print gen. particle information
    process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
        src = cms.InputTag("genParticles"),
        maxEventsToPrint = cms.untracked.int32(100)
    )

    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
    )
    #---------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # define GlobalTag to be used for event reconstruction
    # (only relevant for HPS tau reconstruction algorithm)
    if isMC:
        process.GlobalTag.globaltag = cms.string('START52_V9::All')
    else:
        process.GlobalTag.globaltag = cms.string('GR_R_52_V7::All')
    #--------------------------------------------------------------------------------    

    #--------------------------------------------------------------------------------
    # define skimming criteria
    # (in order to be able to produce Tau Ntuple directly from unskimmed Monte Carlo/datasets;
    #  HLT single jet trigger passed && either two CaloJets or two PFJets of Pt > 10 GeV within |eta| < 2.5)
    process.load('TauAnalysis.TauIdEfficiency.filterTauIdEffSample_cfi')

    process.hltMu.selector.src = cms.InputTag('TriggerResults::%s' % HLTprocessName)

    if isMC or isEmbedded:
        process.dataQualityFilters.remove(process.hltPhysicsDeclared)
        process.dataQualityFilters.remove(process.dcsstatus)
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # produce collections of objects needed as input for PAT-tuple production
    # (e.g. rerun reco::Tau identification algorithms with latest tags)
    #
    
    configurePrePatProduction(process, pfCandidateCollection = pfCandidateCollection, isMC = isMC)
    
    #process.prePatProductionSequence.remove(process.tautagging)
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # import utility function for configurating PAT-tuple production
   
    patTupleConfig = configurePatTupleProductionTauIdEffMeasSpecific(
        process, hltProcess = HLTprocessName, isMC = isMC, isEmbedded = isEmbedded, runSVfit = runSVfit)
    #--------------------------------------------------------------------------------

    return patTupleConfig
