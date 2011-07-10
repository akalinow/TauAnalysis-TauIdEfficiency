import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.tools.configurePrePatProduction import configurePrePatProduction
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProductionTauIdEffMeasSpecific import *
 
def produceTauIdEffMeasPATTuple_base(process, isMC, isEmbedded, HLTprocessName, pfCandidateCollection, applyZrecoilCorrection):

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
            #'file:/data2/friis/CMSSW_4_2_X/skims/06-27-MatthewsZTTEvents/crab_0_110627_082505/ZTTCands_merged_v1.root'
            'file:/data1/veelken/tmp/tauIdEffSample_data_SingleMu_Run2011A_PromptReco_v4_2011Jul06_RECO_995_1_hYS.root'
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
        process.GlobalTag.globaltag = cms.string('START42_V12::All')
    else:
        process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')
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
    
    configurePrePatProduction(process, pfCandidateCollection = pfCandidateCollection, addGenInfo = isMC)
    
    #process.prePatProductionSequence.remove(process.tautagging)
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # import utility function for configurating PAT-tuple production
   
    patTupleConfig = configurePatTupleProductionTauIdEffMeasSpecific(
        process, hltProcess = HLTprocessName, isMC = isMC, isEmbedded = isEmbedded, applyZrecoilCorrection = applyZrecoilCorrection)
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    #
    # configure Jet Energy Corrections
    #
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    process.jec = cms.ESSource("PoolDBESSource",
        DBParameters = cms.PSet(
            messageLevel = cms.untracked.int32(0)
        ),
        timetype = cms.string('runnumber'),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_Jec11V2_AK5PF'),
                label  = cms.untracked.string('AK5PF')
            ),
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_Jec11V2_AK5Calo'),
                label  = cms.untracked.string('AK5Calo')
            )
        ),
        connect = cms.string('sqlite_fip:TauAnalysis/Configuration/data/Jec11V2.db')
    )
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')
    #-------------------------------------------------------------------------------------------------------------------------

    return patTupleConfig
