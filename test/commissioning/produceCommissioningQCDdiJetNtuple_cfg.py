import FWCore.ParameterSet.Config as cms

process = cms.Process("prodCommissioningQDCdiJetNtuple")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#--------------------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0021/F405BC9A-525D-DF11-AB96-002618943811.root',
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0020/EE3E8F74-365D-DF11-AE3D-002618FDA211.root'
        ##'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/data/muTauSkim_1_1.root'
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/mcMinBias/muTauSkim_1_1.root'           
    ),
    skipEvents = cms.untracked.uint32(0)            
)

# To prevent file name CVS battles         
#import os
#if not os.path.exists('/usr/bin/nsls'):
#    process.source.fileNames = cms.untracked.vstring(
#        '/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0021/F405BC9A-525D-DF11-AB96-002618943811.root',
#        '/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0020/EE3E8F74-365D-DF11-AE3D-002618FDA211.root'
#    )

# print event content 
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

isMC = True # use for MC
##isMC = False # use for Data
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define GlobalTag to be used for event reconstruction
# (only relevant for HPS tau reconstruction algorithm)
if isMC:
    process.GlobalTag.globaltag = cms.string('MC_36Y_V7A::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_36X_V11A::All')
#--------------------------------------------------------------------------------    

#--------------------------------------------------------------------------------
# define skimming criteria
# (in order to be able to produce Tau Ntuple directly from unskimmed Monte Carlo/datasets;
#  HLT single jet trigger passed && either two CaloJets or two PFJets of Pt > 10 GeV within |eta| < 2.5)
process.load('TauAnalysis.TauIdEfficiency.filterQCDdiJet_cfi')
if isMC:
    process.dataQualityFilters.remove(process.hltPhysicsDeclared)
    process.dataQualityFilters.remove(process.dcsstatus)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce PAT objects
#
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProduction import configurePatTupleProduction
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildQCDdiJetTauSequence

process.load("PhysicsTools.PatAlgos.cleaningLayer1.tauCleaner_cfi")
process.patTauCleanerPrototype = process.cleanPatTaus.clone(
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(),
    finalCut = cms.string('')
)

retVal = configurePatTupleProduction(process,
                                     patSequenceBuilder = buildQCDdiJetTauSequence,
                                     patTauCleanerPrototype = process.patTauCleanerPrototype,
                                     addGenInfo = isMC)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce Ntuple
#
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrigger_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGlobalVariables_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigCaloTau_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauFixedCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauShrinkingCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPS_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenJets_cfi")

process.ntupleProducer = cms.EDAnalyzer("ObjValEDNtupleProducer",
                                        
    ntupleName = cms.string("tauIdEffNtuple"),
                                        
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc

        # variables indicating decision of HLT trigger paths
        trigger = process.trigger_template,                                    

        # variables specifying x,y,z coordinates of primary event vertex
        vertex = process.vertex_template,                   

        # global variables describing the underlying event/
        # amount of hadronic activity                                            
        jets = process.jets_template,
        met = process.met_template,                                    

        # variables specific to CaloTaus                                            
        caloTaus_part01 = process.caloTaus_recInfo.clone(
            src = cms.InputTag(retVal["caloTauCollection"])                       
        ),

        # variables specific to fixed cone PFTaus                                            
        pfTausFixedCone_part01 = process.pfTausFixedCone_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionFixedCone"])                       
        ),
        pfTausFixedCone_part02 = process.extraTauCandVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionFixedCone"])                       
        ),                                   

        # variables specific to shrinking cone PFTaus                                            
        pfTausShrinkingCone_part01 = process.pfTausShrinkingCone_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])                       
        ),
        pfTausShrinkingCone_part02 = process.extraTauCandVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])                       
        ),                                

        # variables specific to PFTaus reconstructed by hadron + strips (HPS) algorithm                                           
        pfTausHPS_part01 = process.pfTausHPS_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionHPS"])                       
        )
    )
)

if isMC:
    process.caloTaus_genInfo.src = cms.InputTag(retVal["caloTauCollection"])
    setattr(process.ntupleProducer.sources, "caloTaus_part02", process.caloTaus_genInfo)
    process.pfTausFixedCone_genInfo.src = cms.InputTag(retVal["pfTauCollectionFixedCone"])
    setattr(process.ntupleProducer.sources, "pfTausFixedCone_part03", process.pfTausFixedCone_genInfo)
    process.pfTausShrinkingCone_genInfo.src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])
    setattr(process.ntupleProducer.sources, "pfTausShrinkingCone_part03", process.pfTausShrinkingCone_genInfo)
    process.pfTausHPS_genInfo.src = cms.InputTag(retVal["pfTauCollectionHPS"])
    setattr(process.ntupleProducer.sources, "pfTausHPS_part02", process.pfTausHPS_genInfo)
    # add in information about generator level visible taus and all generator level jets
    setattr(process.ntupleProducer.sources, "tauGenJets", process.tauGenJets_genInfo)
    setattr(process.ntupleProducer.sources, "genJets", process.genJets_genInfo)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# updated InputTags for HLT trigger result object
# in case running on reprocessed Spring'10 Monte Carlo samples
if isMC:
    process.hltJet15U.selector.src = cms.InputTag('TriggerResults::REDIGI')
    process.patTrigger.processName = cms.string('REDIGI')
    process.patCaloTausTriggerEvent.processName = cms.string('REDIGI')
    process.patPFTausTriggerEventFixedCone.processName = cms.string('REDIGI')
    process.patPFTausTriggerEventShrinkingCone.processName = cms.string('REDIGI')
    process.patPFTausTriggerEventHPS.processName = cms.string('REDIGI')    
    process.ntupleProducer.sources.trigger.src = cms.InputTag('TriggerResults::REDIGI')
#--------------------------------------------------------------------------------    

#--------------------------------------------------------------------------------
#
# Save ntuple
#
process.ntupleOutputModule = cms.OutputModule("PoolOutputModule",
    cms.PSet(
        outputCommands = cms.untracked.vstring(
            "drop *",
            "keep *_*ntupleProducer*_*_*"
        )                               
    ),
    process.qcdDiJetEventSelection, # comment-out to disable filtering of events put in Tau Ntuple                         
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string("tauIdEffEDNtuple_qcdDiJet.root")      
)
#--------------------------------------------------------------------------------

# recreate collection of CaloTaus
# with four-momenta determined by TCTau instead of (regular) CaloTau algorithm
process.load("RecoTauTag.Configuration.RecoTCTauTag_cff")
process.prePatProductionSequence = cms.Sequence(process.tautagging)

# produce collection of PFTaus reconstructed by hadron + strips (HPS) algorithm
# "on-the-fly", as it is not contained in data taken with CMSSW_3_5_x
process.load("RecoTauTag.Configuration.HPSPFTaus_cfi")
process.prePatProductionSequence += process.produceAndDiscriminateHPSPFTaus

# if running on Monte Carlo, produce ak5GenJets collection "on-the-fly",
# as it is needed for matching reconstructed particles to generator level information by PAT,
# but not contained in Monte Carlo samples produced with CMSSW_3_5_x
if isMC:
    process.load("RecoJets.Configuration.GenJetParticles_cff")
    process.load("RecoJets.JetProducers.ak5GenJets_cfi")
    process.prePatProductionSequenceGen = cms.Sequence(process.genParticlesForJets * process.ak5GenJets)
    process.prePatProductionSequence += process.prePatProductionSequenceGen

process.p = cms.Path(
    process.prePatProductionSequence
   + process.patTupleProductionSequence
   #+ process.printEventContent
   + process.ntupleProducer
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath(process.ntupleOutputModule)

# define order in which different paths are run
process.schedule = cms.Schedule(
   process.p,
   process.caloTauSkimPath,
   process.pfTauSkimPath,
   process.o
)

# print-out all python configuration parameter information
#print process.dumpPython()

