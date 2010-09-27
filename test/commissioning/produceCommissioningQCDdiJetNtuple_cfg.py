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
        ##'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/mcMinBias/muTauSkim_1_1.root'
        ##'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/mcQCDpt15/muTauSkim_1_1.root'
        ##'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/mcMinBias_pythia8/muTauSkim_1_2.root'
        ##'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/dataReReco/muTauSkim_1_1_QBQ.root'
        ##'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/mcWtauNu/muTauSkim_1_1_9yq.root'
        'file:/data1/veelken/CMSSW_3_1_2/skims/selEvents_ZtoMuTau_Ztautau_7TeV_part01.root'
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

isMC = True # use for MC (except for samples from Spring'10 reprocessing)
isSpring10 = True # use for Spring'10 reprocessed MC
##isSpring10 = False # use for non-Spring'10 reprocessed MC
##isMC = False # use for Data
#applyTrackDowngrade = False # default
applyTrackDowngrade = True # to be used for studies of systematic uncertainties only
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define GlobalTag to be used for event reconstruction
# (only relevant for HPS tau reconstruction algorithm)
if isMC:
    process.GlobalTag.globaltag = cms.string('MC_36Y_V10::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_36X_V12::All')
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
# produce collections of objects needed as input for PAT-tuple production
# (e.g. rerun reco::Tau identification algorithms with latest tags)
#
from TauAnalysis.TauIdEfficiency.tools.configurePrePatProduction import configurePrePatProduction

configurePrePatProduction(process, applyTrackDowngrade = applyTrackDowngrade, addGenInfo = isMC)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce PAT objects
#
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProduction import configurePatTupleProduction
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildQCDdiJetTauSequence

process.load("PhysicsTools.PatAlgos.cleaningLayer1.tauCleaner_cfi")

# remove jets outside kinematic range Pt > 10 GeV && |eta| < 2.5 from Tau Ntuple
# (in order to speed-up plotting macros)
patCaloTauCleanerPrototype = process.cleanPatTaus.clone(
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(),
    finalCut = cms.string(
        'caloTauTagInfoRef().jetRef().pt() > 5 & abs(caloTauTagInfoRef().jetRef().eta()) < 2.5'
    )
)

patPFTauCleanerPrototype = process.cleanPatTaus.clone(
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(),
    finalCut = cms.string(
        'pfTauTagInfoRef().pfjetRef().pt() > 5 & abs(pfTauTagInfoRef().pfjetRef().eta()) < 2.5'
    )
)

retVal = configurePatTupleProduction(
    process, patSequenceBuilder = buildQCDdiJetTauSequence,
    patPFTauCleanerPrototype = patPFTauCleanerPrototype,
    patCaloTauCleanerPrototype = patCaloTauCleanerPrototype,
    addGenInfo = isMC
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce Ntuple
#
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrigger_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigCaloTau_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauFixedCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauShrinkingCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPS_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauShrinkingConeEllipticPhotonIso_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGlobalVariables_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrackVariables_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenJets_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPhaseSpaceEventInfo_cfi")

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
        caloTaus_rec01 = process.caloTaus_recInfo.clone(
            src = cms.InputTag(retVal["caloTauCollection"])                       
        ),
        caloTaus_rec02 = process.tauTrackVariables_template.clone(
            src = cms.InputTag(retVal["caloTauCollection"])                       
        ),

        # variables specific to fixed cone PFTaus                                            
        pfTausFixedCone_rec01 = process.pfTausFixedCone_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionFixedCone"])                       
        ),
        pfTausFixedCone_rec02 = process.extraTauCandVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionFixedCone"])                       
        ),
        pfTausFixedCone_rec03 = process.tauTrackVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionFixedCone"])                       
        ),

        # variables specific to shrinking cone PFTaus                                            
        pfTausShrinkingCone_rec01 = process.pfTausShrinkingCone_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])                       
        ),
        pfTausShrinkingCone_rec02 = process.extraTauCandVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])                       
        ),                                
        pfTausShrinkingCone_rec03 = process.tauTrackVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])                       
        ),

        # variables specific to PFTaus reconstructed by hadron + strips (HPS) algorithm                                           
        pfTausHPS_rec01 = process.pfTausHPS_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionHPS"])                       
        ),
        pfTausHPS_rec02 = process.tauTrackVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionHPS"])                       
        ),

        # variables specific to shrinking cone PFTaus
        # reconstructed using ellipse for photon isolation
        pfTausShrinkingConeEllPhotonIso_rec01 = process.pfTausShrinkingConeEllipticPhotonIso_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingConeEllipticPhotonIso"])                       
        ),
        pfTausShrinkingConeEllPhotonIso_rec02 = process.extraTauCandVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingConeEllipticPhotonIso"])                       
        ),                                
        pfTausShrinkingConeEllPhotonIso_rec03 = process.tauTrackVariables_template.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingConeEllipticPhotonIso"])                       
        )
    )
)

if isMC:
    process.caloTaus_genInfo.src = cms.InputTag(retVal["caloTauCollection"])
    setattr(process.ntupleProducer.sources, "caloTaus_gen", process.caloTaus_genInfo)
    process.pfTausFixedCone_genInfo.src = cms.InputTag(retVal["pfTauCollectionFixedCone"])
    setattr(process.ntupleProducer.sources, "pfTausFixedCone_gen", process.pfTausFixedCone_genInfo)
    process.pfTausShrinkingCone_genInfo.src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])
    setattr(process.ntupleProducer.sources, "pfTausShrinkingCone_gen", process.pfTausShrinkingCone_genInfo)
    process.pfTausHPS_genInfo.src = cms.InputTag(retVal["pfTauCollectionHPS"])
    setattr(process.ntupleProducer.sources, "pfTausHPS_gen", process.pfTausHPS_genInfo)
    process.pfTausShrinkingConeEllipticPhotonIso_genInfo.src = cms.InputTag(retVal["pfTauCollectionShrinkingConeEllipticPhotonIso"])
    setattr(process.ntupleProducer.sources, "pfTausShrinkingConeEllPhotonIso_gen", process.pfTausShrinkingConeEllipticPhotonIso_genInfo)
    # add in information about generator level visible taus and all generator level jets
    setattr(process.ntupleProducer.sources, "tauGenJets", process.tauGenJets_genInfo)
    setattr(process.ntupleProducer.sources, "genJets", process.genJets_genInfo)
    setattr(process.ntupleProducer.sources, "genPhaseSpaceEventInfo", process.genPhaseSpaceEventInfo_template)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# updated InputTags for HLT trigger result object
# in case running on reprocessed Spring'10 Monte Carlo samples
if isSpring10:
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

####--------------------------------------------------------------------------------
###
### Print differences between collections of pat::Taus
### reconstructed by "regular" shrinking cone and by HPS (PF)Tau algorithm
###
### NOTE: to be used for debugging purposes only !!
###
##process.DQMStore = cms.Service("DQMStore")
##
##process.patTauCollectionDiffAnalyzer = cms.EDAnalyzer("PATTauCollectionDiffAnalyzer",
##    patTauSource1 = cms.InputTag(retVal["pfTauCollectionShrinkingCone"]),
##    patTauSource2 = cms.InputTag(retVal["pfTauCollectionHPS"]),
##                                          
##    patTauSelection = cms.string("abs(pfTauTagInfoRef.pfjetRef.eta) < 2.5 & pfTauTagInfoRef.pfjetRef.pt > 10. & userFloat('probe') > 0.5"),
##
##    dRmatch = cms.double(0.5),
##
##    dqmDirectory = cms.string('PATTauComparator')
##)                                          
##
##process.savePatTauCollectionDiffAnalyzerPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
##    outputFileName = cms.string('plotsPatTauCollectionDiffAnalyzer.root')
##)
###--------------------------------------------------------------------------------

process.p = cms.Path(
    process.prePatProductionSequence
   + process.patTupleProductionSequence
   ##+ process.patTauCollectionDiffAnalyzer
   ##+ process.savePatTauCollectionDiffAnalyzerPlots
   ##+ process.printEventContent
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

