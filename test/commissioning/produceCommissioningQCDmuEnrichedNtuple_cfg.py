import FWCore.ParameterSet.Config as cms

process = cms.Process("prodCommissioningQDCmuEnrichedNtuple")

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
        '/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0021/F405BC9A-525D-DF11-AB96-002618943811.root',
        '/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0020/EE3E8F74-365D-DF11-AE3D-002618FDA211.root'
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_1_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_2_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_4_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_5_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_6_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_7_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_8_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_9_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_10_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_11_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_12_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_13_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_14_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_15_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_16_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_17_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_18_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_19_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_20_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_21_1.root',
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528sdJetMETTau_noTriggerSel/qcdDiJetSkim_22_1.root'
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
#  HLT muon trigger passed && global muon of Pt > 3 GeV within |eta| < 2.5
#                          && either CaloJet or PFJet of Pt > 10 GeV within |eta| < 2.5)
process.load('TauAnalysis.TauIdEfficiency.filterQCDmuEnriched_cfi')
if isMC:
    process.dataQualityFilters.remove(process.hltPhysicsDeclared)
    process.dataQualityFilters.remove(process.dcsstatus)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce PAT objects
#
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProduction import configurePatTupleProduction
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildQCDmuEnrichedTauSequence

# define muon selection criteria
# (used for "cleaning" of CaloTau/PFTau collection)
#
# NOTE: only muons passing selection will be written to Tau Ntuple
#
##process.load("PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi")
##process.selectedPatMuons.cut = cms.string("isGlobalMuon & pt > 3 & abs(eta) < 2.5")

# "clean" CaloTau/PFTau collections
# (i.e. remove CaloTaus/PFTaus overlapping with muons)
process.load("PhysicsTools.PatAlgos.cleaningLayer1.tauCleaner_cfi")
process.patTauCleanerPrototype = process.cleanPatTaus.clone(
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src                 = cms.InputTag("selectedPatMuons"),
           algorithm           = cms.string("byDeltaR"),
           preselection        = cms.string("isGlobalMuon"),
           deltaR              = cms.double(0.7),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True)
        )
    ),        
    finalCut = cms.string('')
)

retVal = configurePatTupleProduction(process,
                                     patSequenceBuilder = buildQCDmuEnrichedTauSequence,
                                     patTauCleanerPrototype = process.patTauCleanerPrototype,
                                     addGenInfo = isMC)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce Ntuple
#
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrigger_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigMuon_cfi")
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

        # variables specific to Muons
        muons_part01 = process.muons_recInfo,              

        # variables specific to CaloTaus                                            
        caloTaus_part01 = process.caloTaus_recInfo.clone(
            src = cms.InputTag(retVal["caloTauCollection"])                       
        ),

        # variables specific to fixed cone PFTaus                                            
        pfTausFixedCone_part01 = process.pfTausFixedCone_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionFixedCone"])                       
        ),

        # variables specific to shrinking cone PFTaus                                            
        pfTausShrinkingCone_part01 = process.pfTausShrinkingCone_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])                       
        ),

        # variables specific to PFTaus reconstructed by hadron + strips (HPS) algorithm                                           
        pfTausHPS_part01 = process.pfTausHPS_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionHPS"])                       
        )
    )
)

if isMC:
    setattr(process.ntupleProducer.sources, "muons_part02", process.muons_genInfo)
    process.caloTaus_genInfo.src = cms.InputTag(retVal["caloTauCollection"])
    setattr(process.ntupleProducer.sources, "caloTaus_part02", process.caloTaus_genInfo)
    process.pfTausFixedCone_genInfo.src = cms.InputTag(retVal["pfTauCollectionFixedCone"])
    setattr(process.ntupleProducer.sources, "pfTausFixedCone_part02", process.pfTausFixedCone_genInfo)
    process.pfTausShrinkingCone_genInfo.src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])
    setattr(process.ntupleProducer.sources, "pfTausShrinkingCone_part02", process.pfTausShrinkingCone_genInfo)
    process.pfTausHPS_genInfo.src = cms.InputTag(retVal["pfTauCollectionHPS"])
    setattr(process.ntupleProducer.sources, "pfTausHPS_part02", process.pfTausHPS_genInfo)
    # Add in information about generator level visible taus and all generator level jets
    setattr(process.ntupleProducer.sources, "tauGenJets", process.tauGenJets_genInfo)
    setattr(process.ntupleProducer.sources, "genJets", process.genJets_genInfo)
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
    process.qcdMuEnrichedEventSelection, # comment-out to disable filtering of events put in Tau Ntuple                         
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string("tauIdEffEDNtuple_qcdMuEnriched.root")      
)
#--------------------------------------------------------------------------------

# produce collections of PFTaus reconstructed by hadron + strips (HPS) algorithm
# "on-the-fly", as it is not contained in data taken with CMSSW_3_5_x
process.load("RecoTauTag.Configuration.HPSPFTaus_cfi")
process.prePatProductionSequence = cms.Sequence(process.produceAndDiscriminateHPSPFTaus)

# if running on Monte Carlo, produce ak5GenJets collection "on-the-fly",
# as it is needed for matching reconstructed particles to generator level information by PAT,
# but not contained in Monte Carlo samples produced with CMSSW_3_5_x
if isMC:
    process.load("RecoJets.Configuration.GenJetParticles_cff")
    process.load("RecoJets.JetProducers.ak5GenJets_cfi")
    process.prePatProductionSequence += cms.Sequence(process.genParticlesForJets + process.ak5GenJets)

##process.muTauEventDump = cms.EDAnalyzer("EventDumpAnalyzer",
##    plugin = cms.PSet(
##        pluginType = cms.string('MuTauEventDump'),                                    
##                                            
##        l1GtReadoutRecordSource = cms.InputTag(''),
##        l1GtObjectMapRecordSource = cms.InputTag(''),
##        l1BitsToPrint = cms.vstring('L1_SingleMu3', 'L1_SingleMu5'),
##    
##        hltResultsSource = cms.InputTag('TriggerResults::HLT'),
##        hltPathsToPrint = cms.vstring('HLT_Mu3', 'HLT_Mu5'),
##        
##        genParticleSource = cms.InputTag('genParticles'),
##        genJetSource = cms.InputTag('iterativeCone5GenJets'),
##        genTauJetSource = cms.InputTag('tauGenJets'),
##        genEventInfoSource = cms.InputTag('generator'),
##    
##        electronSource = cms.InputTag('cleanPatElectrons'),
##        muonSource = cms.InputTag('cleanPatMuons'),
##        tauSource = cms.InputTag('patPFTausCleanedShrinkingCone'),
##        printTauIdEfficiencies = cms.bool(False),
##        diTauCandidateSource = cms.InputTag('patMuonPFTauPairsShrinkingCone'),
##        muTauZmumuHypothesisSource = cms.InputTag(''),
##        diMuZmumuHypothesisSource = cms.InputTag(''),
##        jetSource = cms.InputTag('patJets'),
##        caloMEtSource = cms.InputTag('patMETs'),
##        pfMEtSource = cms.InputTag('patPFMETs'),
##        genMEtSource = cms.InputTag('genMetTrue'),
##    
##        output = cms.string("std::cout")
##    )
##)

process.p = cms.Path(
    process.prePatProductionSequence
   + process.patTupleProductionSequence
   #+ process.printEventContent
   ##+ process.muTauEventDump
   + process.ntupleProducer
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath(process.ntupleOutputModule)

# print-out all python configuration parameter information
#
# NOTE: need to delete empty sequence produced by call to "switchJetCollection"
#       in order to avoid error when calling "process.dumpPython"
#      ( cf. https://hypernews.cern.ch/HyperNews/CMS/get/physTools/1688/1/1/1/1/1.html )
#
#del process.patJetMETCorrections
#print process.dumpPython()

