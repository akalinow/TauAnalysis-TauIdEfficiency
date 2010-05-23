import FWCore.ParameterSet.Config as cms

process = cms.Process("prodCommissioningDijetNtuple")

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
process.GlobalTag.globaltag = cms.string('MC_36Y_V7A::All')

#--------------------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0021/F405BC9A-525D-DF11-AB96-002618943811.root',
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0020/EE3E8F74-365D-DF11-AE3D-002618FDA211.root'
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_1_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_2_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_3_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_4_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_5_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_6_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_7_2.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_8_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_9_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_10_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_11_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_12_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_13_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_14_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_15_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_16_2.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_17_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_18_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_19_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_20_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_21_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_22_2.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_23_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_24_2.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_25_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_26_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_27_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_28_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_29_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_30_1.root',
        'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/skims/tauCommissioning/run135528/qcdDiJetSkim_31_1.root'
    ),
    skipEvents = cms.untracked.uint32(0)            
)

# print event content 
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

addGenInfo = True # use for MC
##addGenInfo = False # use for Data
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce PAT objects
#
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProduction import configurePatTupleProduction

configurePatTupleProduction(process, addGenInfo = addGenInfo)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce Ntuple
#
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigCaloTau_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauFixedCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauShrinkingCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPS_cfi")

process.ntupleProducer = cms.EDAnalyzer("ObjValEDNtupleProducer",
                                        
    ntupleName = cms.string("tauIdEffNtuple"),
                                        
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc

        # variables specifying x,y,z coordinates of primary event vertex
        vertex = process.vertex_template,                   

        # variables specific to CaloTaus                                            
        caloTaus_part01 = process.caloTaus_recInfo,
        ##caloTaus_part02 = process.caloTaus_genInfo,                                    

        # variables specific to fixed cone PFTaus                                            
        pfTausFixedCone_part01 = process.pfTausFixedCone_recInfo,
        ##pfTausFixedCone_part02 = process.pfTausFixedCone_genInfo,

        # variables specific to shrinking cone PFTaus                                            
        pfTausShrinkingCone_part01 = process.pfTausShrinkingCone_recInfo,
        ##pfTausShrinkingCone_part02 = process.pfTausShrinkingCone_genInfo,

        # variables specific to PFTaus reconstructed by hadron + strips (HPS) algorithm                                           
        pfTausHPS_part01 = process.pfTausHPS_recInfo,
        ##pfTausHPS_part02 = process.pfTausHPS_genInfo                                    
    )
)

if addGenInfo:
    setattr(process.ntupleProducer.sources, "caloTaus_part02", process.caloTaus_genInfo)
    setattr(process.ntupleProducer.sources, "pfTausFixedCone_part02", process.pfTausFixedCone_genInfo)
    setattr(process.ntupleProducer.sources, "pfTausShrinkingCone_part02", process.pfTausShrinkingCone_genInfo)
    setattr(process.ntupleProducer.sources, "pfTausHPS_part02", process.pfTausHPS_genInfo)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# Save ntuple
#
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*ntupleProducer*_*_*"
    ),
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string("tauIdEff_ntuple.root")      
)
#--------------------------------------------------------------------------------

process.p = cms.Path(
    process.patTupleProductionSequence
   #+ process.printEventContent
   + process.ntupleProducer
)

process.end = cms.EndPath(process.out)

# print-out all python configuration parameter information
#print process.dumpPython()
