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

#--------------------------------------------------------------------------------

readFiles = cms.untracked.vstring()

readFiles.extend( [
       '/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0021/F405BC9A-525D-DF11-AB96-002618943811.root',
       '/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0020/EE3E8F74-365D-DF11-AE3D-002618FDA211.root'    
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0054/90776DFB-E149-DF11-9BCC-0017A4771004.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/FAC0F484-D749-DF11-B689-001E0B5FA500.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/F677A203-D549-DF11-A125-00237DA096F8.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/E419AC73-D549-DF11-8778-0017A4771008.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/C41AC50B-D549-DF11-8C61-0017A4771010.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/B8933E15-D649-DF11-8CD5-0017A477100C.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/B247E90D-D549-DF11-AE79-0017A4770424.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/A09DCD0D-D549-DF11-8D26-0017A477082C.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/82E06E82-DB49-DF11-97EE-0017A4771018.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/8265921B-D649-DF11-9E77-001E0B5FC57A.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/825F226F-D549-DF11-8DB4-0017A4770024.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/80D30FC4-D549-DF11-9F07-0017A4770404.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/520E3E1A-D649-DF11-B6E9-001E0B472CEE.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/1C59A56C-D549-DF11-8D44-00237DA1ED1C.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/104F9485-D949-DF11-871B-0017A477082C.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0053/02E1808F-D749-DF11-B166-0017A4770014.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0052/FA0ED97D-D449-DF11-AFF7-0017A477100C.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0052/F4E67FBB-D249-DF11-88C8-0017A477103C.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0052/EAE5CC15-CD49-DF11-821A-0017A4771020.root',
       #'/store/mc/Spring10/QCD_Pt170/GEN-SIM-RECO/START3X_V26_S09-v1/0052/E2EC6C2A-CE49-DF11-9DE8-00237DA10D14.root'
])


process.source = cms.Source("PoolSource",
    fileNames = readFiles,
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
    input = cms.untracked.int32(200)
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

configurePatTupleProduction(process, addGenInfo = isMC)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce Ntuple
#
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrigger_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigCaloTau_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGlobalVariables_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauFixedCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauShrinkingCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPS_cfi")

process.ntupleProducer = cms.EDAnalyzer("ObjValEDNtupleProducer",
                                        
    ntupleName = cms.string("tauIdEffNtuple"),
                                        
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc

        # variables indicating decision of HLT trigger paths
        trigger = process.trigger_template,                                    

        # variables specifying x,y,z coordinates of primary event vertex
        vertex = process.vertex_template,

        # test global variables
        pfTausShrinkingCone_part03 = process.extraVarforTauCands_template,
        jets = process.jets_template,
        met = process.met_template,


        # variables specific to CaloTaus                                            
        caloTaus_part01 = process.caloTaus_recInfo,

        # variables specific to fixed cone PFTaus                                            
        pfTausFixedCone_part01 = process.pfTausFixedCone_recInfo,

        # variables specific to shrinking cone PFTaus                                            
        pfTausShrinkingCone_part01 = process.pfTausShrinkingCone_recInfo,

        # variables specific to PFTaus reconstructed by hadron + strips (HPS) algorithm                                           
        pfTausHPS_part01 = process.pfTausHPS_recInfo,
    )
)

if isMC:
    setattr(process.ntupleProducer.sources, "caloTaus_part02", process.caloTaus_genInfo)
    setattr(process.ntupleProducer.sources, "pfTausFixedCone_part02", process.pfTausFixedCone_genInfo)
    setattr(process.ntupleProducer.sources, "pfTausShrinkingCone_part02", process.pfTausShrinkingCone_genInfo)
    setattr(process.ntupleProducer.sources, "pfTausHPS_part02", process.pfTausHPS_genInfo)
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
            #"keep *"
        )                               
    ),
    process.qcdDiJetEventSelection, # comment-out to disable filtering of events put in Tau Ntuple                         
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string("tauIdEff_ntuple.root")      
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

# print-out all python configuration parameter information
#print process.dumpPython()
