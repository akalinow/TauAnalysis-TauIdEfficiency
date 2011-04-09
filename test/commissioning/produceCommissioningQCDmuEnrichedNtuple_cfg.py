import FWCore.ParameterSet.Config as cms

process = cms.Process("prodCommissioningQDCmuEnrichedNtuple")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.MessageLogger.suppressInfo = cms.untracked.vstring()
process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#--------------------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/data2/veelken/CMSSW_4_1_x/skims/ZtoMuTau/DYtautau_spring11_powhegZ2_1_1_XvY.root'
        #'file:/data2/veelken/CMSSW_4_1_x/skims/ZtoMuTau/data2/veelken/CMSSW_4_1_x/skims/ZtoMuTau/data2011A_tauPlusX_AOD_1_1_MV9.root'
    ),
    skipEvents = cms.untracked.uint32(0)            
)

# print event content 
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

isMC = True # use for MC
##isMC = False # use for Data
##HLTprocessName = "HLT" # use for 2011 Data
HLTprocessName = "REDIGI311X" # use for Spring'11 reprocessed MC
pfCandidateCollection = "particleFlow" # pile-up removal disabled
##pfCandidateCollection = "pfNoPileUp" # pile-up removal enabled
applyEventSelection = True 
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__HLTprocessName = #HLTprocessName#
#__pfCandidateCollection = #pfCandidateCollection#
#__applyEventSelection = #applyEventSelection#
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define GlobalTag to be used for event reconstruction
if isMC:
    process.GlobalTag.globaltag = cms.string('START311_V2::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_311_V2::All')
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
# produce collections of objects needed as input for PAT-tuple production
# (e.g. rerun reco::Tau identification algorithms with latest tags)
#
from TauAnalysis.TauIdEfficiency.tools.configurePrePatProduction import configurePrePatProduction

configurePrePatProduction(process, pfCandidateCollection = pfCandidateCollection, addGenInfo = isMC)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce PAT objects
#
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProduction import configurePatTupleProduction
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildGenericTauSequence

# add muon isolation variables
process.load("PhysicsTools.PFCandProducer.Isolation.pfMuonIsolation_cff")
from PhysicsTools.PFCandProducer.Isolation.tools_cfi import *
process.pfmuIsoDepositPFCandidates = isoDepositReplace("muons", "particleFlow")
process.prePatProductionSequence._seq = process.prePatProductionSequence._seq * process.pfmuIsoDepositPFCandidates

process.load("PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi")
process.patMuons.userIsolation.pfAllParticles = cms.PSet( 
    src = cms.InputTag("pfmuIsoDepositPFCandidates"),
    deltaR = cms.double(0.4)
)

# "clean" CaloTau/PFTau collections
# (i.e. remove CaloTaus/PFTaus overlapping with muons)
process.load("PhysicsTools.PatAlgos.cleaningLayer1.tauCleaner_cfi")

# remove jets outside kinematic range Pt > 10 GeV && |eta| < 2.5 from Tau Ntuple
# (in order to speed-up plotting macros)
patCaloTauCleanerPrototype = process.cleanPatTaus.clone(
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
    finalCut = cms.string(
        'caloTauTagInfoRef().jetRef().pt() > 8.0 & abs(caloTauTagInfoRef().jetRef().eta()) < 2.5'
    )
)

patPFTauCleanerPrototype = patCaloTauCleanerPrototype.clone(
    finalCut = cms.string(
        'pfJetRef().pt() > 8.0 & abs(pfJetRef().eta()) < 2.5'
    )
)

retVal = configurePatTupleProduction(
    process,
    patPFTauCleanerPrototype = patPFTauCleanerPrototype,
    patCaloTauCleanerPrototype = patCaloTauCleanerPrototype,
    hltProcess = HLTprocessName,
    addGenInfo = isMC
)
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
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPSpTaNC_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigDiTauKinematics_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenJets_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPhaseSpaceEventInfo_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPileUpEventInfo_cfi")

process.ntupleProducer = cms.EDProducer("ObjValEDNtupleProducer",
                                        
    ntupleName = cms.string("tauIdEffNtuple"),
                                        
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc

        # variables indicating decision of HLT trigger paths
        trigger = process.trigger_template,

        # variables specifying x,y,z coordinates of primary event vertices
        vertex = process.vertex_template,

        # number of reconstructed primary event vertices
        # with sum(trackPt) exceeding different thresholds
        vertexMultiplicity = process.vertexMultiplicity_template,

        # variables specific to Muons
        muons_rec = process.muons_recInfo,              

        # variables specific to CaloTaus                                            
        caloTaus_rec = process.caloTaus_recInfo.clone(
            src = cms.InputTag(retVal["caloTauCollection"])                       
        ),
        caloTaus_recJetId = process.caloTaus_recJetIdInfo.clone(
            src = cms.InputTag(retVal["caloTauCollection"])                       
        ),                                    
        muCaloTauPairs_rec = process.diTaus_recInfo.clone(
            src = cms.InputTag(retVal["muonCaloTauCollection"])                       
        ),                                             

        # variables specific to fixed cone PFTaus                                            
        pfTausFixedCone_rec = process.pfTausFixedCone_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionFixedCone"])                       
        ),
        pfTausFixedCone_recJetId = process.pfTausFixedCone_recJetIdInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionFixedCone"])                       
        ),                                      
        muPFTauPairsFixedCone_rec = process.diTaus_recInfo.clone(
            src = cms.InputTag(retVal["muonPFTauCollectionFixedCone"])                       
        ),                                           

        # variables specific to shrinking cone PFTaus                                            
        pfTausShrinkingCone_rec = process.pfTausShrinkingCone_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])                       
        ),                                            
        pfTausShrinkingCone_recJetId = process.pfTausShrinkingCone_recJetIdInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])
        ),                                   
        muPFTauPairsShrinkingCone_rec = process.diTaus_recInfo.clone(
            src = cms.InputTag(retVal["muonPFTauCollectionShrinkingCone"])                       
        ),                                            

        # variables specific to PFTaus reconstructed by hadron + strips (HPS) algorithm                                           
        pfTausHPS_rec = process.pfTausHPS_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionHPS"])                       
        ),
        pfTausHPS_recJetId = process.pfTausHPS_recJetIdInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionHPS"])
        ),                                       
        muPFTauPairsHPS_rec = process.diTaus_recInfo.clone(
            src = cms.InputTag(retVal["muonPFTauCollectionHPS"])                       
        ),                                            

        # variables specific to PFTaus reconstructed by HPS + TaNC combined tau id. algorithm
        pfTausHPSpTaNC_rec = process.pfTausHPSpTaNC_recInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionHPSpTaNC"])
        ),
        pfTausHPSpTaNC_recJetId = process.pfTausHPSpTaNC_recJetIdInfo.clone(
            src = cms.InputTag(retVal["pfTauCollectionHPSpTaNC"])
        ),                                       
        muPFTauPairsHPSpTaNC_rec = process.diTaus_recInfo.clone(
            src = cms.InputTag(retVal["muonPFTauCollectionHPSpTaNC"])                       
        )                                           
    )
)

if isMC:
    setattr(process.ntupleProducer.sources, "muons_gen", process.muons_genInfo)
    process.caloTaus_genInfo.src = cms.InputTag(retVal["caloTauCollection"])
    setattr(process.ntupleProducer.sources, "caloTaus_gen", process.caloTaus_genInfo)
    process.pfTausFixedCone_genInfo.src = cms.InputTag(retVal["pfTauCollectionFixedCone"])
    setattr(process.ntupleProducer.sources, "pfTausFixedCone_gen", process.pfTausFixedCone_genInfo)
    process.pfTausShrinkingCone_genInfo.src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])
    setattr(process.ntupleProducer.sources, "pfTausShrinkingCone_gen", process.pfTausShrinkingCone_genInfo)
    process.pfTausHPS_genInfo.src = cms.InputTag(retVal["pfTauCollectionHPS"])
    setattr(process.ntupleProducer.sources, "pfTausHPS_gen", process.pfTausHPS_genInfo)
    process.pfTausHPSpTaNC_genInfo.src = cms.InputTag(retVal["pfTauCollectionHPSpTaNC"])
    setattr(process.ntupleProducer.sources, "pfTausHPSpTaNC_gen", process.pfTausHPSpTaNC_genInfo)    
    # add in information about generator level visible taus and all generator level jets
    setattr(process.ntupleProducer.sources, "tauGenJets", process.tauGenJets_genInfo)
    setattr(process.ntupleProducer.sources, "genJets", process.genJets_genInfo)
    setattr(process.ntupleProducer.sources, "genPhaseSpaceEventInfo", process.genPhaseSpaceEventInfo_template)
    setattr(process.ntupleProducer.sources, "genPileUpEventInfo", process.genPileUpEventInfo_template)    
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# update InputTags for HLT trigger result object
# in case running on reprocessed Monte Carlo samples
#
if HLTprocessName != "HLT":
    process.hltMu.selector.src = cms.InputTag('TriggerResults::' + HLTprocessName)
    process.patTrigger.processName = HLTprocessName
    process.patTriggerEvent.processName = HLTprocessName
    process.patCaloTausTriggerEvent.processName = cms.string(HLTprocessName)
    process.patPFTausTriggerEventFixedCone.processName = cms.string(HLTprocessName)
    process.patPFTausTriggerEventShrinkingCone.processName = cms.string(HLTprocessName)
    process.patPFTausTriggerEventHPS.processName = cms.string(HLTprocessName)
    process.patPFTausTriggerEventHPSpTaNC.processName = cms.string(HLTprocessName)
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
    fileName = cms.untracked.string("tauIdEffEDNtuple_qcdMuEnriched.root")      
)
#--------------------------------------------------------------------------------

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
if applyEventSelection:
    process.schedule = cms.Schedule(
        process.p,
        process.muonCaloTauSkimPath,
        process.muonPFTauSkimPath,
        process.o
    )
else:
    delattr(process.ntupleOutputModule, "SelectEvents")
    process.schedule = cms.Schedule(
        process.p,
        process.o
    )

# print-out all python configuration parameter information
#print process.dumpPython()


