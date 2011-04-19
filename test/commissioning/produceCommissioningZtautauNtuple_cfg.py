import FWCore.ParameterSet.Config as cms

process = cms.Process("prodCommissioningZtautauNtuple")

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
process.load('Configuration/StandardSequences/Reconstruction_cff')
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

# print gen. particle information
process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(100)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

isMC = True # use for MC
##isMC = False # use for Data
addMCEmbeddingFlags = True # use for MC samples generated using MCEmbedding technique
##addMCEmbeddingFlags = False # use for MC not generated using MCEmbedding technique and for Data
##HLTprocessName = "HLT" # use for 2011 Data
HLTprocessName = "REDIGI311X" # use for Spring'11 reprocessed MC
pfCandidateCollection = "particleFlow" # pile-up removal disabled
##pfCandidateCollection = "pfNoPileUp" # pile-up removal enabled
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__HLTprocessName = #HLTprocessName#
#__pfCandidateCollection = #pfCandidateCollection#
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define GlobalTag to be used for event reconstruction
# (only relevant for HPS tau reconstruction algorithm)
if isMC:
    process.GlobalTag.globaltag = cms.string('START311_V2::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_311_V2::All')
#--------------------------------------------------------------------------------    

#--------------------------------------------------------------------------------
# define skimming criteria
# (in order to be able to produce Tau Ntuple directly from unskimmed Monte Carlo/datasets;
#  HLT single jet trigger passed && either two CaloJets or two PFJets of Pt > 10 GeV within |eta| < 2.5)
process.load('TauAnalysis.TauIdEfficiency.filterQCDdiJet_cfi')
process.hltSingleJet.selector.triggerPaths = cms.vstring(
    '*'
)
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
from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder import buildQCDdiJetTauSequence

process.load("PhysicsTools.PatAlgos.cleaningLayer1.tauCleaner_cfi")

# remove jets outside kinematic range Pt > 10 GeV && |eta| < 2.5 from Tau Ntuple
# (in order to speed-up plotting macros)
patCaloTauCleanerPrototype = process.cleanPatTaus.clone(
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(),
    finalCut = cms.string(
        'caloTauTagInfoRef().jetRef().pt() > 8.0 & abs(caloTauTagInfoRef().jetRef().eta()) < 2.5'
    )
)

patPFTauCleanerPrototype = process.cleanPatTaus.clone(
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(),
    finalCut = cms.string(
        'pfJetRef().pt() > 8.0 & abs(pfJetRef().eta()) < 2.5'
    )
)

retVal = configurePatTupleProduction(
    process, patSequenceBuilder = buildQCDdiJetTauSequence,
    patPFTauCleanerPrototype = patPFTauCleanerPrototype,
    patCaloTauCleanerPrototype = patCaloTauCleanerPrototype,
    hltProcess = HLTprocessName,
    isMC = isMC
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
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGlobalVariables_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrackVariables_cfi")
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

        # reweighting factors to be applied to Monte Carlo simulated events
        # in order to match vertex multiplicity distribution in Data                                             
        vertexMultReweight = process.vertexMultReweight_template, 

        # global variables describing the underlying event/
        # amount of hadronic activity                                            
        jets = process.jets_template,
        caloMet_rec = process.caloMet_template,
        pfMet_rec = process.pfMet_template,

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
    # add in information about generator level visible taus and all generator level jets
    setattr(process.ntupleProducer.sources, "tauGenJets", process.tauGenJets_genInfo)
    setattr(process.ntupleProducer.sources, "genJets", process.genJets_genInfo)
    setattr(process.ntupleProducer.sources, "genPhaseSpaceEventInfo", process.genPhaseSpaceEventInfo_template)
    setattr(process.ntupleProducer.sources, "genPileUpEventInfo", process.genPileUpEventInfo_template)

if addMCEmbeddingFlags:
    process.caloTaus_mcEmbeddingInfo.src = cms.InputTag(retVal["caloTauCollection"])
    setattr(process.ntupleProducer.sources, "caloTaus_gen", process.caloTaus_mcEmbeddingInfo)
    process.pfTausFixedCone_mcEmbeddingInfo.src = cms.InputTag(retVal["pfTauCollectionFixedCone"])
    setattr(process.ntupleProducer.sources, "pfTausFixedCone_gen", process.pfTausFixedCone_mcEmbeddingInfo)
    process.pfTausShrinkingCone_mcEmbeddingInfo.src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])
    setattr(process.ntupleProducer.sources, "pfTausShrinkingCone_gen", process.pfTausShrinkingCone_mcEmbeddingInfo)
    process.pfTausHPS_mcEmbeddingInfo.src = cms.InputTag(retVal["pfTauCollectionHPS"])
    setattr(process.ntupleProducer.sources, "pfTausHPS_gen", process.pfTausHPS_mcEmbeddingInfo)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# update InputTags for HLT trigger result object
# in case running on reprocessed Monte Carlo samples
#
if HLTprocessName != "HLT":
    process.hltSingleJet.selector.src = cms.InputTag('TriggerResults::' + HLTprocessName)
    process.patCaloTausTriggerEvent.processName = cms.string(HLTprocessName)
    process.patPFTausTriggerEventFixedCone.processName = cms.string(HLTprocessName)
    process.patPFTausTriggerEventShrinkingCone.processName = cms.string(HLTprocessName)
    process.patPFTausTriggerEventHPS.processName = cms.string(HLTprocessName)    
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
    fileName = cms.untracked.string("tauIdEffEDNtuple_qcdDiJet.root")      
)
#--------------------------------------------------------------------------------

if addMCEmbeddingFlags:

    process.mcEmbeddingTagAndProbeFlags_template = cms.EDProducer(
        src = cms.InputTag(''),
        srcTags = cms.InputTag('goodMuonIsolationTagAndProbeProducer', 'p4Tags'),
        srcProbes = cms.InputTag('goodMuonIsolationTagAndProbeProducer', 'p4Probes'),
        dRmatch_ = cms.double(0.2)
    )

    process.mcEmbeddingTagAndProbeFlags_caloTau = process.mcEmbeddingTagAndProbeFlags_template.clone(
        src = cms.InputTag(retVal["caloTauCollection"])
    )
    process.mcEmbeddingTagAndProbeFlags_pfTauFixedCone = process.mcEmbeddingTagAndProbeFlags_template.clone(
        src = cms.InputTag(retVal["pfTauCollectionFixedCone"])
    )
    process.mcEmbeddingTagAndProbeFlags_pfTauShrinkingCone = process.mcEmbeddingTagAndProbeFlags_template.clone(
        src = cms.InputTag(retVal["pfTauCollectionShrinkingCone"])
    )
    process.mcEmbeddingTagAndProbeFlags_pfTauHPS = process.mcEmbeddingTagAndProbeFlags_template.clone(
        src = cms.InputTag(retVal["pfTauCollectionHPS"])
    )
    process.patTupleProductionSequence._seq = process.patTupleProductionSequence._seq
      * process.mcEmbeddingTagAndProbeFlags_caloTau
      * process.mcEmbeddingTagAndProbeFlags_pfTauFixedCone
      * process.mcEmbeddingTagAndProbeFlags_pfTauShrinkingCone
      * process.mcEmbeddingTagAndProbeFlags_pfTauHPS

process.p = cms.Path(
    process.prePatProductionSequence
   + process.patTupleProductionSequence
   #+ process.printEventContent
   #+ process.printGenParticleList
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


