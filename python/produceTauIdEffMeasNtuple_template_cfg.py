import FWCore.ParameterSet.Config as cms

process = cms.Process("prodTauIdEffMeasNtuple")

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
        'file:/data2/friis/CMSSW_4_2_X/skims/06-27-MatthewsZTTEvents/crab_0_110627_082505/ZTTCands_merged_v1.root'
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

isMC = True # use for MC
##isMC = False # use for Data
##HLTprocessName = "HLT" # use for 2011 Data
HLTprocessName = "HLT" # use for Summer'11 MC
pfCandidateCollection = "particleFlow" # pile-up removal disabled
##pfCandidateCollection = "pfNoPileUp" # pile-up removal enabled
applyZrecoilCorrection = False
#applyZrecoilCorrection = True
#---------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__HLTprocessName = #HLTprocessName#
#__pfCandidateCollection = #pfCandidateCollection#
#__applyZrecoilCorrection = #applyZrecoilCorrection#
#
#--------------------------------------------------------------------------------

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

#process.prePatProductionSequence.remove(process.tautagging)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for configurating PAT-tuple production
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProductionTauIdEffMeasSpecific import *

patTupleConfig = configurePatTupleProductionTauIdEffMeasSpecific(
    process, hltProcess = HLTprocessName, addGenInfo = isMC, applyZrecoilCorrection = applyZrecoilCorrection)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce Ntuple
#
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigRunLumiSectionEventNumber_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrigger_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigCaloTau_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauFixedCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauShrinkingCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPS_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPSpTaNC_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTauIdEffMeas_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGlobalVariables_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenJets_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPhaseSpaceEventInfo_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPileUpEventInfo_cfi")

process.ntupleProducer = cms.EDProducer("ObjValEDNtupleProducer",
                                        
    ntupleName = cms.string("tauIdEffNtuple"),
                                        
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc

        # run number, luminosity section, event number variables
        runLumiSectionEventNumber = process.runLumiSectionEventNumber_template,

        # variables indicating decision of HLT trigger paths
        trigger = process.trigger_template,

        # variables specifying x,y,z coordinates of primary event vertices
        vertex = process.vertex_template,

        # number of reconstructed primary event vertices
        # with sum(trackPt) exceeding different thresholds
        vertexMultiplicity = process.vertexMultiplicity_template,

        # number of reconstructed global muons in the event
        numGlobalMuons = process.numGlobalMuons_template,
        numStandAloneMuons = process.numStandAloneMuons_template,                                    

        # global variables describing the underlying event/
        # amount of hadronic activity                                                                                            
        caloMet_rec = process.caloMet_template,
        pfMet_rec = process.pfMet_template,

        # variables specific to tau id. efficiency measurement
        tauIdEffMeas01 = process.tauIdEffMeas_template01.clone(
            src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative')
        ),
        tauIdEffMeas05 = process.tauIdEffMeas_template05.clone(
            src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative')
        )                                                                               
    )
)

patTauTemplates = {
    "caloTau"            : process.caloTaus_recInfo,
    "pfTauFixedCone"     : process.pfTausFixedCone_recInfo,
    "pfTauShrinkingCone" : process.pfTausShrinkingCone_recInfo,
    "pfTauHPS"           : process.pfTausHPS_recInfo,
    "pfTauHPSpTaNC"      : process.pfTausHPSpTaNC_recInfo
}

muTauTemplates = {
    "caloTau"            : [ process.tauIdEffMeas_template03caloTau, process.tauIdEffMeas_template04caloTau ], 
    "pfTauFixedCone"     : [ process.tauIdEffMeas_template03pfTau,   process.tauIdEffMeas_template04pfTau   ],
    "pfTauShrinkingCone" : [ process.tauIdEffMeas_template03pfTau,   process.tauIdEffMeas_template04pfTau   ],
    "pfTauHPS"           : [ process.tauIdEffMeas_template03pfTau,   process.tauIdEffMeas_template04pfTau   ],
    "pfTauHPSpTaNC"      : [ process.tauIdEffMeas_template03pfTau,   process.tauIdEffMeas_template04pfTau   ]
}

# add branches for collections of Tau-jet candidates and Muon + Tau-Jet candidate pairs
# (including collecions shifted Up/Down in energy/momentum and with Z-recoil corrections applied)
for algorithm in patTupleConfig["algorithms"]:

    for patJetCollection in patTupleConfig[algorithm]["patJetCollections"]:
        attrJetsName = "tauIdEffMeasJets_%s" % patJetCollection
        attrJets = process.jets_template.clone(
            src = cms.InputTag(patJetCollection),
            srcNotToBeFiltered = cms.VInputTag()
        )
        setattr(process.ntupleProducer.sources, attrJetsName, attrJets)

    for patTauCollection in patTupleConfig[algorithm]["patTauCollections"]:                
        attr02Name = "tauIdEffMeas02_%s" % patTauCollection
        attr02 = patTauTemplates[algorithm].clone(
            src = cms.InputTag(patTauCollection),
            vector = cms.bool(True),
            indices = cms.vuint32([0])
        )
        setattr(process.ntupleProducer.sources, attr02Name, attr02)

    for muTauPairCollection in patTupleConfig[algorithm]["muTauPairCollections"]:
        attr03Name = "tauIdEffMeas03_%s" % muTauPairCollection
        attr03 = muTauTemplates[algorithm][0].clone(
            src = cms.InputTag(muTauPairCollection)
        )
        setattr(process.ntupleProducer.sources, attr03Name, attr03)
        attr04Name = "tauIdEffMeas04_%s" % muTauPairCollection
        attr04 = muTauTemplates[algorithm][1].clone(
            src = cms.InputTag(muTauPairCollection)
        )
        setattr(process.ntupleProducer.sources, attr04Name, attr04)

if isMC:
    # add in information about generator level visible taus and all generator level jets
    setattr(process.ntupleProducer.sources, "tauGenJets", process.tauGenJets_genInfo)
    setattr(process.ntupleProducer.sources, "genJets", process.genJets_genInfo)
    setattr(process.ntupleProducer.sources, "genPhaseSpaceEventInfo", process.genPhaseSpaceEventInfo_template)
    setattr(process.ntupleProducer.sources, "genPileUpEventInfo", process.genPileUpEventInfo_template)
    # add reweighting factors to be applied to Monte Carlo simulated events
    # in order to match vertex multiplicity distribution in Data                                             
    setattr(process.ntupleProducer.sources, "vertexMultReweight", process.vertexMultReweight_template)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#                                   Put the correct jet energy correction
#------------------------------------------------------------------------------------------------------------------------
process.load('CondCore.DBCommon.CondDBSetup_cfi')
process.jec = cms.ESSource("PoolDBESSource",
    process.CondDBSetup,
    ## DBParameters = cms.PSet(
    ##     messageLevel = cms.untracked.int32(0)
    ##     ),
    ## timetype = cms.string('runnumber'),
    toGet = cms.VPSet(
         cms.PSet(record = cms.string("JetCorrectionsRecord"),
                  tag = cms.string("JetCorrectorParametersCollection_Jec10V3_AK5Calo"),#JetCorrectorParametersCollection_Jec11_V1_AK5Calo
                  label=cms.untracked.string("AK5Calo")),
         cms.PSet(record = cms.string("JetCorrectionsRecord"),
                  tag = cms.string("JetCorrectorParametersCollection_Jec10V3_AK5PF"),
                  label=cms.untracked.string("AK5PF")),                                   
         cms.PSet(record = cms.string("JetCorrectionsRecord"),
                  tag = cms.string("JetCorrectorParametersCollection_Jec10V3_AK5PFchs"),
                  label=cms.untracked.string("AK5PF"))
    ),
    ## here you add as many jet types as you need (AK5Calo, AK5JPT, AK7PF, AK7Calo, KT4PF, KT4Calo, KT6PF, KT6Calo)
    connect = cms.string('sqlite_fip:TauAnalysis/Configuration/data/Jec10V3.db')
    #connect = cms.string("frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS")
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')
#-------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
#
# Save ntuple
#
process.ntupleOutputModule = cms.OutputModule("PoolOutputModule",
    cms.PSet(
        outputCommands = cms.untracked.vstring(
            "drop *",
            "keep *_*ntupleProducer*_*_*",
            "keep edmMergeableCounter_*_*_*"
        )                               
    ),
    process.tauIdEffSampleEventSelection,
    fileName = cms.untracked.string("tauIdEffMeasEDNtuple.root")
)                                              
#--------------------------------------------------------------------------------

process.p = cms.Path(
    process.prePatProductionSequence
   + process.patDefaultSequence
   + process.producePatTupleTauIdEffMeasSpecific
   #+ process.printEventContent
   #+ process.printGenParticleList
   + process.ntupleProducer
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath(process.ntupleOutputModule)

#--------------------------------------------------------------------------------
#                   Modify the content of the extractors
#--------------------------------------------------------------------------------
import TauAnalysis.Configuration.pathModifiers as pathModifiers
pathModifiers.ExtractorAddColumn(process.ntupleProducer.sources,
                                 'PATTauVectorValExtractor', 'productionVertexZ', cms.string("vertex().z()"), True)
pathModifiers.ExtractorAddColumn(process.ntupleProducer.sources,
                                 'PATMuTauPairValExtractor', 'tauVertexZ', cms.string("leg2().vertex().z()"), True)
pathModifiers.ExtractorAddColumn(process.ntupleProducer.sources,
                                 'PATMuTauPairValExtractor', 'muVertexZ', cms.string("leg1().vertex().z()"), True)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# CV: do **not** apply HLT trigger conditions to CMSSW_4_1_x MC,
#     weight simulated events by trigger efficiencies measured in Data instead
#
if isMC:
    process.commonSkimSequence.remove(process.hltMu)
#--------------------------------------------------------------------------------    

# define order in which different paths are run
process.schedule = cms.Schedule(
    process.p,
    process.muonCaloTauSkimPath,
    process.muonPFTauFixedConeSkimPath,
    process.muonPFTauShrinkingConeSkimPath,
    process.muonPFTauHPSskimPath,
    process.muonPFTauHPSpTaNCskimPath,
    process.o
)

# print-out all python configuration parameter information
#print process.dumpPython()
processDumpFile = open('OldNtupleDump.py' , 'w')
print >> processDumpFile, process.dumpPython()
