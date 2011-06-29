import FWCore.ParameterSet.Config as cms

process = cms.Process("prodTauIdEffMeasNtuple")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.MessageLogger.suppressInfo = cms.untracked.vstring()
process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#--------------------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/m/mverzett/tagprobe/skims/TauIdEffMeas_2011May13/tauIdEffSample_Ztautau_powheg_2011May20_RECO_266_1_qvJ.root'
        #'rfio:/castor/cern.ch/user/m/mverzett/tagprobe/skims/TauIdEffMeas_2011May13/tauIdEffSample_DYtautauM10to20_powheg_2011May20_RECO_14_1_yrh.root'
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
    input = cms.untracked.int32(-1)
)

##isMC = True # use for MC
isMC = False # use for Data
#HLTprocessName = "HLT" # use for 2011 Data
HLTprocessName = "REDIGI311X" # use for Spring'11 reprocessed MC
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
    process.GlobalTag.globaltag = cms.string('START311_V2::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_P_V14::All')
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
    process, hltProcess = HLTprocessName, isMC = isMC, applyZrecoilCorrection = applyZrecoilCorrection)
process.patTrigger.addL1Algos = False
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
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
#                   Pileup reweight (from old NTuple)
#-------------------------------------------------------------------------------------------------------------------------
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.ntupleProducer = cms.EDProducer("ObjValEDNtupleProducer",
                                        ntupleName = cms.string("tauIdEffNtuple"),
                                        sources = cms.PSet(                                            
                                            )
                                        )
# in order to match vertex multiplicity distribution in Data                                             
setattr(process.ntupleProducer.sources, "vertexMultReweight", process.vertexMultReweight_template)
#-------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# Save ntuple
#
process.ntupleOutputModule = cms.OutputModule("PoolOutputModule",
                                              cms.PSet(
                                                  outputCommands = cms.untracked.vstring(
                                                      'drop *',
                                                      'keep edmMergeableCounter_*_*_*',
                                                      'keep *TriggerEvent_*_*_*',
                                                      'keep patMuons_selectedPatMuonsForTauIdEffTrkIPcumulative_*_*',
                                                      'keep *_patMuonsStandAlone_*_*',
                                                      'keep patTaus_selectedPatPFTausShrinkingConePFRelIsoCumulative_*_*',
                                                      'keep patTaus_selectedPatPFTausHPSpTaNCPFRelIsoCumulative_*_*',
                                                      'keep patTaus_selectedPatPFTausHPSPFRelIsoCumulative_*_*',
                                                      'keep patTaus_selectedPatPFTausHPSpTaNCPFRelIsoCumulative_*_*',
                                                      'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauShrinkingConePairsForTauIdEffCumulative_*_*',
                                                      'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative_*_*',
                                                      'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpairsForTauIdEffCumulative_*_*',
                                                      'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative_*_*',
                                                      'keep *Vertex*_*_*_*',
                                                      'keep *_ak5PFJets_*_*',
                                                      'keep *_ak5CaloJets_*_*',
                                                      'keep *PFCandidate*_*_*_*'
                                                      )               
                                                  ),
                                              process.tauIdEffSampleEventSelection,
                                              fileName = cms.untracked.string("tauIdEffMeasEDNtuple.root")      
                                              )
if isMC:
    process.ntupleOutputModule.outputCommands.extend(
      cms.untracked.vstring(
            "keep *_*_*vtxMultReweight*_*",
            'keep *GenJet*_*_*_*',
            'keep *recoGenMET_*_*_*',
            'keep *_genParticles_*_*',
            'keep *_generator_*_*',
            'keep *_*GenJets_*_*',
            'keep *_addPileupInfo_*_*'
            )
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

import TauAnalysis.Configuration.pathModifiers as pathModifiers
pathModifiers.PathModifier(process.p, 'doSVreco', cms.bool(False), True)
pathModifiers.PathModifier(process.muonPFTauFixedConeSkimPath, 'doSVreco', cms.bool(False),True)
pathModifiers.PathModifier(process.muonPFTauShrinkingConeSkimPath, 'doSVreco', cms.bool(False),True)
pathModifiers.PathModifier(process.muonPFTauHPSskimPath, 'doSVreco', cms.bool(False),True)
pathModifiers.PathModifier(process.muonPFTauHPSpTaNCskimPath, 'doSVreco', cms.bool(False),True)
## process.allMuTauPairs.doSVreco = cms.bool(False)

# define order in which different paths are run
process.schedule = cms.Schedule(
    process.p,
    process.muonPFTauFixedConeSkimPath,
    process.muonPFTauShrinkingConeSkimPath,
    process.muonPFTauHPSskimPath,
    process.muonPFTauHPSpTaNCskimPath,
    process.o
)
